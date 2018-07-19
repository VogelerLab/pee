######################################################################################################### 
##                                                                                                    #\\
##                                         LANDTRENDR LIBRARY                                         #\\
##                                                                                                    #\\
#########################################################################################################


# date: 2018-06-11
# author: Zhiqiang Yang  | zhiqiang.yang@oregonstate.edu
#         Justin Braaten | jstnbraaten@gmail.com
#         Robert Kennedy | rkennedy@coas.oregonstate.edu
# website: https:#github.com/eMapR/LT-GEE

import ee


#########################################################################################################
###### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS ##### 
#########################################################################################################


#------ L8 to L7 HARMONIZATION FUNCTION -----
# slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K.,
# Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization
# of Landsat-7 to Landsat-8 reflective wavelength and normalized
# difference vegetation index continuity, Remote Sensing of Environment,
# 185, 57-70.(http:#dx.doi.org/10.1016/j.rse.2015.12.024) Table 2 -
# reduced major axis (RMA) regression coefficients


def harmonizationRoy(oli):
    # create an image of slopes per band for L8 TO L7 regression line - David
    # Roy
    slopes = ee.Image.constant(
        [0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949])
    # create an image of y-intercepts per band for L8 TO L7 regression line -
    # David Roy
    itcp = ee.Image.constant(
        [-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029])
    y = (oli.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'], ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])  # select OLI bands 2-7 and rename them to match L7 band names
         # ...resample the L8 bands using bicubic
         .resample('bicubic')
         # ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
         .subtract(itcp.multiply(10000)).divide(slopes)
         .set('system:time_start', oli.get('system:time_start')))                      # ...set the output system:time_start metadata to the input image time_start otherwise it is null
    # return the image as short to match the type of the other data
    return y.toShort()


#------ FILTER A COLLECTION FUNCTION -----
def filterCollection(year, startDay, endDay, sensor, aoi):
    return (ee.ImageCollection('LANDSAT/' + sensor + '/C01/T1_SR')
            .filterBounds(aoi)
            .filterDate(year + '-' + startDay, year + '-' + endDay))


#------ BUILD A COLLECTION FOR A GIVEN SENSOR AND YEAR -----
def buildSensorYearCollection(year, startDay, endDay, sensor, aoi):
    startMonth = int(startDay[:2])
    endMonth = int(endDay[:2])
    if startMonth > endMonth:
        oldYear = str(int(year) - 1)
        newYear = year
        oldYearStartDay = startDay
        oldYearEndDay = '12-31'
        newYearStartDay = '01-01'
        newYearEndDay = endDay

        oldYearCollection = filterCollection(
            oldYear, oldYearStartDay, oldYearEndDay, sensor, aoi)
        newYearCollection = filterCollection(
            newYear, newYearStartDay, newYearEndDay, sensor, aoi)

        srCollection = ee.ImageCollection(
            oldYearCollection.merge(newYearCollection))

    else:
        srCollection = filterCollection(year, startDay, endDay, sensor, aoi)

    return srCollection




#------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
def getSRcollection(year, startDay, endDay, sensor, aoi):
    # get a landsat collection for given year, day range, and sensor
    srCollection = buildSensorYearCollection(
        year, startDay, endDay, sensor, aoi)

    # srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR') # get surface reflectance images
    #                     .filterBounds(aoi)                                  # ...filter them by intersection with AOI
    #                     .filterDate(year+'-'+startDay, year+'-'+endDay)    # ...filter them by year and day range

    # apply the harmonization function to LC08 (if LC08), subset bands,
    # unmask, and resample
    def l7_l8_harmonize(img):
        dat = (ee.Image(
            ee.Algorithms.If(
                # condition - if image is OLI
                img.get("SATELLITE").compareTo('LANDSAT_8').eq(0),
                # true - then apply the L8 TO L7 alignment function after
                # unmasking pixels that were previosuly masked (why/when are
                # pixels masked)
                harmonizationRoy(img.unmask()),
                # false - else select out the reflectance bands from the
                # non-OLI image
                img.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])
                # ...unmask any previously masked pixels
                .unmask()
                # ...resample by bicubic
                .resample('bicubic')
                .set('system:time_start', img.get('system:time_start'))))         # ...set the output system:time_start metadata to the input image time_start otherwise it is null
               )
        return dat
     
    # apply the mask to the image and return it
    def mask_img(img):
        # make a cloud, cloud shadow, and snow mask from fmask band
        # select out the fmask band
        qa = img.select('pixel_qa')
        mask = qa.bitwiseAnd(8).eq(0).And(qa.bitwiseAnd(16).eq(0)).And(qa.bitwiseAnd(32).eq(0))

        # apply the mask - 0's in mask will be excluded from computation and
        # set to opacity=0 in display
        return img.mask(mask)

    srCollection = srCollection.map(l7_l8_harmonize)
    srCollection = srCollection.map(mask_img)

    return srCollection  # return the prepared collection


#------ FUNCTION TO COMBINE LT05, LE07, & LC08 COLLECTIONS -----
def getCombinedSRcollection(year, startDay, endDay, aoi):
    # get TM collection for a given year, date range, and area
    lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi)
    # get ETM+ collection for a given year, date range, and area
    le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi)
    # get OLI collection for a given year, date range, and area
    lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi)
    # merge the individual sensor collections into one imageCollection object
    mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8))
    return mergedCollection                                              # return the Imagecollection


#------ FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE PER YEAR BY MEDOID -----


# make a medoid composite with equal weight among indices
def medoidMosaic(inCollection, dummyCollection):
    """
    LT expects only a single image per year in a time series, there are lost of ways to
    do best available pixel compositing - we have found that a mediod composite requires little logic
    is robust, and fast

    Medoids are representative objects of a data set or a cluster with a data set whose average 
    dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to 
    means or centroids, but medoids are always members of the data set.
    """

    # fill in missing years with the dummy collection
    # get the number of images
    imageCount = inCollection.toList(1).length()
    # if the number of images in this year is 0, then use the dummy
    # collection, otherwise use the SR collection
    finalCollection = ee.ImageCollection(ee.Algorithms.If(
        imageCount.gt(0), inCollection, dummyCollection))

    # calculate median across images in collection per band
    # calculate the median of the annual image collection - returns a
    # single 6 band image - the collection median per band
    median = finalCollection.median()

    # calculate the different between the median and the observation per image
    # per band
    def get_dif_from_med(img):
        # get the difference between each image/band and the corresponding band
        # median and take to power of 2 to make negatives positive and make
        # greater differences weight more
        diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2))
        # per image in collection, sum the powered difference across the bands
        # - set this as the first band add the SR bands to it - now a 7 band
        # image collection
        return diff.reduce('sum').addBands(img)

    difFromMedian = finalCollection.map(get_dif_from_med)

    # get the medoid by selecting the image pixel with the smallest difference
    # between median and observation per band
    # find the powered difference that is the least - what image object is the
    # closest to the median of teh collection - and then subset the SR bands
    # and name them - leave behind the powered difference band
    return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1, 2, 3, 4, 5, 6], ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])


#------ FUNCTION TO APPLY MEDOID COMPOSITING FUNCTION TO A COLLECTION -------------------------------------------
# create a temp iable to hold the upcoming annual mosiac
def buildMosaic(year, startDay, endDay, aoi, dummyCollection):
    # add the year to each medoid image - the data is hard-coded Aug 1st
    from datetime import date
    from time import mktime
    millis = mktime(date(year, 8, 1).timetuple()) * 1000.
    
    # get the SR collection
    collection = getCombinedSRcollection(year, startDay, endDay, aoi)  
    
    # apply the medoidMosaic function to reduce the collection to single image per year by medoid
    img = medoidMosaic(collection, dummyCollection).set('system:time_start', millis)  
    
    return ee.Image(img)


#------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
def buildSRcollection(startYear, endYear, startDay, endDay, aoi):
    # add the year to each image - the data is hard-coded Aug 1st
    from datetime import date
    from time import mktime
    
    # make an image collection from an image with 6 bands all set to 0 and
    # then make them masked values
    dummyCollection = ee.ImageCollection([ee.Image([0, 0, 0, 0, 0, 0]).mask(ee.Image(0))])
    # create empty array to fill
    imgs = []
    for i in range(startYear, endYear+1):                                           # for each year from hard defined start to end build medoid composite and then add to empty img array
        millis = mktime(date(i, 8, 1).timetuple()) * 1000.
        tmp = buildMosaic(i, startDay, endDay, aoi, dummyCollection)                    # build the medoid mosaic for a given year
        imgs = imgs.concat(tmp.set('system:time_start', millis))       # concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1st
  
    return ee.ImageCollection(imgs)                                                       # return the array img array as an image collection


#########################################################################################################
###### UNPACKING LT-GEE OUTPUT STRUCTURE FUNCTIONS ##### 
#########################################################################################################

# ----- FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS -----
def getLTvertStack(LTresult, run_params):
    # make empty array to hold another array whose length will y depending on
    # maxSegments parameter
    emptyArray = [0] * (run_params["maxSegments"] + 1)
    # make empty array to hold band names whose length will y depending on
    # maxSegments parameter
    vertLabels = ["vert_" + i for i in range(run_params["maxSegments"] + 1)]

    # make an image to fill holes in result 'LandTrendr' array where vertices
    # found is not equal to maxSegments parameter plus 1
    zeros = ee.Image(ee.Array([emptyArray, emptyArray, emptyArray]))

    # labels for 2 dimensions of the array that will be cast to each other
    # in the final step of creating the vertice output
    lbls = [['yrs_', 'src_', 'fit_'], vertLabels, ]

    # slices out the 4th row of a 4 row x N col (N = number of years in
    # annual stack) matrix, which identifies vertices - contains only 0s
    # and 1s, where 1 is a vertex (referring to spectral-temporal
    # segmentation) year and 0 is not
    vmask = LTresult.arraySlice(0, 3, 4)

    ltVertStack = (LTresult.arrayMask(vmask)       # uses the sliced out isVert row as a mask to only include vertice in this data - after this a pixel will only contain as many "bands" are there are vertices for that pixel - min of 2 to max of 7.
                   # ...from the vertOnly data subset slice out the vert year row, raw spectral row, and fitted spectral row
                   .arraySlice(0, 0, 3)
                   # ...adds the 3 row x 7 col 'zeros' matrix as a band to the vertOnly array - this is an intermediate step to the goal of filling in the vertOnly data so that there are 7 vertice slots represented in the data - right now there is a mix of lengths from 2 to 7
                   .addBands(zeros)
                   # ...concatenates the 3 row x 7 col 'zeros' matrix band to the vertOnly data so that there are at least 7 vertice slots represented - in most cases there are now > 7 slots filled but those will be truncated in the next step
                   .toArray(1)
                   # ...before this line runs the array has 3 rows and between 9 and 14 cols depending on how many vertices were found during segmentation for a given pixel. this step truncates the cols at 7 (the max verts allowed) so we are left with a 3 row X 7 col array
                   .arraySlice(1, 0, run_params["maxSegments"] + 1)
                   .arrayFlatten(lbls, ''))      # ...this takes the 2-d array and makes it 1-d by stacking the unique sets of rows and cols into bands. there will be 7 bands (vertices) for vertYear, followed by 7 bands (vertices) for rawVert, followed by 7 bands (vertices) for fittedVert, according to the 'lbls' list

    return ltVertStack                               # return the stack
  




# #######################################################################################
# ###### INDEX CALCULATION FUNCTIONS ####################################################
# #######################################################################################

# TODO: convert to python

## TASSELLED CAP
# tcTransform = function(img){ 
#   b = ee.Image(img).select(["B1", "B2", "B3", "B4", "B5", "B7"]) # select the image bands
#   brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]) # set brt coeffs - make an image object from a list of values - each of list element represents a band
#   grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]) # set grn coeffs - make an image object from a list of values - each of list element represents a band
#   wet_coeffs = ee.Image.constant([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]) # set wet coeffs - make an image object from a list of values - each of list element represents a band
#  
#   sum = ee.Reducer.sum() # create a sum reducer to be applyed in the next steps of summing the TC-coef-weighted bands
#   brightness = b.multiply(brt_coeffs).reduce(sum) # multiply the image bands by the brt coef and then sum the bands
#   greenness = b.multiply(grn_coeffs).reduce(sum) # multiply the image bands by the grn coef and then sum the bands
#   wetness = b.multiply(wet_coeffs).reduce(sum) # multiply the image bands by the wet coef and then sum the bands
#   angle = (greenness.divide(brightness)).atan().multiply(180/Math.PI).multiply(100)
#   tc = brightness.addBands(greenness)
#                     .addBands(wetness)
#                     .addBands(angle)
#                     .select([0,1,2,3], ['TCB','TCG','TCW','TCA']) #stack TCG and TCW behind TCB with .addBands, use select() to name the bands
#                     .set('system:time_start', img.get('system:time_start'))
#  return tc
#}
#
## NBR
# nbrTransform = function(img) {
#     nbr = img.normalizedDifference(['B4', 'B7']) # calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
#                 .multiply(1000) # scale results by 1000
#                 .select([0], ['NBR']) # name the band
#                 .set('system:time_start', img.get('system:time_start'))
#    return nbr
#}
#
## NDVI
# ndviTransform = function(img){ 
#   ndvi = img.normalizedDifference(['B4', 'B3']) # calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
#                .multiply(1000) # scale results by 1000
#                .select([0], ['NDVI']) # name the band
#                .set('system:time_start', img.get('system:time_start'))
#  return ndvi
#}
#                
## NDSI
# ndsiTransform = function(img){ 
#   ndsi = img.normalizedDifference(['B2', 'B5']) # calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
#                .multiply(1000) # scale results by 1000
#                .select([0], ['NDSI']) # name the band
#                .set('system:time_start', img.get('system:time_start'))
#  return ndsi
#}
#
## NDMI
# ndmiTransform = function(img) {
#     ndmi = img.normalizedDifference(['B4', 'B5']) # calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
#                 .multiply(1000) # scale results by 1000
#                 .select([0], ['NDMI']) # name the band
#                 .set('system:time_start', img.get('system:time_start'))
#    return ndmi
#}


## CALCULATE A GIVEN INDEX
# calcIndex = function(img, index, flip){
#  # make sure index string in upper case
#  index = index.toUpperCase()
#  
#  # figure out if we need to calc tc
#   tcList = ['TCB', 'TCG', 'TCW', 'TCA']
#   doTC = tcList.indexOf(index)
#  if(doTC >= 0){
#     tc = tcTransform(img)
#  }
#  
#  # need to flip some indices if this is intended for segmentation
#   indexFlip = 1
#  if(flip == 1){
#    indexFlip = -1
#  }
#  
#  # need to cast raw bands to float to make sure that we don't get errors regarding incompatible bands
#  # ...derivations are already float because of division or multiplying by decimal
#   indexImg
#  switch (index){
#    case 'B1':
#      indexImg = img.select(['B1']).float()#.multiply(indexFlip)
#      break
#    case 'B2':
#      indexImg = img.select(['B2']).float()#.multiply(indexFlip)
#      break
#    case 'B3':
#      indexImg = img.select(['B3']).float()#.multiply(indexFlip)
#      break
#    case 'B4':
#      indexImg = img.select(['B4']).multiply(indexFlip).float()
#      break
#    case 'B5':
#      indexImg = img.select(['B5']).float()#.multiply(indexFlip)
#      break
#    case 'B7':
#      indexImg = img.select(['B7']).float()#.multiply(indexFlip)
#      break
#    case 'NBR':
#      indexImg = nbrTransform(img).multiply(indexFlip)
#      break
#    case 'NDMI':
#      indexImg = ndmiTransform(img).multiply(indexFlip)
#      break
#    case 'NDVI':
#      indexImg = ndviTransform(img).multiply(indexFlip)
#      break
#    case 'NDSI':
#      indexImg = ndsiTransform(img).multiply(indexFlip)
#      break
#    case 'TCB':
#      indexImg = tc.select(['TCB'])#.multiply(indexFlip)
#      break
#    case 'TCG':
#      indexImg = tc.select(['TCG']).multiply(indexFlip)
#      break
#    case 'TCW':
#      indexImg = tc.select(['TCW']).multiply(indexFlip)
#      break
#    case 'TCA':
#      indexImg = tc.select(['TCA']).multiply(indexFlip)
#      break
#    default:
#      print('The index you provided is not supported')
#  }
#
#  return indexImg.set('system:time_start', img.get('system:time_start'))
#}

#exports.calcIndex = calcIndex

# MAKE AN LT STACK
# makeLTstack = function(img){
#   allStack = calcIndex(img, index, 1) # calc index for segmentation
#   ftvimg
#  for( ftv in ftvList){
#    ftvimg = calcIndex(img, ftvList[ftv], 0) # calc index for FTV
#                      .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()])
#    
#    allStack = allStack.addBands(ftvimg)
#                       .set('system:time_start', img.get('system:time_start'))
#  }
#  
#  return allStack
#}
#
#
#
## arrange collection as an annual stack
## TODO: would be nice to name the bands with the original band name and the year 
# TScollectionToStack = function(collection, startYear, endYear){
#   collectionArray = collection.toArrayPerBand()
#   nBands = ee.Image(collection.first()).bandNames().getInfo().length#
#   years = ltgee.getYearBandNames(startYear, endYear)
#   allStack = ee.Image()
#  for ( i = 0 i < nBands i++){
#     bandTS = collectionArray.select([i]).arrayFlatten([years])
#    allStack = ee.Image.cat([allStack, bandTS])
#  }
#                     
#  return allStack.slice(1,null)
#}
#
#exports.TScollectionToStack = TScollectionToStack
#
#
#
#
# standardize = function(collection){
#   mean = collection.reduce(ee.Reducer.mean())
#   stdDev = collection.reduce(ee.Reducer.stdDev())
#  
#   meanAdj = collection.map(function(img){
#    return img.subtract(mean).set('system:time_start', img.get('system:time_start'))
#  })
#  
#  return meanAdj.map(function(img){
#    return img.divide(stdDev).set('system:time_start', img.get('system:time_start'))
#  })
#}
#
#
## STANDARDIZE TASSELED CAP BRIGHTNESS GREENNESS WETNESS AND REDUCE THE COLLECTION
# makeTCcomposite = function(annualSRcollection, reducer){
#   TCcomposite = annualSRcollection.map(function(img){
#     tcb = calcIndex(img, 'TCB', 1)#.unmask(0)
#     tcg = calcIndex(img, 'TCG', 1)#.unmask(0)
#     tcw = calcIndex(img, 'TCW', 1)#.unmask(0)
#    return tcb.addBands(tcg)
#              .addBands(tcw)
#              .set('system:time_start', img.get('system:time_start'))
#  })
#  
#   tcb = TCcomposite.select(['TCB'])
#   tcg = TCcomposite.select(['TCG'])
#   tcw = TCcomposite.select(['TCW'])
#    
#  # standardize the TC bands
#   tcbStandard = standardize(tcb)
#   tcgStandard = standardize(tcg)
#   tcwStandard = standardize(tcw)
#  
#  # combine the standardized TC band collections into a single collection
#   tcStandard = tcbStandard.combine(tcgStandard).combine(tcwStandard)
#  
#  
#  TCcomposite = tcStandard.map(function(img){
#     imgCollection = ee.ImageCollection.fromImages(
#      [
#        img.select(['TCB'],['Z']),
#        img.select(['TCG'],['Z']),
#        img.select(['TCW'],['Z'])
#      ]
#    )
#     reducedImg
#    switch(reducer){
#      case 'mean':
#        reducedImg = imgCollection.mean()
#        break
#      case 'max':
#        reducedImg = imgCollection.max()
#        break
#      case 'sum':
#        reducedImg = imgCollection.sum()
#        break
#      default:
#        print('The reducer you provided is not supported')
#    }
#    
#    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'))            
#  })  
#  
#  
#  return TCcomposite
#}
#
#
## STANDARDIZE B5, TCB, TCG, NBR AND REDUCE THE COLLECTION
# makeEnsemblecomposite = function(annualSRcollection, reducer){
#  # make a collection of the ensemble indices stacked as bands
#   TCcomposite = annualSRcollection.map(function(img){
#     b5   = calcIndex(img, 'B5', 1)
#     b7   = calcIndex(img, 'B7', 1)
#     tcw  = calcIndex(img, 'TCW', 1)
#     tca  = calcIndex(img, 'TCA', 1)
#     ndmi = calcIndex(img, 'NDMI', 1)
#     nbr  = calcIndex(img, 'NBR', 1)
#    return b5.addBands(b7)
#             .addBands(tcw)
#             .addBands(tca)
#             .addBands(ndmi)
#             .addBands(nbr)
#             .set('system:time_start', img.get('system:time_start'))
#  })
#  
#  # make subset collections of each index
#   b5 = TCcomposite.select('B5')
#   b7 = TCcomposite.select('B7')
#   tcw = TCcomposite.select('TCW')
#   tca = TCcomposite.select('TCA')
#   ndmi = TCcomposite.select('NDMI')
#   nbr = TCcomposite.select('NBR')
#    
#  # standardize each index - get z-score
#   b5Standard = standardize(b5)
#   b7Standard = standardize(b7)
#   tcwStandard = standardize(tcw)
#   tcaStandard = standardize(tca)
#   ndmiStandard = standardize(ndmi)
#   nbrStandard = standardize(nbr)
#    
#  # combine the standardized TC band collections into a single collection
#   tcStandard = b5Standard.combine(b7Standard).combine(tcwStandard).combine(tcaStandard)
#                             .combine(ndmiStandard).combine(nbrStandard)
#  
#  # reduce the collection to a single value
#  TCcomposite = tcStandard.map(function(img){
#     imgCollection = ee.ImageCollection.fromImages(
#      [
#        img.select(['B5'],['Z']),
#        img.select(['B7'],['Z']),
#        img.select(['TCW'],['Z']),
#        img.select(['TCA'],['Z']),
#        img.select(['NDMI'],['Z']),
#        img.select(['NBR'],['Z']),
#      ]
#    )
#     reducedImg
#    switch(reducer){
#      case 'mean':
#        reducedImg = imgCollection.mean()
#        break
#      case 'max':
#        reducedImg = imgCollection.max()
#        break
#      case 'sum':
#        reducedImg = imgCollection.sum()
#        break
#      default:
#        print('The reducer you provided is not supported')
#    }
#    
#    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'))            
#  })  
#   
#  return TCcomposite
#}
#
#
## STANDARDIZE B5, TCB, TCG, NBR AND REDUCE THE COLLECTION
# makeEnsemblecomposite1 = function(annualSRcollection, reducer){
#  # make a collection of the ensemble indices stacked as bands
#   TCcomposite = annualSRcollection.map(function(img){
#     b5   = calcIndex(img, 'B5', 1)
#     tcb  = calcIndex(img, 'TCB', 1)
#     tcg  = calcIndex(img, 'TCG', 1)
#     nbr  = calcIndex(img, 'NBR', 1)
#    return b5.addBands(tcb)
#             .addBands(tcg)
#             .addBands(nbr)
#             .set('system:time_start', img.get('system:time_start'))
#  })
#  
#  # make subset collections of each index
#   b5 = TCcomposite.select('B5')
#   tcb = TCcomposite.select('TCB')
#   tcg = TCcomposite.select('TCG')
#   nbr = TCcomposite.select('NBR')
#    
#  # standardize each index - get z-score
#   b5Standard = standardize(b5)
#   tcbStandard = standardize(tcb)
#   tcgStandard = standardize(tcg)
#   nbrStandard = standardize(nbr)
#    
#  # combine the standardized TC band collections into a single collection
#   tcStandard = b5Standard.combine(tcbStandard).combine(tcgStandard).combine(nbrStandard)
#  
#  # reduce the collection to a single value
#  TCcomposite = tcStandard.map(function(img){
#     imgCollection = ee.ImageCollection.fromImages(
#      [
#        img.select(['B5'],['Z']),#.pow(ee.Image(1)).multiply(img.select('B5').gte(0).where(img.select('B5').lt(0),-1)),
#        img.select(['TCB'],['Z']),#.pow(ee.Image(1.5)).multiply(img.select('TCB').gte(0).where(img.select('TCB').lt(0),-1)),
#        img.select(['TCG'],['Z']),#.pow(ee.Image(1.5)).multiply(img.select('TCG').gte(0).where(img.select('TCG').lt(0),-1)),
#        img.select(['NBR'],['Z'])#.pow(ee.Image(1.5)).multiply(img.select('NBR').gte(0).where(img.select('NBR').lt(0),-1))
#      ]
#    )
#     reducedImg
#    switch(reducer){
#      case 'mean':
#        reducedImg = imgCollection.mean()
#        break
#      case 'max':
#        reducedImg = imgCollection.max()
#        break
#      case 'sum':
#        reducedImg = imgCollection.sum()
#        break
#      default:
#        print('The reducer you provided is not supported')
#    }
#    
#    return reducedImg.multiply(1000).set('system:time_start', img.get('system:time_start'))            
#  })  
#   
#  return TCcomposite
#}
#
#
#
## STANDARDIZE A INDEX INDEX - all disturbances are up
# standardizeIndex = function(collection, index){
#   zCollection = collection.map(function(img){
#    return calcIndex(img, index, 1)  
#  })
#
#  zCollection = standardize(zCollection)
#
#  zCollection = zCollection.map(function(img){
#    return img.multiply(1000).set('system:time_start', img.get('system:time_start'))  
#  })
#
#  return zCollection
#}
#
#
## BUILD AN LT COLLECTION
# buildLTcollection = function(collection, index, ftvList){
#   LTcollection
#  switch(index){
#    # tasseled cap composite
#    case 'TCC':
#      LTcollection = makeTCcomposite(collection, 'mean') 
#      break
#    case 'TCM':
#      LTcollection = makeTCcomposite(collection, 'max') 
#      break
#    case 'TCS':
#      LTcollection = makeTCcomposite(collection, 'sum') 
#      break
#    
#    # 6-band composite - Based on Figure 3 of the linked paper: https:#larse.forestry.oregonstate.edu/sites/larse/files/pub_pdfs/Cohen_et_al_2018.pdf
#    case 'ENC':
#      LTcollection = makeEnsemblecomposite(collection, 'mean') 
#      break
#    case 'ENM':
#      LTcollection = makeEnsemblecomposite(collection, 'max') 
#      break
#    case 'ENS':
#      LTcollection = makeEnsemblecomposite(collection, 'sum') 
#      break
#    
#    # 6-band composite - Based on Table 5 of the linked paper: https:#larse.forestry.oregonstate.edu/sites/larse/files/pub_pdfs/Cohen_et_al_2018.pdf
#    case 'ENC1':
#      LTcollection = makeEnsemblecomposite1(collection, 'mean') 
#      break
#    case 'ENM1':
#      LTcollection = makeEnsemblecomposite1(collection, 'max') 
#      break
#    case 'ENS1':
#      LTcollection = makeEnsemblecomposite1(collection, 'sum') 
#      break
#    
#    # standardized versions of indices: mean 0 stdDev 1  
#    case 'B5z':
#      LTcollection = standardizeIndex(collection, 'B5') 
#      break
#    case 'B7z':
#      LTcollection = standardizeIndex(collection, 'B7') 
#      break
#    case 'TCWz':
#      LTcollection = standardizeIndex(collection, 'TCW') 
#      break
#    case 'TCAz':
#      LTcollection = standardizeIndex(collection, 'TCA') 
#      break
#    case 'NDMIz':
#      LTcollection = standardizeIndex(collection, 'NDMI') 
#      break
#    case 'NBRz':
#      LTcollection = standardizeIndex(collection, 'NBR') 
#      break      
#      
#    default:
#      LTcollection = collection.map(function(img){
#               allStack = calcIndex(img, index, 1)
#               ftvimg
#              for( ftv in ftvList){
#                ftvimg = calcIndex(img, ftvList[ftv], 0)
#                                     .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()])
#      
#                  allStack = allStack.addBands(ftvimg)
#                                     .set('system:time_start', img.get('system:time_start'))
#              }
#              return allStack
#      })
#    }
#  
#  return LTcollection
#}
#
#exports.buildLTcollection = buildLTcollection
#
#
## THIS BUILDS AN ANNUAL-SR-TRANSFORMED COLLECTION FOR USE BY TIMESYNC-LEGACY
# buildTSindexCollection = function(collection, ftvList){
#  return collection.map(function(img){
#             allStack = ee.Image()
#             ftvimg
#            for( ftv in ftvList){
#              ftvimg = calcIndex(img, ftvList[ftv], 0)
#                                   .select([ftvList[ftv]],['ftv_'+ftvList[ftv].toLowerCase()])
#                                   .unmask(0)
#    
#                allStack = allStack.addBands(ftvimg)
#                                   .set('system:time_start', img.get('system:time_start'))
#            }
#            return allStack.slice(1,null)
#         })
#}
#
#
#
## TRANSFORM AN ANNUAL SR COLLECTION TO AN ANNUAL COLLECTION OF SELECTED INDICES OR BANDS
# transformSRcollection = function(srCollection, bandList){
#  return srCollection.map(function(img){
#     allStack = calcIndex(img, bandList[0], 0)
#    for( band=1 band < bandList.length band++){
#       bandImg = calcIndex(img, bandList[band], 0)
#      allStack = allStack.addBands(bandImg)
#                         .set('system:time_start', img.get('system:time_start'))
#    }
#    return allStack
#  })
#}
#
#exports.transformSRcollection = transformSRcollection
#
#
#
#exports.runLT = function(startYear, endYear, startDay, endDay, aoi, index, ftvList, runParams){
#   annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi)
#   annualLTcollection = buildLTcollection(annualSRcollection, index, ftvList)
#  runParams.timeSeries = annualLTcollection
#  return ee.Algorithms.TemporalSegmentation.LandTrendr(runParams)
#}
#
#
#
# getFittedData = function(lt, startYear, endYear, index){
#   bandNames = getYearBandNames(startYear, endYear)
#  return lt.select('ftv_'+index.toLowerCase()+'_fit')
#           .arrayFlatten([bandNames])
#}
#
#exports.getFittedData = getFittedData
#
#
#
## MAKE RGB COMPOSITE
# makeRGBcomposite = function(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev){
#   annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi)
#   annualLTcollection = buildLTcollection(annualSRcollection, index, [index])
#  runParams.timeSeries = annualLTcollection
#   lt = ee.Algorithms.TemporalSegmentation.LandTrendr(runParams)
#   bandNames = getYearBandNames(startYear, endYear)
#   ftvStack = getFittedData(lt, startYear, endYear, index)
#  return ftvStack.select([redYear.toString(),greenYear.toString(),blueYear.toString()]) 
#}
#
#exports.makeRGBcomposite = makeRGBcomposite
#
#
#
#exports.mapRGBcomposite = function(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev){
#   rgb = makeRGBcomposite(index, startYear, endYear, startDay, endDay, redYear, greenYear, blueYear, aoi, runParams, nStdev)
#   stretch = stdDevStretch(rgb, aoi, nStdev)
#   rgbVis = rgb.visualize({min: stretch[0], max: stretch[1], gamma: 1}).clip(aoi) 
#  return rgbVis
#  #Map.centerObject(aoi, 10)
#  #Map.addLayer(rgbVis)
#}
#
#
#
## FLATTEN SOURCE COLLECTION INTO A GIANT STACK FOR EXPORT - USED BY TIMESYNC
#exports.timesyncLegacyStack = function(startYear, endYear, startDay, endDay, aoi){
#   annualSRcollection = buildSRcollection(startYear, endYear, startDay, endDay, aoi)
#   annualTScollection = buildTSindexCollection(annualSRcollection, ['B1','B2','B3','B4','B5','B7','TCB','TCG','TCW'])
#  
#  # THE FOLLOWING COULD BE A GENERIC FUNCTION
#   annualTSArray = annualTScollection.toArrayPerBand()
#   years = []
#   allStack = ee.Image()
#  for ( i = startYear i <= endYear ++i) years.push(i.toString())
#  for (i = 0 i <= 8 i++){
#     bandTS = annualTSArray.select([i]).arrayFlatten([years])
#    allStack = ee.Image.cat([allStack, bandTS])
#  }
#                     
#  return allStack.slice(1,null).toShort()
#  #######################  
#}







## #######################################################################################
## ###### PIXEL TIME SERIES FUNCTIONS ####################################################
## #######################################################################################
#
## TODO: convert to python
#
## RETURN LT RESULTS FOR A SINGLE PIXEL AS AN OBJECT
# ltPixelTimeSeries = function(img, pixel) {
#  return img.reduceRegion({
#   reducer: 'first',
#   geometry: pixel,
#   scale: 30
#  }).getInfo()
#}
#
#
#
## PARSE OBJECT RETURNED FROM 'getPoint' TO ARRAY OF SOURCE AND FITTED
#exports.ltPixelTimeSeriesArray = function(lt, pixel, indexFlip){
#   pixelTS = ltPixelTimeSeries(lt, pixel)
#   data = [['Year', 'Original', 'Fitted']]
#   len = pixelTS.LandTrendr[0].length
#  for ( i = 0 i < len i++) {
#    data = data.concat([[pixelTS.LandTrendr[0][i], pixelTS.LandTrendr[1][i]*indexFlip, pixelTS.LandTrendr[2][i]*indexFlip]])
#  }
#  return {ts:data, rmse:pixelTS.rmse}
#}




# #######################################################################################
# ###### LT ARRAY FUNCTIONS ###############################################################
# #######################################################################################

# get segment info array
def getSegmentData(lt, index):
    lt = lt.select('LandTrendr')                # select the LandTrendr band
    # slice out the 'Is Vertex' row - yes(1)/no(0)
    vertexMask = lt.arraySlice(0, 3, 4)
    # use the 'Is Vertex' row as a mask for all rows
    vertices = lt.arrayMask(vertexMask)
    # slice out the vertices as the start of segments
    left = vertices.arraySlice(1, 0, -1)
    # slice out the vertices as the end of segments
    right = vertices.arraySlice(1, 1)
    # get year dimension of LT data from the segment start vertices
    startYear = left.arraySlice(0, 0, 1)
    # get spectral index dimension of LT data from the segment start
    # vertices
    startVal = left.arraySlice(0, 2, 3)
    # get year dimension of LT data from the segment end vertices
    endYear = right.arraySlice(0, 0, 1)
    # get spectral index dimension of LT data from the segment end vertices
    endVal = right.arraySlice(0, 2, 3)
    # subtract the segment start year from the segment end year to
    # calculate the duration of segments
    dur = endYear.subtract(startYear)
    # substract the segment start index value from the segment end index
    # value to calculate the delta of segments
    mag = endVal.subtract(startVal)
    # calculate the rate of spectral change
    rate = mag.divide(dur)

    if indexFlipper(index) == -1:
        startVal = startVal.multiply(-1)
        endVal = endVal.multiply(-1)

    return (ee.Image.cat([startYear.add(1), endYear, startVal, endVal, mag, dur, rate])
            .toArray(0)
            .mask(vertexMask.mask()))



# #######################################################################################
# ###### HELPER FUNCTIONS ###############################################################
# #######################################################################################

# GET A SERIES OF BANDS NAMES AS YEARS
def getYearBandNames(startYear, endYear):
    years = [str(i) for i in range(startYear, endYear + 1)]
    return years


## STANDARD DEVIATION STRETCH
# stdDevStretch = function(img, aoi, nStdev){
#   mean = img.reduceRegion({reducer:ee.Reducer.mean(), geometry:aoi, scale:900, bestEffort:true, maxPixels:1e4})
#              .toArray()
#              .reduce(ee.Reducer.mean(), [0])
#
#   stdDev = img.reduceRegion({reducer:ee.Reducer.stdDev(), geometry:aoi, scale:900, bestEffort:true, maxPixels:1e4})
#                .toArray()
#                .reduce(ee.Reducer.mean(), [0])
#                .multiply(nStdev)
#
#   max = mean.add(stdDev).getInfo()[0]
#   min = mean.subtract(stdDev).getInfo()[0]
#  return [min, max]
#}

# INDEX FLIPPER
def indexFlipper(index):
    indexObj = {'NBR': -1,
                'NDVI': -1,
                'NDSI': -1,
                'TCB': 1,
                'TCG': -1,
                'TCW': -1,
                'B1': 1,
                'B2': 1,
                'B3': 1,
                'B4': -1,
                'B5': 1,
                'B7': 1,
                'ENC': 1,
                'ENC1': 1,
                'TCC': 1,
                'NBRz': 1,
                'B5z': 1}
    
    return indexObj[index]

def makeBoolean(value):
    if not isinstance(value, bool):
        value = value.lower() != 'false'
    return value







## #######################################################################################
## ###### UI STUFF #######################################################################
## #######################################################################################
#
## TODO: delete or copy import parts such as default params
#
## LT PARAMS
#exports.paramPanel = function(){
#
# runParams = [
#  {label: 'Max Segments:', value: 6},
#  {label: 'Spike Threshold:', value: 0.9},
#  {label: 'Vertex Count Overshoot:', value: 3},
#  {label: 'Prevent One Year Recovery:', value: true},
#  {label: 'Recovery Threshold:', value: 0.25},
#  {label: 'p-value Threshold:', value: 0.05},
#  {label: 'Best Model Proportion:', value: 0.75},
#  {label: 'Min Observations Needed:', value: 6},
#]
#
# paramBoxes = []
# paramPanels = [ui.Label('Define Segmentation Parameters',{fontWeight: 'bold'})]
#runParams.forEach(function(param, index){
#   paramLabel = ui.Label(param.label)
#   paramBox = ui.Textbox({value:param.value})
#  paramBox.style().set('stretch', 'horizontal')
#   paramPanel = ui.Panel([paramLabel,paramBox], ui.Panel.Layout.Flow('horizontal'))
#  paramBoxes.push(paramBox)
#  paramPanels.push(paramPanel)
#})
#
#  return ui.Panel(paramPanels, null, {stretch: 'horizontal'})
#}
#
#exports.getParams = function(paramPanel){
#   prevOneYrRec = paramPanel.widgets().get(4).widgets().get(1).getValue()
#  
#  prevOneYrRec = makeBoolean(prevOneYrRec)
#  
#  #if(typeof(prevOneYrRec) !== "boolean"){
#  #  prevOneYrRec = prevOneYrRec.toLowerCase() != 'false'
#  #}
#  
#  return { 
#    maxSegments:              parseInt(paramPanel.widgets().get(1).widgets().get(1).getValue()),
#    spikeThreshold:         parseFloat(paramPanel.widgets().get(2).widgets().get(1).getValue()),
#    vertexCountOvershoot:     parseInt(paramPanel.widgets().get(3).widgets().get(1).getValue()),
#    preventOneYearRecovery:                                                        prevOneYrRec,
#    recoveryThreshold:      parseFloat(paramPanel.widgets().get(5).widgets().get(1).getValue()),
#    pvalThreshold:          parseFloat(paramPanel.widgets().get(6).widgets().get(1).getValue()),
#    bestModelProportion:    parseFloat(paramPanel.widgets().get(7).widgets().get(1).getValue()),
#    minObservationsNeeded:    parseInt(paramPanel.widgets().get(8).widgets().get(1).getValue())
#  }
#}
#
#
#
## INDEX PANEL
#exports.indexSelectPanel = function(){
#   indexLabel = ui.Label('Select Index',{fontWeight: 'bold'})
#   indexList = ['NBR','NDVI','NDSI','TCB','TCG','TCW','B1','B2','B3','B4','B5','B7','ENC','ENC1','TCC']
#   indexSelect = ui.Select({items:indexList, value:'NBR', style:{stretch: 'horizontal'}})
#  return ui.Panel([indexLabel,indexSelect], null, {stretch: 'horizontal'})
#}
#
#exports.getIndexSelect = function(indexSelectPanel){
#  return indexSelectPanel.widgets().get(1).getValue()
#}
#
#
#
## YEAR PANEL
#exports.yearPanel = function(){
#   yearSectionLabel = ui.Label('Define Year Range',{fontWeight: 'bold'})
#  
#   startYearLabel = ui.Label('Start Year:')
#   startYearslider = ui.Slider({min:1984, max:2018, value:1984, step:1})
#  startYearslider.style().set('stretch', 'horizontal')
#  
#   endYearLabel = ui.Label('End Year:')
#   endYearslider = ui.Slider({min:1984, max:2018, value:2017, step:1})
#  endYearslider.style().set('stretch', 'horizontal')
#  
#  return ui.Panel(
#    [
#      yearSectionLabel,
#      ui.Panel([startYearLabel, startYearslider], ui.Panel.Layout.Flow('horizontal'), {stretch: 'horizontal'}), #
#      ui.Panel([endYearLabel  , endYearslider], ui.Panel.Layout.Flow('horizontal'), {stretch: 'horizontal'})
#    ] 
#  )
#}
#
#exports.getYears = function(yearPanel){
#  return {
#    startYear:yearPanel.widgets().get(1).widgets().get(1).getValue(),
#    endYear:yearPanel.widgets().get(2).widgets().get(1).getValue()
#  }
#}
#
#
#
#
#
## DATE PANEL
#exports.datePanel = function(){
#   dateSectionLabel = ui.Label('Define Date Range (month-day)',{fontWeight: 'bold'})
#   startDayLabel = ui.Label('Start Date:')
#   startDayBox = ui.Textbox({value:'06-10'})
#  startDayBox.style().set('stretch', 'horizontal')
#  
#   endDayLabel = ui.Label('End Date:')
#   endDayBox = ui.Textbox({value:'09-20'})
#  endDayBox.style().set('stretch', 'horizontal')
#  
#  return ui.Panel(
#    [
#      dateSectionLabel,
#      ui.Panel(
#        [startDayLabel, startDayBox, endDayLabel, endDayBox],
#        ui.Panel.Layout.Flow('horizontal'), {stretch: 'horizontal'}
#      )
#    ]
#  )
#}
#
#exports.getDays = function(datePanel){
#  return {
#    startDay:datePanel.widgets().get(1).widgets().get(1).getValue(),
#    endDay:datePanel.widgets().get(1).widgets().get(3).getValue()
#  }
#}
#
#
## COORDINATE PANEL
#exports.coordsPanel = function(){
#   coordSectionLabel = ui.Label('Define Pixel Coordinates (optional)',{fontWeight: 'bold'})
#  
#   latLabel = ui.Label('Latitude:')
#   latBox = ui.Textbox({value:43.7929})
#  latBox.style().set('stretch', 'horizontal')
#  
#   lonLabel = ui.Label('Longitude:')
#   lonBox = ui.Textbox({value:-122.8848})
#  lonBox.style().set('stretch', 'horizontal')
#  
#  return ui.Panel(
#    [
#      coordSectionLabel,
#      ui.Panel([lonLabel, lonBox, latLabel, latBox],ui.Panel.Layout.Flow('horizontal'))
#    ],
#    null,
#    {stretch: 'horizontal'}
#  )
#}
#
#exports.getCoords = function(coordsPanel){
#    return {
#      lon:parseFloat(coordsPanel.widgets().get(1).widgets().get(1).getValue()),
#      lat:parseFloat(coordsPanel.widgets().get(1).widgets().get(3).getValue())
#  }
#}
#
#
#
## BUFFER PANEL
#exports.bufferPanel = function(){
#   bufferSectionLabel = ui.Label('Define a Buffer Around Point (km)',{fontWeight: 'bold'})
#   bufferBoxLabel = ui.Label('Buffer:')
#   bufferBox = ui.Textbox({value: 50, style:{stretch: 'horizontal'}})
#  return ui.Panel(
#    [
#      bufferSectionLabel,
#      ui.Panel([bufferBoxLabel,bufferBox], ui.Panel.Layout.Flow('horizontal'), {stretch: 'horizontal'})
#    ]
#  )
#}
#
#exports.getBuffer = function(bufferPanel){
#  return bufferPanel.widgets().get(1).widgets().get(1).getValue()
#}
#
## SUBMIT BUTTON
#exports.submitButton = function(){
#  return ui.Button({label: 'Submit', style:{stretch: 'horizontal'}})
#}
#
#
#
#
#
#exports.disturbanceMap = function(distParams){
#/*
#  distParams.type = 'Greatest'
#  distParams.segInfo = result from setInfo function
#  distParams.year = {checked:false, start:null, end:null}
#  distParams.mag = {checked:false, year1:null, year20:null}
#  distParams.preval = {checked:false, value:null}
#  distParams.mmu = {checked:false, value:null}
#*/  
#
#   sortByThis
#   segInfoSorted
#  switch (distParams.type){
#    case 'Greatest':
#      sortByThis = distParams.segInfo.arraySlice(0,4,5).toArray(0).multiply(-1) # need to flip the delta here, since arraySort is working by ascending order
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by magnitude
#
#      break
#    case 'Least':
#      sortByThis = distParams.segInfo.arraySlice(0,4,5).toArray(0)
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by magnitude
#
#      break
#    case 'Newest':
#      sortByThis = distParams.segInfo.arraySlice(0,0,1).toArray(0).multiply(-1) # need to flip the delta here, since arraySort is working by ascending order
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by startYear
#      
#      break
#    case 'Oldest':
#      sortByThis = distParams.segInfo.arraySlice(0,0,1).toArray(0)
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by startYear
#      
#      break
#    case 'Longest':
#      sortByThis = distParams.segInfo.arraySlice(0,5,6).toArray(0).multiply(-1) # need to flip the delta here, since arraySort is working by ascending order
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by duration
#      break
#    case 'Fastest':
#      sortByThis = distParams.segInfo.arraySlice(0,5,6).toArray(0)
#      segInfoSorted = distParams.segInfo.arraySort(sortByThis) # sort the array by duration
#      break
#  }
#
#   distArray = segInfoSorted.arraySlice(1, 0, 1) # get the first segment in the sorted array
#  # make an image from the array of attributes for the greatest disturbance
#   distImg = ee.Image.cat(distArray.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod']]).toShort(),
#                             distArray.arraySlice(0,4,5).arrayProject([1]).arrayFlatten([['mag']]),
#                             distArray.arraySlice(0,5,6).arrayProject([1]).arrayFlatten([['dur']]),
#                             distArray.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['preval']]),
#                             distArray.arraySlice(0,6,null).arrayProject([1]).arrayFlatten([['rate']]))
#  
#  
#  
#  # start a mask based on magnitudes greater than 0
#   mask = distImg.select('mag').gt(0)#.mask()
#  mask = ee.Image(1).mask(mask)
#  
#  # add additional masks depending on selection
#  # filter by year
#   addMask
#  if(makeBoolean(distParams.year.checked) === true){
#     yod = distImg.select('yod')
#    addMask = yod.gte(distParams.year.start).and(yod.lte(distParams.year.end))
#    mask = mask.updateMask(addMask)
#  }
#
#  # filter by mag
#  if(makeBoolean(distParams.mag.checked) === true){
#    addMask = distImg.select(['dur'])                        
#                     .multiply((distParams.mag.year20 - distParams.mag.year1) / 19.0)
#                     .add(distParams.mag.year1)                                     
#                     .lte(distImg.select(['mag']))
#    mask = mask.updateMask(addMask)
#  }
#  
#  # filter by preval
#  if(makeBoolean(distParams.preval.checked) === true){
#    if(indexFlipper(distParams.index) == -1){
#      addMask = distImg.select(['preval']).gte(distParams.preval.value)
#    } else{
#      addMask = distImg.select(['preval']).lte(distParams.preval.value)
#    }
#    mask = mask.updateMask(addMask)
#  }
#  
#  # filter by mmu    # this works 2018/06/05 - JDB
#  if(makeBoolean(distParams.mmu.checked) === true){ 
#    if(distParams.mmu.value > 1){
#      addMask = distImg.select(['yod'])           # patchify based on disturbances having the same year of detection
#                       .connectedPixelCount(distParams.mmu.value, true) # count the number of pixel in a candidate patch
#                       .gte(distParams.mmu.value)                      # are the the number of pixels per candidate patch greater than user-defined minimum mapping unit?
#      mask = mask.updateMask(addMask)     # mask the pixels/patches that are less than minimum mapping unit
#    }
#  }
#  
#  # apply the filter mask
#  return distImg.mask(mask)
#}






