# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:24:04 2017

@author: filip107
"""
import ee
import ee.mapclient
ee.Initialize()
from math import pi

# --------------- Prep ------------------------------------------------------#

def rename_bands(img):
    """
    Rename bands for an S2 image.
    
    This is should be applicable for the original S2 collections and the recommended
    L2A collection "COPERNICUS/S2_SR_HARMONIZED", which has additional bands.
    """
    rename_dict = ee.Dictionary(
                  {'B1':'cb',
                   'B2':'blue',
                   'B3':'green',
                   'B4':'red',
                   'B5':'re1',
                   'B6':'re2',
                   'B7':'re3',
                   'B8':'nir',
                   'B8A':'re4',
                   'B9':'waterVapor',
                   'B11':'swir1',
                   'B12':'swir2',
                   'AOT':'aot',
                   'WVP':'wvp',
                   'SCL':'scl',
                   'TCI_R':'tci_r',
                   'TCI_G':'tci_g',
                   'TCI_B':'tci_b',
                   'MSK_CLDPRB':'msk_cldprb',
                   'MSK_SNWPRB':'msk_snwprb',
                   'QA10':'qa10',
                   'QA20':'qa20',
                   'QA60':'qa60',
                   'MSK_CLASSI_OPAQUE':'msk_opaque',
                   'MSK_CLASSI_CIRRUS':'msk_cirrus',
                   'MSK_CLASSI_SNOW_ICE':'msk_snowice'
                   })
    bands = img.bandNames()
    newbands = rename_dict.values(bands)
    return img.select(bands, newbands)


def prep_l1c_img(img):
    """Rename and rescale all bands to TOA reflectance"""
    t = rename_bands(img).multiply(0.0001) # Rescale to 0-1
    out = t.copyProperties(img).copyProperties(img,['system:time_start']).set('date', img.date().format('YYYY-MM-dd'))
    return out

def prep_l2a_img(img):
    """Rename and rescale all (available) bands to surface reflectance"""
    img = rename_bands(img)
    
    optbands = ['cb', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'waterVapor', 'cirrus', 'swir1', 'swir2']
    optbands = img.bandNames().filter(ee.Filter.inList('item', optbands))
    opt = img.select(optbands).multiply(0.0001)
    
    atmbands = ['aot', 'wvp']
    atmbands = img.bandNames().filter(ee.Filter.inList('item', atmbands))
    atm = img.select(atmbands).multiply(0.001)
    
    qabands = ['qa60', 'scl', 'msk_cldprb', 'msk_snwprb']
    qabands = img.bandNames().filter(ee.Filter.inList('item', qabands))
    qa = img.select(qabands)
    
    t = ee.Image([opt, atm, qa])
    out = (t.copyProperties(img)
            .copyProperties(img,['system:time_start', 'system:time_end', 'system:footprint'])
            .set('date', img.date().format('YYYY-MM-dd')))
    return out

def prep_sr_harm_img(img):
    """Rename and rescale all (available) bands to surface reflectance for "COPERNICUS/S2_SR_HARMONIZED" collection."""
    img = rename_bands(img)
    
    optbands = ['cb', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'waterVapor', 'cirrus', 'swir1', 'swir2']
    optbands = img.bandNames().filter(ee.Filter.inList('item', optbands))
    opt = img.select(optbands).multiply(0.0001)
    
    tcibands = ['tci_r', 'tci_g', 'tci_b']
    tcibands = img.bandNames().filter(ee.Filter.inList('item', tcibands))
    tci = img.select(tcibands)
    
    atmbands = ['aot', 'wvp']
    atmbands = img.bandNames().filter(ee.Filter.inList('item', atmbands))
    atm = img.select(atmbands).multiply(0.001)
    
    qabands = ['qa60', 'scl', 'msk_cldprb', 'msk_snwprb', 'msk_opaque', 'msk_cirrus', 'msk_snowice']
    qabands = img.bandNames().filter(ee.Filter.inList('item', qabands))
    qa = img.select(qabands)
    
    t = ee.Image([opt, tci, atm, qa])
    out = (t.copyProperties(img)
            .copyProperties(img,['system:time_start', 'system:time_end', 'system:footprint'])
            .set('date', img.date().format('YYYY-MM-dd')))
    return out

#--------------- Cloud Masking ----------------------------------------------#

# CloudScore originally written by Matt Hancher and adapted for S2 data by Ian Housman
# User Params
cloudThresh = 20 # Ranges from 1-100.Lower value will mask more pixels out. Generally 10-30 works well with 20 being used most commonly 
cloudHeights = ee.List.sequence(200,10000,500) # Height of clouds to use to project cloud shadows (min of 500 used by Ian)
irSumThresh = 0.35 # Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
dilatePixels = 2 # Pixels to dilate around clouds (3)
contractPixels = 1 # Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors (2)


def _rescale(img, exp, thresholds):
    return (img.expression(exp, {'img': img})
               .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]))


def sentinelCloudScore(img):
    """
    Cloud masking algorithm for Sentinel2
    Built on ideas from Landsat cloudScore algorithm
    Currently in beta and may need tweaking for individual study areas
    
    Author: Matt Hancher, Ian Housman
    """
    # Compute several indicators of cloudyness and take the minimum of them.
    score = ee.Image(1)
    blueCirrusScore = ee.Image(0)

    # Clouds are reasonably bright in the blue or cirrus bands.
    # Use .max as a pseudo OR conditional
    blueCirrusScore = blueCirrusScore.max(_rescale(img, 'img.blue', [0.1, 0.5]))
    blueCirrusScore = blueCirrusScore.max(_rescale(img, 'img.cb', [0.1, 0.5]))
    blueCirrusScore = blueCirrusScore.max(_rescale(img, 'img.cirrus', [0.1, 0.3]))
  
    score = score.min(blueCirrusScore)

    # Clouds are reasonably bright in all visible bands.
    score = score.min(_rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]))
    
    # Clouds are reasonably bright in all infrared bands.
    score = score.min(_rescale(img, 'img.nir + img.swir1 + img.swir2', [0.3, 0.8]))

#    # Clouds are moist
#    ndmi = img.normalizedDifference(['nir', 'swir1'])
#    score = score.min(_rescale(ndmi, 'img', [-0.1, 0.1]))

    # Also mask snow (consider doing snow masking seperately as high ndsi & high vis)
    ndsi = img.normalizedDifference(['green', 'swir1'])
    score = score.max(_rescale(ndsi, 'img', [0, 1.0]))

    score = score.multiply(100).byte()

    return img.addBands(score.select([0], ['cloudScore']))  # rename(ee.List(['cloudScore'])))


def maskS2clouds(image):
    """ Function to mask clouds using the Sentinel-2 QA band.
    """
    qa = image.select('QA60').int16()

    # Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = 2**10
    cirrusBitMask = 2**11

    # Both flags should be set to zero, indicating clear conditions.
    mask = (qa.bitwiseAnd(cloudBitMask).eq(0).And(
            qa.bitwiseAnd(cirrusBitMask).eq(0)))

    # Return the masked and scaled data.
    return image.updateMask(mask)


def tdom2(imgs, sum_bands=['nir', 'swir1'], zscore_thresh=-1, sum_thresh=0.35,
         erode_pix=1.5, dilate_pix=3.5, mean_img=None, stddev_img=None,
         stat_imgs=None
        ):
    """ Apply Temporal Dark Outlier Mask (TDOM) to an image or image collection
    
    Identify dark outliers for an image(s) using a time series of images or 
    precomputed mean and standard deviations from a time series. Then mask out
    the dark outlier pixels. This is most useful for identifying cloud shadows 
    which may be missed by Fmask when used with nir and swir1 as the sum bands. 
    If the time series is too short this will tend to have higher 
    commission/omission error. The pixel-wise mean and standard deviation are
    taken in order of preference from mean_img & stddev_img > stat_imgs > imgs.
    
    Parameters
    ----------
    img: ee.Image or ee.ImageCollection
        Image or images for which to mask out dark outliers. If an image collection,
        it is used to calculate the mean and standard deviation for outlier
        detection only if precomputed mean_img & stddev_img or stat_imgs are
        not provided.
        
    sum_bands: list
        list of bands names to use identifying outliers
        
    zscore_thresh: float
        Number of standard deviations below the mean to consider as low outlier.
        
    sum_thresh: float
        Maximum brightness from sum of sum_bands to consider as a dark pixel.
        Default of 0.35 is based on use of reflectance rescaled 0-1 with default
        nir and swir1 bands.
        
    erode_pix: float
        Number of pixels to apply erosion to the detected dark pixels
    
    dilate_pix: float
        Number of pixels of dilation to apply to the mask. Dilation is applied
        after erosion so a morphological opening is applied to the mask if both 
        parameters are not zero.
        
    mean_img: ee.Image
        Image with mean of time series for the bands used in outlier detection.
        Must have the same band names as sum_bands.
        
    stddev_img: ee.Image
        Image with standard deviation of time series for the bands used in
        outlier detection. Must have the same band names as sum_bands.
        
    stat_imgs: ee.ImageCollection
        Images to use for calculating the pixel-wise means and standard deviations
        to use in outlier detection. Only used if mean_img and stddev_img not
        provided.
        
        
    Returns
    -------
    ee.image or ee.ImageCollection
        Collection of images with dark outliers masked out  
        
    Authors
    -------
    Original concept written by Carson Stam and adapted by Ian Housman.
    Ported to python and modified by Steven Filippelli.
        
    TODO: add error handling for imgs not provided, and other bad params
    """
    
    # Get mean and stddev of a time series for sum bands
    if not(mean_img and stddev_img):
        if stat_imgs:
            stddev_img = stat_imgs.select(sum_bands).reduce(ee.Reducer.stdDev())
            mean_img = stat_imgs.select(sum_bands).mean()
        
        elif type(imgs) is ee.ImageCollection:
            stddev_img = imgs.select(sum_bands).reduce(ee.Reducer.stdDev())
            mean_img = imgs.select(sum_bands).mean()
        
        else:
            # TODO: raise error since stats can't be computed from anything
            return imgs
    
    # function to mask dark outliers for a single image
    def dark_outliers(img):
        z = img.select(sum_bands).subtract(mean_img).divide(stddev_img)
        img_sum = img.select(sum_bands).reduce(ee.Reducer.sum())
        tdom_mask = z.lt(zscore_thresh).reduce(ee.Reducer.sum()).eq(len(sum_bands)).And(img_sum.lt(sum_thresh))
        
        # Erosion then dilate is morphological opening if both params not 0
        # TODO: I think fast distance transform may be faster for erosion & dilation
        if erode_pix:
            tdom_mask = tdom_mask.focal_min(erode_pix)
        if dilate_pix:
            tdom_mask = tdom_mask.focal_max(dilate_pix)
        
        return img.updateMask(tdom_mask.Not())
    
    # apply mask to single image or image collection
    if type(imgs) is ee.ImageCollection:
        imgs = imgs.map(dark_outliers)
    else:
        imgs = dark_outliers(imgs)
    
    return imgs


def projectShadows(cloudMask,image,cloudHeights):
    """
    Implementation of Basic cloud shadow shift

    Author: Gennadii Donchyts
    License: Apache 2.0
    """
    meanAzimuth = image.get('MEAN_SOLAR_AZIMUTH_ANGLE')
    meanZenith = image.get('MEAN_SOLAR_ZENITH_ANGLE')

    # Find dark pixels
    darkPixels = (image.select(['nir', 'swir1', 'swir2']).reduce(ee.Reducer.sum()).lt(irSumThresh)
                  .focal_min(contractPixels).focal_max(dilatePixels))  # .gte(1)

    # Get scale of image
    nominalScale = cloudMask.projection().nominalScale()

    # Find where cloud shadows should be based on solar geometry
    # Convert to radians
    azR = ee.Number(meanAzimuth).multiply(pi).divide(180.0).add(ee.Number(0.5).multiply(pi))
    zenR = ee.Number(0.5).multiply(pi).subtract(ee.Number(meanZenith).multiply(pi).divide(180.0))
    
    #Find the shadows
    def get_shadows(cloudHeight):
        cloudHeight = ee.Number(cloudHeight)
        
        shadowCastedDistance = zenR.tan().multiply(cloudHeight)#Distance shadow is cast
        x = azR.cos().multiply(shadowCastedDistance).divide(nominalScale).round()#X distance of shadow
        y = azR.sin().multiply(shadowCastedDistance).divide(nominalScale).round()#Y distance of shadow
        height_shadow = cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y))
        return height_shadow

    shadows = cloudHeights.map(get_shadows)

    shadowMask = ee.ImageCollection.fromImages(shadows).max()
 
    #Create shadow mask
    shadowMask = shadowMask.And(cloudMask.Not())
    shadowMask = shadowMask.And(darkPixels)
    
    cloudShadowMask = shadowMask.Or(cloudMask)
    
    image = image.updateMask(cloudShadowMask.Not()).addBands(shadowMask.rename(['cloudShadowMask']))
    return image


def uniqueValues(collection,field):
    """ Function to find unique values of a field in a collection
    Author: Ian Housman
    """
    values = ee.Dictionary(collection.reduceColumns(ee.Reducer.frequencyHistogram(),[field]).get('histogram')).keys()
    return values


def dailyMosaics(imgs):
    """ Function to simplify data into daily mosaics
    Author: Ian Housman
    """
    # Simplify date to exclude time of day
    def simplify_date(img):
        d = ee.Date(img.get('system:time_start'))
        day = d.get('day')
        m = d.get('month')
        y = d.get('year')
        simpleDate = ee.Date.fromYMD(y,m,day)
        return img.set('simpleTime',simpleDate.millis())

    imgs = imgs.map(simplify_date)

    # Find the unique days
    days = uniqueValues(imgs,'simpleTime')
    
    def unique_days(d):
        d = ee.Number.parse(d)
        d = ee.Date(d)
        t = imgs.filterDate(d,d.advance(1,'day'))
        f = ee.Image(t.first())
        t = t.mosaic()
        t = t.set('system:time_start',d.millis())
        t = t.copyProperties(f)
        return t
  
    imgs = days.map(unique_days)
    imgs = ee.ImageCollection.fromImages(imgs)
    
    return imgs


def dailyMosaics2(imgs):
    """ 
    Method for daily mosaics taken from - https://code.earthengine.google.com/ec5f79e5bd195ec0da0ddbdc2cd11088
    This method may be more efficient than GTAC's method (dailyMosaic) according to a brief Profiler test
    """
    # set date on each image
    def set_date(img):
        date = img.date().format('YYYY-MM-dd')
        return img.set('date', date)
    
    imgs = imgs.map(set_date)
  
    # 'distinct' removes duplicates from a collection based on a property.
    distinctDates_S2_sr = imgs.distinct('date').sort('date')

    # define the filter
    filt = ee.Filter.equals(leftField='date', rightField='date')

    # 'ee.Join.saveAll' Returns a join that pairs each element from the first collection with a group of matching elements from the second collection
    # the matching images are stored in a new property called 'date_match'
    join = ee.Join.saveAll('date_match')

    # 'apply' Joins to collections.
    joinCol_S2_sr = join.apply(distinctDates_S2_sr, imgs, filt)

    # This function mosaics image acquired on the same day (same image swath)
    def mosaic_collection(img):
        orig = img

        # create a collection of the date-matching images
        col = ee.ImageCollection.fromImages(img.get('date_match')) 

        # extract collection properties to assign to the mosaic
        time_start = col.aggregate_min('system:time_start')
        time_end = col.aggregate_max('system:time_end')
        index_list = col.aggregate_array('system:index')
        index_list = index_list.join(',')
        scene_count = col.size()

        # get the unified geometry of the collection (outer boundary)
        col_geo = col.geometry().dissolve()

        # clip the mosaic to set a geometry to it
        mosaic = col.mosaic().clip(col_geo).copyProperties(img, ["system:time_start", "system:index", "date", "month", "SENSING_ORBIT_NUMBER", "PROCESSING_BASELINE", "SPACECRAFT_NAME", "MEAN_SOLAR_ZENITH_ANGLE", "MEAN_SOLAR_AZIMUTH_ANGLE"])

        # set the extracted properties to the mosaic
        mosaic = mosaic.set({'system:time_start': time_start,
                             'system:time_end': time_end,
                             'index_list': index_list,
                             'scene_count': scene_count})

        # # (don't think this is needed) reset the projection to epsg:32632 as mosaic changes it to epsg:4326 (otherwise the registration fails)
        # mosaic = ee.Image(mosaic).setDefaultProjection('epsg:32613', null, 10) #'epsg:32636'

        return mosaic
  
    imgs = ee.ImageCollection(joinCol_S2_sr.map(mosaic_collection))
    return imgs


def s2_masking(img):
    """ Function for wrapping the entire process to be applied across collection
    """
    img = sentinelCloudScore(img)
    cloudMask = (img.select(['cloudScore']).gt(cloudThresh)
                    .focal_min(contractPixels).focal_max(dilatePixels))
    img = projectShadows(cloudMask, img, cloudHeights)

    return img


def bustClouds(img):
    """ Mask clouds in Sentinel-2 image using cloudScoring method"    
    """
    img = sentinelCloudScore(img)
    img = img.updateMask(img.select(['cloudScore']).gt(cloudThresh).focal_min(contractPixels).focal_max(dilatePixels).Not())
    return img


def s2_bustclouds_tdom(imgs):
    """ Function to mask clouds with cloudScore method and shadows with TDOM
    
    Convert a sentinel-2 collection into a collection of daily mosaics then
    mask clouds and cloud shadows.
    
    Author: Ian Housman
    """
    imgs = dailyMosaics(imgs)
    imgs = imgs.map(bustClouds)
    imgs = simpleTDOM2(imgs)
    return imgs


###############################################################################
#------------- Sentinel-2 spectral index functions ----------------------------

def sri(i):
    """
    Simple Ratio Index (Jordan 1969) https://doi.org/10.2307/1936256
    """
    return (i.expression("nir / red", 
                        {'nir' : i.select('nir'),
                         'red' : i.select('red')
                        })
            .select([0], ['sri']))

def ndvi(i):
    """
    Normalized Difference Vegetation Index (Rouse et al, 1974) Monitoring vegetation systems in the Great Plains with ERTS, in: Proceedings. Presented at the 3rd Earth Resource Technology Satellite (ERTS) Symposium, pp. 48–62
    """
    return (i.normalizedDifference(["nir", "red"])
             .select([0], ["ndvi"]))


def evi(i):
    """
    Enhanced Vegetation Index (Liu and Huete, 1995) https://doi.org/10.1109/TGRS.1995.8746027
    """
    return (i.expression("2.5* ((nir - red) / (nir + 6.0 * red - 7.5 * blue + 1.))",
                      {'nir':i.select('nir'),
                       'red':i.select('red'),
                       'blue':i.select('blue')
                      })
           .select([0], ["evi"]))


def savi(i, L=0.5):
    """
    Soil Adjusted Vegetation Index (Huete 1988) https://doi.org/10.1016/0034-4257(88)90106-X
    """
    return (i.expression("((nir - red) / (nir + red + L)) * (1 + L)",
                      {'nir':i.select('nir'),
                       'red':i.select('red'),
                       'L':L
                      })
          .select([0], ["savi"]))


def msavi(i):
    """
    Modified Soil Adjusted Vegetation Index (Qi et al., 1994) https://doi.org/10.1016/0034-4257(94)90134-1
    """
    return (i.expression("(2. * nir + 1. - ((2. * nir + 1.)**2 - 8. * (nir - red))**0.5) / 2.",
                      {'nir':i.select('nir'),
                       'red':i.select('red')
                      })
           .select([0], ["msavi"]))


def satvi(i, L=0.5):
    """
    Soil Adjusted Total Vegetation Index (Marsett et al. 2006) https://doi.org/10.2111/05-201R.1
    Hill 2013 RSE adapted it for S2 and found it corresponded to a gradient of increasing tree cover.
    """
    return (i.expression("((swir1 - red) / (swir1 + red + L)) * (1 + L) - (swir2 / 2)",
                    {'swir1':i.select('swir1'),
                     'red':i.select('red'),
                     'swir2':i.select('swir2'),
                     'L':L
                      })
          .select([0], ["satvi"]))


def ndmi(i):
    """
    Normalized Difference Moisture Index (Gao 1996) https://doi.org/10.1016/S0034-4257(96)00067-3
    """
    return (i.normalizedDifference(["nir", "swir1"])
           .select([0], ["ndmi"]))


def nbr(i):
    """
    Normalized Burn Ratio (Key and Benson 2006) No. RMRS-GTR-164-CD
    """
    return (i.normalizedDifference(["nir", "swir2"])
           .select([0], ["nbr"]))


def nbr2(i):
    """
    Normalized Burn Ratio 2 (Product Guide: Landsat Surface Reflectance-Derived Spectral Indices; 3.6 Version; USGS)
    I can't find a basis in the scientific literature.
    """
    return (i.normalizedDifference(["swir1", "swir2"])
          .select([0], ["nbr2"]))


def cri1(i):
    """
    Carotenoid Reflectance Index 1 (Gitelson et al 2002). https://doi.org/10.1016/S0034-4257(01)00289-9 
    Hill 2013 RSE found it discriminated treed landscape from others in texas savannahs
    """
    return (i.expression("(1 / blue) - (1 / green)",
                      {'green':i.select('green'),
                       'blue':i.select('blue')
                      })
          .select([0], ["cri1"]))

# --- Red-edge indices---# 
    
def s2rep(i):
    """ Red-edge Position Index (Frampton et al. 2013; https://doi.org/10.1016/j.isprsjprs.2013.04.007)
    Indicator of crop N and growth status
    """
    return (i.expression(".705 + 0.035 * (0.5 * (re3 + red) - re1) / (re2 - re1)",
                      {'re1':i.select('re1'),
                       're2':i.select('re2'),
                       're3':i.select('re3'),
                       'red':i.select('red')
                      })
           .select([0], ["s2rep"]))


def ndvi705(i):
    """ Red-edge Normalized Difference Vegetation Index (Gitelson and Merzlyak 1994) https://doi.org/10.1016/S0176-1617(11)81633-0
    Proposed for assessing chlorophyll.
    """
    return (i.normalizedDifference(["re2", "re1"])
             .select([0], ["ndvi705"]))


def ndi45(i):
    """NDI45. Delegido 2011. https://doi.org/10.3390/s110707063
    Frampton 2013 found correlation to crop LAI.
    """
    return (i.normalizedDifference(["re1", "red"])
             .select([0], ["ndi45"]))
    

def mndvi705(i):
    """ Modified Red Edge Normalized Difference Vegetation Index (Sims et al. 2002) https://doi.org/10.1016/S0034-4257(02)00010-X
    Had high correlation with leaf chlorophyll content.    
    """
    return (i.expression("(re2 - re1) / (re2 + re1 - (2 * cb))",
                      {'cb':i.select('cb'),
                       're1':i.select('re1'), 
                       're2':i.select('re2')
                      })
          .select([0], ["mndvi705"]))
    
    
def tcari1(i):
    # Transformed Chlorophyll Absorption in Reflectance Index (Haboudane et al. 2002) https://doi.org/10.1016/S0034-4257(02)00018-4
    return (i.expression("3 * ((re1 - red) - 0.2 * (re1 - green) * (re1 / red))",
                      {'green':i.select('green'),
                       'red':i.select('red'), 
                       're1':i.select('re1')
                      })
          .select([0], ["tcari1"]))
    
    
def tcari2(i):
    # TCARI2 (Wu et al. 2008) https://doi.org/10.1016/j.agrformet.2008.03.005
    # Korhonen 2017 found it highly correlated with canopy cover.
    return (i.expression("3 * ((re2 - re1) - 0.2 * (re2 - green) * (re2 / re1))",
                      {'green':i.select('green'),
                       're1':i.select('re1'), 
                       're2':i.select('re2')
                      })
          .select([0], ["tcari2"]))

def ireci(i):
    """Inverted Red-Edge Chlorophyll Index (Frampton et al. 2013; https://doi.org/10.1016/j.isprsjprs.2013.04.007). 
    Frampton 2013 found strong correlation with crop LAI for S2.
    """
    return (i.expression("(re3 - red)/(re1 / re2)",
                      {'re1':i.select('re1'),
                       're2':i.select('re2'),
                       're3':i.select('re3'),
                       'red':i.select('red')
                      })
          .select([0], ["ireci"]))


def ari1(i):
    """
    Anthocyanin Reflectance Index (Gitelson et al 2001). https://doi.org/10.1562/0031-8655(2001)0740038OPANEO2.0.CO2
    Hill 2013 RSE found it discriminated treed landscape from others in texas savannahs
    """
    return (i.expression("(1 / green) - (1 / re1)",
                      {'green':i.select('green'),
                       're1':i.select('re1')
                      })
          .select([0], ["ari1"]))


def specixs(image, ixlist='all', exclude=None):
    """ Get Sentinel-2 spectral indices

    Parameters
    ----------
    image: ee.Image
        Sentinel-2 image with bands labeled by color rather than band number with rename_img.

    indices: ee.List
        List of desired band indices to return in a new image.

     Returns
     -------
     ee.Image
         Image with given band indices.

    Author
    ------
    Steven Filippelli

    Date Modified: 2017-06-19
    """

    ixbands = ee.Dictionary({"sri" : sri(image),
                             "ndvi" : ndvi(image),
                             "evi" : evi(image),
                             "savi" : savi(image),
                             "msavi" : msavi(image),
                             "satvi" : satvi(image),
                             "ndmi" : ndmi(image),
                             "nbr" : nbr(image),
                             "nbr2" : nbr2(image),
                             "cri1" : cri1(image),
                             "s2rep": s2rep(image),
                             "ndi45": ndi45(image),
                             "ndvi705" : ndvi705(image),
#                              "mndvi705" : mndvi705(image),
                             "tcari1" : tcari1(image),
                             "tcari2" : tcari2(image),
                             "ireci" : ireci(image),
                             "ari1" : ari1(image)
                            })
    
    if ixlist=='all':
        ixlist = ixbands.keys()
        
    if exclude:
        ixlist = ixlist.removeAll(exclude)
    
    ixbands = ixbands.values(ixlist)

    # convert list of images to multiband image
    empty = ee.Image().select()

    ixbands = ee.Image(ixbands.iterate(lambda image, result: 
                                          ee.Image(result).addBands(image)
                                     , empty))
                              
    return ixbands

