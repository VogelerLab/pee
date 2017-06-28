# -*- coding: utf-8 -*-

import numpy as np
import ee
import ee.mapclient
ee.Initialize()


def get_landsat_collection(studyArea, startDate, endDate, startJulian, endJulian):
    """ Return a Landsat collection.

    Get collection of all landsat images that intersect the given feature,
    feature collection, or geometry for the given time frame.

    Parameters
    ----------
    studyArea: ee.FeatureCollection, ee.Feature or ee.Geometry
        Area over which to retrieve Landsat images.

    startYear: int
        First year to include imagery from

    endYear: int
        Last year to include imagery from

    startJulian: int
        Starting Julian Date- Supports wrapping for tropics and 
        southern hemisphere

    endJulian: int
        Ending Julian date- Supports wrapping for tropics and 
        southern hemisphere
        
    Returns
    -------
    ee.ImageCollection
        Collection of all Landsat images intersecting study area within the 
        given time frame.
        
        
    Authors
    -------
    Main code body written by: Ian Housman, Karis Tenneson, and Carson Stam
    Adapted to snippet by Joshua Goldstein
    Updated 8-Sept-16 by Joshua Goldstein
    Updated by Steven Filippelli 10-May-2017 to remove cloud scoring and use
        SR product instead of TOA
    
    Date Modified: 2017-06-17
    """
    sensorBandDictLandsatSR =ee.Dictionary({'L8' : ee.List([1,2,3,4,5,6,7,8]), # no coastal for L8
                        'L7' : ee.List([0,1,2,3,4,5,6,7]),
                        'L5' : ee.List([0,1,2,3,4,5,6,7]),
                        'L4' : ee.List([0,1,2,3,4,5,6,7])})

    bandNamesLandsatSR = ee.List(['blue','green','red','nir','swir1','swir2', 'cfmask', 'cfmask_conf'])

    l4SRs = (ee.ImageCollection('LANDSAT/LT4_SR')
      .filterDate(startDate,endDate)
      .filter(ee.Filter.calendarRange(startJulian,endJulian))
      .filterBounds(studyArea)
      .select(sensorBandDictLandsatSR.get('L4'),bandNamesLandsatSR))


    # TODO: clip L5 edges
    l5SRs = (ee.ImageCollection('LANDSAT/LT5_SR')
            .filterDate(startDate,endDate)
            .filter(ee.Filter.calendarRange(startJulian,endJulian))
            .filterBounds(studyArea)
            .select(sensorBandDictLandsatSR.get('L5'),bandNamesLandsatSR))

    l8SRs = (ee.ImageCollection('LANDSAT/LC8_SR')
            .filterDate(startDate,endDate)
            .filter(ee.Filter.calendarRange(startJulian,endJulian))
            .filterBounds(studyArea)
            .select(sensorBandDictLandsatSR.get('L8'),bandNamesLandsatSR))

    l7SRs = (ee.ImageCollection('LANDSAT/LE7_SR')
            .filterDate(startDate,endDate)
            .filter(ee.Filter.calendarRange(startJulian,endJulian))
            .filterBounds(studyArea)
            .select(sensorBandDictLandsatSR.get('L7'),bandNamesLandsatSR))

    # get merged image collection without messing up the system:index
    clist = ee.List([l4SRs, l5SRs, l8SRs, l7SRs])
    def get_images(collection):
        return (ee.ImageCollection(collection)
                .toList(ee.ImageCollection(collection).size().max(1)))
        
    ilist = clist.map(get_images).flatten()
    
    out = ee.ImageCollection(ilist)

    return out


def filter_reduce_seasons(images,year):
    """ Filter collection and reduce an image collection by year and seasons.
    
    Parameters
    ----------
    images: ee.ImageCollection
    
    year: int
    
    Returns
    -------
    ee.Image
        Image with bands for median of 'year', seasons within 'year' for
        each of the original bands in 'images', and differences between
        seasons.
        
    Authors
    -------
    Steven Filippelli
    """
    
    # Filter and reduce by year and seasons
    annual_med = images.median()
    spring_med = images.filterDate(str(year)+'-03-20', str(year)+'-06-20').median()
    summer_med = images.filterDate(str(year)+'-06-21', str(year)+'-09-22').median()
    fall_med = images.filterDate(str(year)+'-09-23', str(year)+'-12-20').median()
    #TODO: Use continuous winter(s). Winter 2010-11 and 11-12. Rather than splitting winter by other seasons.
    winter_med = (images.filter(ee.Filter.Or(ee.Filter.date(str(year)+'-01-01', str(year)+'-03-19'),
                                             ee.Filter.date(str(year)+'-12-21', str(year)+'-12-31')))
                        .median())

    # Get difference between seasons
    spring_summer_med = spring_med.subtract(summer_med)
    spring_fall_med = spring_med.subtract(fall_med)
    spring_winter_med = spring_med.subtract(winter_med)
    summer_fall_med = summer_med.subtract(fall_med)
    summer_winter_med = summer_med.subtract(winter_med)
    fall_winter_med = fall_med.subtract(winter_med)


    # Rename bands
    def appendToBandNames(i, string):
        names = i.bandNames()
        names = names.map(lambda name: ee.String(name).cat(string))
        return i.rename(names)
    
    annual_med = appendToBandNames(annual_med, "_annual_med")
    spring_med = appendToBandNames(spring_med, "_spring_med")
    summer_med = appendToBandNames(summer_med, "_summer_med")
    fall_med = appendToBandNames(fall_med, "_fall_med")
    winter_med = appendToBandNames(winter_med, "_winter_med")
    
    spring_summer_med = appendToBandNames(spring_summer_med, "_spring_summer_med")
    spring_fall_med = appendToBandNames(spring_fall_med, "_spring_fall_med")
    spring_winter_med = appendToBandNames(spring_winter_med, "_spring_winter_med")
    summer_fall_med = appendToBandNames(summer_fall_med, "_summer_fall_med")
    summer_winter_med = appendToBandNames(summer_winter_med, "_summer_winter_med")
    fall_winter_med = appendToBandNames(fall_winter_med, "_fall_winter_med")
    
    # Combine bands into a single image
    image_list = ee.List([annual_med, spring_med, summer_med, fall_med, winter_med,
                        spring_summer_med, spring_fall_med, spring_winter_med,
                        summer_fall_med, summer_winter_med, fall_winter_med])
  
    empty = ee.Image().select()
  
    merged = ee.Image(image_list.iterate(lambda image, result: 
                                         ee.Image(result).addBands(image)
                                        ,empty))
    
    return merged


def landsat_harmonics(collection, band, peak_norm=None, rgb=False):
    """ Return amplitude and phase image for a band in an image collection.
    
    Parameters
    ----------
    collection: ee.ImageCollection
        Earth engine image collection of a Landsat time series.
        
    band: str
        Name of band to select from each image in "images" for calculating
        harmonic phase and amplitude over time.
        
    rgb: bool
        Add band to returned image to display phase as hue and amplitude as
        brightness.
    
    Returns
    -------
    ee.Image
        Earth engine image with harmonic phase and amplitude of the selected
        band for the given image collection
        
    Authors
    -------
    Nick Clinton. Adapted by Steven Filippeli.
    
    Date Modified: 2017-06-14
    """
    def get_time(image):
        # Compute time in fractional years since the epoch.
        timeField = 'system:time_start'
        date = ee.Date(image.get(timeField))
        years = date.difference(ee.Date('1970-01-01'), 'year')
        return image.addBands(ee.Image(years).rename(['t']).float())

    def get_constant(image):
        return image.addBands(ee.Image.constant(1))

    # Get images
    images = collection.map(get_time).map(get_constant)

    # Name of the dependent variable.
    dependent = ee.String(band)

    # Use these independent variables in the harmonic regression.
    harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin'])

    # Add harmonic terms as new image bands.
    def get_cos_sin(image):
        timeRadians = image.select('t').multiply(2 * np.pi)
        return (image
                .addBands(timeRadians.cos().rename(['cos']))
                .addBands(timeRadians.sin().rename(['sin'])))

    harmonicLandsat = images.map(get_cos_sin)

    # The output of the regression reduction is a 4x1 array image.
    harmonicTrend = (harmonicLandsat
      .select(harmonicIndependents.add(dependent))
      .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1)))

    # Turn the array image into a multi-band image of coefficients.
    harmonicTrendCoefficients = (harmonicTrend.select('coefficients')
      .arrayProject([0])
      .arrayFlatten([harmonicIndependents]))

    # Compute phase and amplitude.
    phase = (harmonicTrendCoefficients.select('cos')
                .atan2(harmonicTrendCoefficients
                       .select('sin')).rename([band+'_phase']))

    amplitude = (harmonicTrendCoefficients.select('cos').hypot(
                    harmonicTrendCoefficients.select('sin'))
                    .rename([band+'_amplitude']))

    harmonic = phase.addBands(amplitude)

    # Use the HSV to RGB transform to display phase and amplitude
    if rgb:
        rgb = phase.unitScale(-np.pi, np.pi).addBands(
                  amplitude.multiply(2.5)).addBands(
                  ee.Image(1)).hsvToRgb()

        harmonic = harmonic.addBands(rgb)

    return harmonic


###############################################################################
#-------------Landsat spectral index functions --------------------------------
def ndvi(i):
    return (i.normalizedDifference(["nir", "red"])
             .select([0], ["ndvi"]))


def evi(i):
    return (i.expression("2.5* ((nir - red) / (nir + 6.0 * red - 7.5 * blue + 1.))",
                      {'nir':i.select('nir'),
                       'red':i.select('red'),
                       'blue':i.select('blue')
                      })
           .select([0], ["evi"]))


def savi(i):
  # Assuming L=0.5
    return (i.expression("((nir - red) / (nir + red + 0.5)) * 1.5",
                      {'nir':i.select('nir'),
                       'red':i.select('red')
                      })
          .select([0], ["savi"]))


def msavi(i):
    return (i.expression("(2. * nir + 1. - ((2. * nir + 1.)**2 - 8. * (nir - red))**0.5) / 2.",
                      {'nir':i.select('nir'),
                       'red':i.select('red')
                      })
           .select([0], ["msavi"]))


def ndmi(i):
    return (i.normalizedDifference(["nir", "swir1"])
           .select([0], ["ndmi"]))


def nbr(i):
    return (i.normalizedDifference(["nir", "swir2"])
           .select([0], ["nbr"]))


def nbr2(i):
    return (i.normalizedDifference(["swir1", "swir2"])
          .select([0], ["nbr2"]))


# ---Tasseled cap band coefficients--------
# Baig 2014 coeffs - TOA refl (http:#www.tandfonline.com/doi/pdf/10.1080/2150704X.2014.915434)
l8_tc_coeffs = [ee.Image([0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872]),
              ee.Image([-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608]),
              ee.Image([ 0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559]),
              ee.Image([-0.8239, 0.0849, 0.4396, -0.058, 0.2013, -0.2773]),
              ee.Image([-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085]),
              ee.Image([0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252])]
# Huang 2002 coeffs - TOA refl (http:#landcover.usgs.gov/pdf/tasseled.pdf)
l7_tc_coeffs = [ee.Image([0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596]),
              ee.Image([-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630]),
              ee.Image([0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388]),
              ee.Image([0.0805, -0.0498, 0.1950, -0.1327, 0.5752, -0.7775]),
              ee.Image([-0.7252, -0.0202, 0.6683, 0.0631, -0.1494, -0.0274]),
              ee.Image([0.4000, -0.8172, 0.3832, 0.0602, -0.1095, 0.0985])]
  
# Crist 1985 coeffs - TOA refl (http:#www.gis.usu.edu/~doug/RS5750/assign/OLD/RSE(17)-301.pdf)
l5_tc_coeffs = [ee.Image([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]),
              ee.Image([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]),
              ee.Image([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]),
              ee.Image([-0.2117, -0.0284, 0.1302, -0.1007, 0.6529, -0.7078]),
              ee.Image([-0.8669, -0.1835, 0.3856, 0.0408, -0.1132, 0.2272]),
              ee.Image([0.3677, -0.8200, 0.4354, 0.0518, -0.0066, -0.0104])]
            
tc_coeff_dict =ee.Dictionary({"LANDSAT_5": l5_tc_coeffs,
                              "LANDSAT_7": l7_tc_coeffs,
                              "LANDSAT_8": l8_tc_coeffs})
tc_names = ['brightness','greenness', 'wetness', 'tct4','tct5','tct6']


def tasseled_cap(image):
    """ Function for computing tasseled cap transform

    Parameters
    ----------
    image: ee.Image
        Image to perform tasseled cap on. Bands must be named by color not num.

    Returns
    -------
    ee.Image
        Image with tasseled cap bands.
        
    Authors
    -------
    Ian Housman, adapted by Steven Filippelli
        
    TODO: tcap bands lose dimensions. Fix this.
    """
    
    image = ee.Image(image).select(["blue", "green", "red", "nir", "swir1", "swir2"])
    # Nested function to do the multiplication and addition
    def mult_sum(matrix):
        return image.multiply(matrix).reduce(ee.call("Reducer.sum")).float()

    # Find the coeffs
    coeffs = ee.List(tc_coeff_dict.get(image.get('satellite')))
                  
    # Compute the tc transformation and give it useful names
    tco = coeffs.map(mult_sum)
    empty = ee.Image().select()
    tco = ee.Image(tco.iterate(lambda image, result: 
                              ee.Image(result).addBands(image),
                              empty))
  
    tco = tco.select([0,1,2,3,4,5], tc_names)
    
    return tco


def greenness(i):
    return tasseled_cap(i).select("greenness")


def wetness(i):
    return tasseled_cap(i).select("wetness")


def brightness(i):
    return tasseled_cap(i).select("brightness")


def landsat_rename_bands(image):
    # Function for rename the bands of TM, ETM+, and OLI to aid
    # easy interpretation and calculation of spectral indices
    band_dicts = ee.Dictionary({
    'LANDSAT_8': {'B1':'coastal', 'B2':'blue', 'B3':'green', 'B4':'red', 'B5':'nir', 'B6':'swir1', 'B7':'swir2'},
    'LANDSAT_7': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'},
    'LANDSAT_5': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'},
    'LANDSAT_4': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'}
    })
    band_dict = ee.Dictionary(band_dicts.get(image.get('satellite')))

    def swap_name(b):
        return ee.Algorithms.If(band_dict.contains(b), band_dict.get(b), b)
  
    new_bands = image.bandNames().map(swap_name)
    return image.rename(new_bands)


def landsat_spectral_indices(image, ixlist):
    """ Get landsat spectral indices

    Parameters
    ----------
    image: ee.Image
        Landsat 5-8 image with bands labeled by color rather than band number.
        An image collection with landsat bands relabeled can be obtained from
        the "get_landsat_sr_collection" script.

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

    # check if bands were already renamed, and if not then rename them
    image = landsat_rename_bands(image)

    ixbands = ee.Dictionary(
            {"ndvi":ndvi(image),
             "evi":evi(image),
             "savi":savi(image),
             "msavi":msavi(image),
             "ndmi":ndmi(image),
             "nbr":nbr(image),
             "nbr2":nbr2(image),
             "tasseled_cap":tasseled_cap(image),
             "greenness":greenness(image),
             "wetness":wetness(image),
             "brightness":brightness(image)
            })
    
    ixbands = ixbands.values(ixlist)

    # convert list of images to multiband image
    empty = ee.Image().select()

    ixbands = ee.Image(ixbands.iterate(lambda image, result: 
                                          ee.Image(result).addBands(image)
                                     , empty))
                              
    return ixbands
