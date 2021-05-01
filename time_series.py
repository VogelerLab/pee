# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:40:56 2017

@author: filip107
"""

# Scripts related to processing a time series which aren't specific to any
# particular sensor

# TODO: Get rid of or simplify reduce seasons functions

import ee
# import ee.mapclient
ee.Initialize()
import numpy as np


def appendToBandNames(img, string):
    """ Append a "string" to each band name for the image ("img")
    """
    bands = img.bandNames()
    bands = bands.map(lambda b: ee.String(b).cat(string))
    return i.rename(bands)


def datetime_band(img, unit='day', inUnit='year', base=1):
    """ Return a band that contains datetime info based on system:time_start
    
    img: ee.Image
        Image to append datetime band to
        
    unit: str
        datetime unit for returned band passed to ee.Date.getRelative
        
    inUnit: str
        datetime unit that 'unit' is relative to. Passed to ee.Date.getRelative
    
    base: int
        base number to add to the value. getRelative() is 0-based, but a base of 1
        is likely more commonly desired (i.e. 1st day of year, not 0th day of year).
    """
    band = ee.Image(img.date().getRelative(unit, inUnit).add(base)).rename('date').toInt16() # todo: consider rename as unit+'_of_'+inUnit
    return band


def quality_mosaic(imgs, band, reducer=ee.Reducer.min, date_band=False, unit='day', inUnit='year'):
    """ Mosaic according to a reducer applied to a given band and pass corresponding bands
    for the selected pixel. This is equivalent to qualityMosaic but allows use 
    of the minimum or maximum reducer. 
    
    TODO: Update to use other reducers such as percentile or median. This 
    would probably require use of conversion to an array and sorting which will
    be less efficent.
    
    TODO: Maybe return the image in the original band order
    
    imgs: ee.ImageCollection
        A time series of images with the same bands
        
    band: str
        Name of band to apply reducer to
        
    reducer: ee.Reducer
        A reducer to apply to band for selecting a pixel. Currently only
        ee.Reducer.min and ee.Reducer.max will work.
    
    date_band: bool
        Optionally return a band with the datetime of the image selected from
        medoid calculation based on system:time_start. Default is day of year.
        
    unit: str
        Datetime unit to get for for date_band passed to ee.Date.getRelative
    
    inUnit:str
        Relative unit passed to ee.Date.getRelative.
        
    returns: ee.Image
        Mosaic containing all original bands and an optional 'date' band with 
        a pixel from the time series selected according to the reducer applied
        to the given band.
    
    """
#     ## getInfo version (which uses less memory?!)
#     # sort band as first 
#     bands = imgs.first().bandNames().getInfo()
#     bands = [band] + [b for b in bands if b!=band]
#     imgs = imgs.select(bands)
    
#    # Add optional date band
#     if date_band:
#         bands = bands.add('date')
#         imgs = imgs.map(lambda i: i.addBands(datetime_band(i, unit, inUnit)))
    
#     # select the pixel according to the reducer and band
#     img = (ee.ImageCollection(imgs)
#                 .reduce(reducer(len(bands)))
#                 .select(list(range(len(bands))), bands)
#            )
    
    # ee version (no getInfo)
    # sort band as first
    bands = imgs.first().bandNames()
    bands = ee.List.cat([band], bands.remove(band))
    imgs = imgs.select(bands)
    
   # Add optional date band
    if date_band:
        bands = bands.add('date')
        imgs = imgs.map(lambda i: i.addBands(datetime_band(i, unit, inUnit)))
    
    # select the pixel according to the reducer and band
    img = (ee.ImageCollection(imgs)
                .reduce(reducer(bands.length()))
                .select(ee.List.sequence(0, bands.length().subtract(1)), bands)
           )
    
    return img


def medoid(imgs, med_bands=None, date_band=False, unit='day', inUnit='year'):
    """ Calculate the medoid for an image collection
    
    imgs: ee.ImageCollection
        A time series of images with the same bands
        
    med_bands: list
        List of band names to use for calculating the medoid
    
    date_band: bool
        Optionally return a band with the datetime of the image selected from
        medoid calculation based on system:time_start. Default is day of year.
        
    unit: str
        Datetime unit to get for for date_band passed to ee.Date.getRelative
    
    inUnit:str
        Relative unit passed to ee.Date.getRelative.
        
    returns: ee.Image
        Medoid mosaic containing all original bands, including those not used
        for medoid calculation, and an optional 'date' band.
    """
    bands = imgs.first().bandNames()
    # get median of selected bands
    if med_bands is None:
        med_bands = bands
    
    median = imgs.select(med_bands).median()
    
     # Add band with sum of squared differences from the median
    def diff_from_median(img, med_bands):
        """Sum of squared differences from the median"""
        diff = ee.Image(img).select(med_bands).subtract(median).pow(ee.Image.constant(2))      
        img = diff.reduce('sum').addBands(img).copyProperties(img, img.propertyNames())
        return img
    
    med_dif = imgs.map(lambda i: diff_from_median(i, med_bands))
    
   # Add optional date band
    if date_band:
        bands = bands.add('date')
        med_dif = med_dif.map(lambda i: i.addBands(datetime_band(i, unit, inUnit)))
    
    # select the pixel with the min dif, then drop dif band and restore original band names
    img = (ee.ImageCollection(med_dif)
                .reduce(ee.Reducer.min(bands.length().add(1)))
                .select(ee.List.sequence(1, bands.length()), bands)
           )
    
    return img


def lt_flipbands(img, flip_bands=True):
    """ Orient values of given bands in preparation for LandTrendr segmentation.
    For LandTrendr bands must be oriented such that a disturbance would cause an
    INCREASE in value.
    
    img: ee.Image
        The image to have bands flipped
    
    flip_bands: bool or list
        List of band names to flip (i.e. multiply by -1)
        If True, flip bands in imgs that contain commonly flipped band names.
        Be careful that these common names are not a part of any band names that 
        should not to be flipped.
        
    returns: ee.Image
        Image with requested bands flipped
    """
    # TODO: May be faster to just multiply all bands image by a constant image
    #       of 1 and -1 depending on bands to flip.
    bands = img.bandNames()
    
    # check substring of bands for common bands to flip (e.g. fall_ndvi)
    if not isinstance(flip_bands, list):
        flip_bands = ee.List(["nir", "sri", "ndvi", "evi", "savi", "msavi", "ndmi", "nbr", "nbr2", "cri1", "satvi", "tcg", "tcw", "tca"])
        flip_bands = flip_bands.map(lambda fb: bands.filter(ee.Filter.stringContains(leftField='item', rightValue=fb))).flatten()
    
    flipped = img.select(flip_bands).multiply(-1).copyProperties(img, img.propertyNames())
    img = img.addBands(flipped, overwrite=True)
    
    return img


def lt_fitted(imgs, flip_bands=True, fit_band=None, maxSegments=6, **kwargs):
    """ Get LandTrendr fitted values for a collection of annual images
    
    imgs: ee.ImageCollection
        A collection of annual images with system:time_start properties that are a
        year apart. Bands must be oriented so a disturbance would cause an increase
        in value or they must flipped using flip_bands.
    
    flip_bands: bool or list
        If True reorient bands for LandTrendr using lt_flipbands, selecting
        bands which contain the acronym of commonly flipped bands. Provide a list 
        of band names to flip specific bands. Flipped bands are flipped back to 
        their original orientation after running LandTrendr.
    
    fit_band: str
        Name of the band to use for temporal segmentation. If None, fit each 
        band independently.
        
    maxSegments: int
        Maximum number of segments for LandTrendr
        
    **kwargs:
        Arguments passed to ee.Algorithms.LandTrendr
    """
    bands = imgs.first().bandNames()
    fitnames = bands.map(lambda b: ee.String(b).cat('_fit'))
    
#     imgs = imgs.map(lambda i: i.toInt16())
        
    if flip_bands:
        imgs = imgs.map(lambda i: lt_flipbands(i, flip_bands=flip_bands))
    
    if fit_band:
        # put copy of fitting band as first band
        imgs = imgs.map(lambda i: i.select([fit_band], [fit_band+'_ftv']).addBands(i))
        
        # run LandTrendr and rename bands
        ltimg = ee.Algorithms.TemporalSegmentation.LandTrendr(
            timeSeries=imgs, maxSegments=maxSegments, **kwargs
        )
        ltimg = ltimg.select(fitnames, bands)        
        
    else:
        # run LT and get fitted for each band
        def lt_fitband(band):
            """ run LandTrendr for a single band and get fitted"""
            band = ee.String(band)
            # Duplicate band to fit
            imgs_band = imgs.map(lambda i: i.select([band, band], [band.cat('_ftv'), band]))
            ltimg = ee.Algorithms.TemporalSegmentation.LandTrendr(
                timeSeries=imgs_band, maxSegments=maxSegments, **kwargs
            )
            ltimg = ltimg.select([band.cat('_fit')], [band])
            
            return ltimg            
        
        lt_list = bands.map(lt_fitband)
           
        # convert list of bands into a single image (ee.Image.cat doesn't work)
        def combineBands(image, result):
            return ee.Image(result).addBands(image)

        empty = ee.Image().select()
        ltimg = ee.Image(lt_list.iterate(combineBands, empty))

    # Convert band arrays back to image collection of years
    # TODO: Faster to first set year on fitted images then join to imgs to get properties?
    years = (imgs.aggregate_array('system:time_start')
            .map(lambda x: ee.Date(x).get('year')).distinct().sort()
        )
    y1 = ee.Number(years.get(0))
    
    def fitted_year(y):
        ydate = ee.Date.fromYMD(ee.Number(y), 1, 1)
        yfit = ltimg.arrayGet(ee.Number(y).subtract(y1))
        yimg = ee.Image(imgs.filterDate(ydate, ydate.advance(1, 'year')).first())
        yfit = yfit.copyProperties(yimg, yimg.propertyNames())
        return yfit
    
    imgs_fit = ee.ImageCollection(years.map(fitted_year))
    
    if flip_bands:
        # TODO: probably faster to flip band used for vertices calc and leave "fitted" copy unflipped rather than flipping both twice.
        imgs_fit = imgs_fit.map(lambda i: lt_flipbands(i, flip_bands=flip_bands))
    
    return imgs_fit
  

def annual_composites(aoi, starty, endy, startdoy, enddoy, coll_func, comp_func,
                      coll_kwargs={}, comp_kwargs={}, fill=False):
    """ Create an annual composite image collection using a function to prepare
    an image collection and a function to reduce a collection to an image.
    
    aoi: ee.Geometry
        Area of interest over which to grab images
    
    starty: int
        Start year of output image collection
        
    endy: int
        End year of output image collection (inclusive)
        
    startdoy: int
        Start day of year for filtering the base image collection
        
    enddoy: int
        End day of year for filtering the base image collection (inclusive).
        Wraps to next year if less than startdoy.
        
    coll_func: function
        Function for constructing an ee.ImageCollection that takes an aoi (ee.Geometry),
        start date (ee.Date), and end date (ee.Date) as it's first three arguments.
        
    comp_func: function
        Function that reduces an image collection to a single image taking the
        input image collection as it's first argument.
        
    coll_kwargs: dict
        Dictionary of additional keyword arguments to pass to coll_func
        
    comp_kwargs: dict
        Dictionary of additional keyword arguments to pass to comp_func
        
    fill: bool
        If True fill years between starty and endy with an empty image if there
        are no valid images returned by coll_func for a particular year
        
    Returns: ee.ImageCollection
        An image collection of annual compsite images with year and system:time_start
        properties based on the startdoy.
    """
    
    def year_comp(y):
        # TODO: maybe add fill here instead
        start = ee.Date.fromYMD(y,1,1).advance(startdoy-1, 'day')
        endy = y if enddoy > startdoy else ee.Number(y).add(1)
        end = ee.Date.fromYMD(endy,1,1).advance(enddoy-1, 'day')
        
        imgs = ee.ImageCollection(coll_func(aoi, start, end, **coll_kwargs))
        
        img = ee.Algorithms.If(imgs.size(), 
                               (ee.Image(comp_func(imgs, **comp_kwargs))
                                .set({'system:time_start':start.millis(), 
                                    'year':y
                                    })),
                               None        #ee.Image(0).updateMask(0)
                             )
        
        return img
                       
    years = ee.List.sequence(starty, endy)
    imgs = ee.ImageCollection(years.map(year_comp))
    
    if fill:
        # fill years with no valid images with an empty image (needed for LandTrendr fit across empty years)
        val_years = (imgs.aggregate_array('system:time_start')
            .map(lambda x: ee.Date(x).get('year')).distinct().sort()
        )
        missing = years.removeAll(val_years)
        names = imgs.first().bandNames()
        dtypes = imgs.first().bandTypes()
        empty = (ee.Image.constant(ee.List.repeat(0, names.size()))
                     .rename(names)
                     .cast(dtypes, names)
                     .mask(ee.Image(0)))
        fillers = missing.map(lambda y: empty.set({
            'system:time_start':ee.Date.fromYMD(y, 1, 1).advance(startdoy-1, 'day').millis(),
            'year':y})
                             )
        fillers = ee.ImageCollection(fillers)
        
        imgs = ee.ImageCollection.merge(imgs, fillers).sort('system:time_start')
    
    return imgs


def fit_harmonic(collection, band, peak_norm=None, rgb=False):
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

    harmonicImage = images.map(get_cos_sin)

    # The output of the regression reduction is a 4x1 array image.
    harmonicTrend = (harmonicImage
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
