# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:40:56 2017

@author: filip107
"""

# Scripts related to processing a time series which aren't specific to any
# particular sensor

import ee
import ee.mapclient
ee.Initialize()
import numpy as np


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