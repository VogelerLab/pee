# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:54:31 2017

@author: filip107
"""
import ee
import ee.mapclient
ee.Initialize()
import numpy as np

def landsat_harmonics(collection, band, peak_norm=None, rgb=False):
    """
    collection: ee.ImageCollection
        Earth engine image collection of a Landsat time series.
        
    band: str
        Name of band to select from each image in "images" for calculating
        harmonic phase and amplitude over time.
        
    rgb: bool
        Add band to returned image to display phase as hue and amplitude as
        brightness.
    
    Returns: ee.Image
        Earth engine image with harmonic phase and amplitude of the selected
        band for the given image collection
        
    Author: Nick Clinton. Adapted by Steven Filippeli.
    Date: 2017-06-14
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
    
    
    
