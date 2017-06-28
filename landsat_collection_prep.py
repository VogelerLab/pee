# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 10:35:29 2017

@author: filip107
"""
import ee
import ee.mapclient
ee.Initialize()


####################################
# Filter collection by time periods and reduce collection to image.
def filter_reduce_seasons(images,year):
    """
    Input
    images: ee.ImageCollection.
    year: integer
    """
    
    #----------- Filter and reduce collection---------------------------------------------
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
    
    
    #----------- Rename bands ----------------------------------------------------------
    def appendToBandNames(i, string):
        names = i.bandNames()
        names = names.map(lambda name: ee.String(name).cat(string))
        return i.rename(names)
    
    # Rename annual bands
    annual_med = appendToBandNames(annual_med, "_annual_med")
    
    # Rename season bands
    spring_med = appendToBandNames(spring_med, "_spring_med")
    summer_med = appendToBandNames(summer_med, "_summer_med")
    fall_med = appendToBandNames(fall_med, "_fall_med")
    winter_med = appendToBandNames(winter_med, "_winter_med")
    
    # Rename season difference bands
    spring_summer_med = appendToBandNames(spring_summer_med, "_spring_summer_med")
    spring_fall_med = appendToBandNames(spring_fall_med, "_spring_fall_med")
    spring_winter_med = appendToBandNames(spring_winter_med, "_spring_winter_med")
    summer_fall_med = appendToBandNames(summer_fall_med, "_summer_fall_med")
    summer_winter_med = appendToBandNames(summer_winter_med, "_summer_winter_med")
    fall_winter_med = appendToBandNames(fall_winter_med, "_fall_winter_med")
    
    #---------- Combine images --------------------------------------------------------
    image_list = ee.List([annual_med, spring_med, summer_med, fall_med, winter_med,
                        spring_summer_med, spring_fall_med, spring_winter_med,
                        summer_fall_med, summer_winter_med, fall_winter_med])
  
    empty = ee.Image().select()
  
    merged = ee.Image(image_list.iterate(lambda image, result: 
                                         ee.Image(result).addBands(image)
                                        ,empty))
    
    return merged