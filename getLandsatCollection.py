import ee
import ee.mapclient
ee.Initialize()

#######################################/
#Function for acquiring Landsat SR image collection
#
# Remote Sensing Applications Center, USFS
# Main code body written by: Ian Housman, Karis Tenneson, and Carson Stam
# Adapted to snippet by Joshua Goldstein
# updated 8-Sept-16 by Joshua Goldstein
# Updated by Steven Filippelli 10-May-2017 to remove cloud scoring and use SR product instead of TOA

def getLandsatCollection(studyArea,startDate,endDate,startJulian,endJulian):
    """ 
    DESCRIPTION: Get collection of all landsat images that intersect the given feature,
                 feature collection, or geometry for the given time frame.

    User inputs:
    studyArea: Study area as feature, feature collection, or geometry

    startYear: First year to include imagery from

    endYear: Last year to include imagery from

    startJulian: Starting Julian Date- Supports wrapping for tropics and 
        southern hemisphere

    endJulian: Ending Julian date- Supports wrapping for tropics and 
        southern hemisphere
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