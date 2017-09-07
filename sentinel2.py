# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:24:04 2017

@author: filip107
"""
import ee
import ee.mapclient
ee.Initialize()
from math import pi


#--------------- Cloud Masking ----------------------------------------------#

# CloudScore originally written by Matt Hancher and adapted for S2 data by Ian Housman
# User Params
cloudThresh = 20 # Ranges from 1-100.Lower value will mask more pixels out. Generally 10-30 works well with 20 being used most commonly 
cloudHeights = ee.List.sequence(200,10000,500) # Height of clouds to use to project cloud shadows
irSumThresh = 0.35 # Sum of IR bands to ☻include as shadows within TDOM and the shadow shift method (lower number masks out less)
dilatePixels = 2 # Pixels to dilate around clouds
contractPixels = 1 # Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors


def rescale(img, exp, thresholds):
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

    # Clouds are reasonably bright in the blue and cirrus bands.
    score = score.min(rescale(img, 'img.blue', [0.1, 0.5]))
    score = score.min(rescale(img, 'img.cb', [0.1, 0.3]))
    score = score.min(rescale(img, 'img.cb + img.cirrus', [0.15, 0.2]))

    # Clouds are reasonably bright in all visible bands.
    score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]))

    # Clouds are moist
    ndmi = img.normalizedDifference(['nir', 'swir1'])
    score = score.min(rescale(ndmi, 'img', [-0.1, 0.1]))

    # However, clouds are not snow.
    ndsi = img.normalizedDifference(['green', 'swir1'])
    score = score.min(rescale(ndsi, 'img', [0.8, 0.6]))

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


def simpleTDOM2(c):
  """ Cloud shadow masking with TDOM
  
  Author: Ian Housman
  """
  shadowSumBands = ['nir','swir1']
  irSumThresh = 0.4
  zShadowThresh = -1.2
  # Get some pixel-wise stats for the time series
  irStdDev = c.select(shadowSumBands).reduce(ee.Reducer.stdDev())
  irMean = c.select(shadowSumBands).mean()
  bandNames = ee.Image(c.first()).bandNames()
#  print('bandNames',bandNames)
  # Mask out dark dark outliers
  def dark_outliers(img):
    z = img.select(shadowSumBands).subtract(irMean).divide(irStdDev)
    irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum())
    m = z.lt(zShadowThresh).reduce(ee.Reducer.sum()).eq(2).And(irSum.lt(irSumThresh)).Not()
    return img.updateMask(img.mask().And(m))

  c = c.map(dark_outliers)

  return c.select(bandNames)


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


def s2_masking(img, shadows='TDOM'):
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
    return i.select("nir").divide(i.select("red"))

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

def re_ndvi(i):
    return (i.normalizedDifference(["nir", "re2"])
             .select([0], ["re_ndvi"]))
    
    
def re_pos(i):
    return (i.expression(".705 + 0.035 * (0.5 * (re3 + red) - re1) / (re2 - re1)",
                      {'re1':i.select('re1'),
                       're2':i.select('re2'),
                       're3':i.select('re3'),
                       'red':i.select('red')
                      })
           .select([0], ["re_pos"]))
    
def ndvi705(i):
    return (i.normalizedDifference(["re2", "re1"])
             .select([0], ["ndvi705"]))
    

def mndvi705(i):
    return (i.expression("(re2 - re1) / (re2 + re1 - (2 * cb))",
                      {'cb':i.select('cb'),
                       're1':i.select('re1'), 
                       're2':i.select('re2')
                      })
          .select([0], ["mndvi705"]))
    

def ndvig(i):
    return (i.expression("green * (nir - red) / (nir + red)",
                      {'green':i.select('green'),
                       'red':i.select('red'), 
                       'nir':i.select('nir')
                      })
          .select([0], ["ndvig"]))
    
    
def tcari1(i):
    return (i.expression("3 * ((re1 - red) - 0.2 * (re1 - green) * (re1 / red))",
                      {'green':i.select('green'),
                       'red':i.select('red'), 
                       're1':i.select('re1')
                      })
          .select([0], ["tcari1"]))
    
    
def tcari2(i):
    return (i.expression("3 * ((re2 - re1) - 0.2 * (re2 - green) * (re2 / re1))",
                      {'green':i.select('green'),
                       're1':i.select('re1'), 
                       're2':i.select('re2')
                      })
          .select([0], ["tcari2"]))


#def sentinel2_rename_bands(image):
#    # Function for rename the bands of TM, ETM+, and OLI to aid
#    # easy interpretation and calculation of spectral indices
#    band_dicts = ee.Dictionary({
#    'LANDSAT_8': {'B1':'coastal', 'B2':'blue', 'B3':'green', 'B4':'red', 'B5':'nir', 'B6':'swir1', 'B7':'swir2'},
#    'LANDSAT_7': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'},
#    'LANDSAT_5': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'},
#    'LANDSAT_4': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 'B5':'swir1', 'B7':'swir2'}
#    })
#    band_dict = ee.Dictionary(band_dicts.get(image.get('satellite')))
#
#    def swap_name(b):
#        return ee.Algorithms.If(band_dict.contains(b), band_dict.get(b), b)
#  
#    new_bands = image.bandNames().map(swap_name)
#    return image.rename(new_bands)


def sentinel2_spectral_indices(image, ixlist):
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
#    image = landsat_rename_bands(image)

    ixbands = ee.Dictionary({"ndvi" : ndvi(image),
                             "evi" : evi(image),
                             "savi" : savi(image),
                             "msavi" : msavi(image),
                             "ndmi" : ndmi(image),
                             "nbr" : nbr(image),
                             "nbr2" : nbr2(image),
                             "re_ndvi" : re_ndvi(image),
                             "re_pos" : re_pos(image),
                             "ndvi705" : ndvi705(image),
                             "mndvi705" : mndvi705(image),
                             "ndvig" : ndvig(image),
                             "tcari1" : tcari1(image),
                             "tcari2" : tcari2(image)
                            })
    
    ixbands = ixbands.values(ixlist)

    # convert list of images to multiband image
    empty = ee.Image().select()

    ixbands = ee.Image(ixbands.iterate(lambda image, result: 
                                          ee.Image(result).addBands(image)
                                     , empty))
                              
    return ixbands

