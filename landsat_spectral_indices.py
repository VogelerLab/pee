# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:22:07 2017

@author: filip107
"""
import ee
import ee.mapclient
ee.Initialize()

# Get landsat spectral indices 
# Inputs:
#   image: landsat 5-8 image with bands labeled by color rather than band number.
#          An image collection with landsat bands relabeled can be obtained from
#          the "get_landsat_sr_collection" script
#
#   indices: an ee.List of desired band indices to return in a new image.
#
# Returns: image with given band indices

# spectral index functions
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

#####################################
#Function for computing tasseled cap transform
#Author:  Ian Housman, adapted by Steven Filippelli
# 
#TODO: tcap bands lose dimensions. Fix this.

#Baig 2014 coeffs - TOA refl (http:#www.tandfonline.com/doi/pdf/10.1080/2150704X.2014.915434)
l8_tc_coeffs = [ee.Image([0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872]),
              ee.Image([-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608]),
              ee.Image([ 0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559]),
              ee.Image([-0.8239, 0.0849, 0.4396, -0.058, 0.2013, -0.2773]),
              ee.Image([-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085]),
              ee.Image([0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252])]
#Huang 2002 coeffs - TOA refl (http:#landcover.usgs.gov/pdf/tasseled.pdf)
l7_tc_coeffs = [ee.Image([0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596]),
              ee.Image([-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630]),
              ee.Image([0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388]),
              ee.Image([0.0805, -0.0498, 0.1950, -0.1327, 0.5752, -0.7775]),
              ee.Image([-0.7252, -0.0202, 0.6683, 0.0631, -0.1494, -0.0274]),
              ee.Image([0.4000, -0.8172, 0.3832, 0.0602, -0.1095, 0.0985])]
  
#Crist 1985 coeffs - TOA refl (http:#www.gis.usu.edu/~doug/RS5750/assign/OLD/RSE(17)-301.pdf)
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
  image = ee.Image(image).select(["blue", "green", "red", "nir", "swir1", "swir2"])
  #Nested function to do the multiplication and addition
  def mult_sum(matrix):
    return image.multiply(matrix).reduce(ee.call("Reducer.sum")).float()

  #Find the coeffs
  coeffs = ee.List(tc_coeff_dict.get(image.get('satellite')))
#  print(coeffs.getInfo())
                  
  #Compute the tc transformation and give it useful names
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


# get spectral indices from given list
def landsat_spectral_indices(image, ixlist):  
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