# -*- coding: utf-8 -*- 
# Functions for preparing Landsat data in GEE. 

# There are a few more sophisticated public libraries that may be worth using instead:
# https://github.com/fitoprincipe/geedatasets
# https://github.com/george-azzari/eetc/tree/master/gee_tools/datasources

# TODO: Add function to get a scene's snow and cloud cover from pqa band (see redcedar_landtrendr_v4)

import numpy as np
import ee
import time_series
ee.Initialize()


__srbands = ee.Dictionary({
    'LANDSAT_4':{'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'therm':['therm1'],
                 'atmos':['sr_atmos_opacity'],
                 'qa':['sr_cloud_qa', 'pixel_qa', 'radsat_qa'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']
                },
    'LANDSAT_5':{'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'therm':['therm1'],
                 'atmos':['sr_atmos_opacity'],
                 'qa':['sr_cloud_qa', 'pixel_qa', 'radsat_qa'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']
                },
    'LANDSAT_7':{'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'therm':['therm1'],
                 'atmos':['sr_atmos_opacity'],
                 'qa':['sr_cloud_qa', 'pixel_qa', 'radsat_qa'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']
                },
    'LANDSAT_8':{'opt':['cb','blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'therm':['therm1', 'therm2'],
                 'atmos':[],
                 'qa':['sr_aerosol', 'pixel_qa', 'radsat_qa'],
                 'all':['cb', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'therm1', 'therm2', 'sr_aerosol', 'pixel_qa', 'radsat_qa']
                }   
    })

__atmos_scale = ee.Dictionary({
    'LANDSAT_4':ee.Image(0.001),
    'LANDSAT_5':ee.Image(0.001),
    'LANDSAT_7':ee.Image(0.001),
    'LANDSAT_8':ee.Image().select([])
})


def swap_bandnames(image):
    """Function to rename the band without assuming all bands present or original
       band order.
    """
    band_dicts = ee.Dictionary({
    'LANDSAT_8': {'B1':'cb', 'B2':'blue', 'B3':'green', 'B4':'red', 'B5':'nir', 
                  'B6':'swir1', 'B7':'swir2', 'B10':'therm1', 'B11':'therm2'},
    'LANDSAT_7': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 
                  'B5':'swir1', 'B7':'swir2'},
    'LANDSAT_5': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 
                  'B5':'swir1', 'B7':'swir2'},
    'LANDSAT_4': {'B1':'blue', 'B2':'green', 'B3':'red', 'B4':'nir', 
                  'B5':'swir1', 'B7':'swir2'}
    })
    band_dict = ee.Dictionary(band_dicts.get(image.get('SATELLITE')))

    def swap_name(b):
        # TODO: try to use filter instead of if
        return ee.Algorithms.If(band_dict.contains(b), band_dict.get(b), b)
  
    new_bands = image.bandNames().map(swap_name)
    return image.rename(new_bands)


def sr_rename(img):
    """ Rename bands based on satellite, assuming all bands present in original order.
    """
    bnames = ee.Dictionary(__srbands.get(img.get('SATELLITE'))).get('all')
    return img.rename(bnames)


def sr_rescale(img):
    """ Rescale Landsat SR images from any satellite.
    """
    sat = img.get('SATELLITE')
    bdict = ee.Dictionary(__srbands.get(sat))
    opt = img.select(bdict.get('opt'))
    therm = img.select(bdict.get('therm'))
    atmos = img.select(bdict.get('atmos'))
    qa = img.select(bdict.get('qa'))
    
    # scale and cast to float (may reduce memory bottlenecks?)
    opt = opt.multiply(0.0001).toFloat()
    atmos = atmos.multiply(__atmos_scale.get(sat)).toFloat() #0.001 for L57, empty for L8
    therm = therm.multiply(0.1).toFloat()

    scaled = ee.Image(ee.Image.cat([opt, therm, qa, atmos])).copyProperties(img)
    scaled = scaled.set('system:time_start', img.get('system:time_start'))
    return scaled


def l8sr_harmonize(img):
    """Harmonize L8 to L5/7 following Roy 20?? for surface reflectance images
       Expects images in the original scale (i.e.0-10000) with all renamed bands
    """
    # RMA slopes and intercepts per band
    slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949])
#     itcp_coefs = [-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029] # img scaled 0-1
    itcp_coefs = [-95.0, -16.0, -22.0, -21.0, -30.0,  29.0] # img scaled 0-10,000
    itcp = ee.Image.constant(itcp_coefs)
#     # Least squares slopes and intercepts per band
#     slopes = ee.Image.constant([0.885, 0.9317, 0.9372, 0.8339, 0.8639, 0.9165])
#     itcp = ee.Image.constant([0.0183, 0.0123, 0.0123, 0.0448, 0.0306, 0.0116])

    # apply reverse of regression (Roy gives ETM to OLI)
    y = (img.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            .subtract(itcp).divide(slopes)
            .toInt16()
        )
    new = ee.Image.cat([img.select(['cb']),
                        y,
                        img.select('therm1', 'therm2', 'sr_aerosol', 'pixel_qa', 'radsat_qa')])
    
    img = (new.copyProperties(img)
             .set('system:time_start', img.get('system:time_start'))
          )
#     # SF: now treat as L5/7 for calculating TCAP
#     img = img.set('SATELLITE', 'LANDSAT_7')
    return img


def pqa_mask(img, fill=0, clear=1, water=0, shadow=0, snow=0, cloud=0, 
             cloud_conf=1, cirrus_conf=1, terrain=0):
    """Keep pixels which meet all the given critera, and mask out others.
    Values should match the LEDAPS and LaSRC product guide table with 
    0=='No' and 1=='Yes', None='accept either' and conf values: 0,1,2,3=='None', 'low', 'med', 'high'
    https://www.usgs.gov/media/files/land-surface-reflectance-code-lasrc-product-guide
    https://www.usgs.gov/media/files/landsat-4-7-surface-reflectance-code-ledaps-product-guide
    
    'Less than or equal to' operation applied to conf values. For example, 
    cloud_conf=2 means keep medium and low confidence clouds. 
    
    Defaults retain clear, non-fill pixels with no water, snow, clouds, terrain
    occlusion, and low confidence clouds and cirrus.
    
    L4-7 does not have cirrus confidence or terrain occlusion. 
    Warning: Terrain occlusion set to 1 will completely mask out all pixels in L4-7.
    
    returns: ee.Image
        Mask which can be applied using ee.Image.updateMask()
    """    
    pixqa = img.select('pixel_qa')
    
    mask = img.mask() # use existing mask or ee.Image(1)?
    if fill is not None:
        mask = mask.And(pixqa.bitwiseAnd(1).eq(fill))
    if clear is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<1).eq(clear<<1))
    if water is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<2).eq(water<<2))
    if shadow is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<3).eq(shadow<<3))
    if snow is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<4).eq(snow<<4))
    if cloud is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<5).eq(cloud<<5))
    if cloud_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<6).lte(cloud_conf<<6))
    if cirrus_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<8).lte(cirrus_conf<<8))
    if terrain is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<10).eq(terrain<<10))
        
    return mask
    

def sr_mask(img, opt_lo=0, opt_hi=10000, edges=True, aerosol=False, **kwargs):
    """ Update pixel mask based on surface reflectance QA bands.
    opt_low: float or list of floats
        The lowest valid pixel value for all optical bands (float) or each band (list)
    
    opt_hi: float or list of floats
        The highest valid pixel value for all optical (float) or each band (list)
        
    aerosol: bool
        Remove high aerosols in L8. No aerosol band in L4-7 so throws an error.
        
    **kwargs:
        Options passed to pqa_mask function
        
    returns: ee.Image
        Mask
    """
    # TODO: Allow selection of bands to return and to consider when masking
    
    sat = img.get('SATELLITE')
    
    # Apply pqa mask
    pqa = pqa_mask(img, **kwargs)  
    
    # Keep only valid values in the optical range
    # TODO: Mask for thermal too?
    opt = img.select(ee.Dictionary(__srbands.get(sat)).get('opt'))
    valid = opt.lt(opt_lo).Or(opt.gt(opt_hi)).reduce(ee.Reducer.max()).eq(0)
    
    # Remove edge pixels that don't occur in all bands (common in L5)
    if edges:
        edge = opt.mask().reduce(ee.Reducer.min())
    else:
        edge = ee.Image(1)
    
    # Keep pixels without radiometric saturation in any band
    # Note: oversaturation can rollover into valid range, and pixels outside 
    #       valid range don't necessarily register as oversaturated
    # TODO: blue and coastal often oversatured, consider just clamping 0-1
    rad = img.select('radsat_qa').eq(0)
    
    # Mask out high aerosols in Landsat 8. Not possible in L4-7.
    if aerosol:
        # TODO: avoid using ee.Algorithm.If
        has_aero = img.bandNames().contains('sr_aerosol')
        aero = ee.Algorithms.If(has_aero,
                                img.select('sr_aerosol').lt(192),
                                ee.Image(1))
    else:
        aero = ee.Image(1)
    
    # combine masks
    mask = pqa.And(edge).And(valid).And(rad).And(aero)
    
    return mask


def sr_collection(aoi, start, end, startdoy=1, enddoy=366, bands=None, sats=["L4", "L5", "L7", "L8"], harmonize=True,
                  rescale=True, cloud_cover=70, mask_func=sr_mask, slc_on=False, mask_kwargs={}, exclude=None):
    """ Return a Landsat collection with renamed bands and other optional preparations applied.

    Get collection of all landsat images that intersect the given feature,
    feature collection, or geometry for the given time frame. Images are
    prepared by renaming and optionally rescaling bands, harmonizing landsat 8,
    filtering for percent cloud cover, and applying a masking function.

    Parameters
    ----------
    aoi: ee.FeatureCollection, ee.Feature or ee.Geometry
        Area over which to retrieve Landsat images.

    start: ee.Date or date string readable by ee.imageCollection.filterDate
        Start date for filtering, inclusive

    end: ee.Date or date string readable by ee.imageCollection.filterDate
        End date for filtering, exclusive
        
    startdoy: int
        Start day of year passed to ee.Filter.calendarRange. startdoy may be more
        than enddoy to allow for wrapping between years (e.g. Dec 20 to Jan 10)
        
    enddoy: int
        End day of year passed to ee.Filter.calendarRange
        
    bands: list
        list of band names to keep. Returns all bands if None.
        
    sats: list
        list of abbreviations of landsat satellites to include
        
    rescale: bool
        whether to rescale bands to reflectance (i.e. 0-1).
        
    harmonize: bool
        whether to harmonize L8 to L57 following Roy et al. 20??
        
    cloud_cover: int
        Maximum of percent of clouds allowed in an image
        
    mask_func: function
        a pixel masking function to apply to all images.
        
    slc_on: bool
        Keep Landsat 7 images only when SLC is on, January 1998 - May 2003.
        
    exclude: list
        List of system:index to exclude from the collection.
        
    **kwargs:
        arguments passed to mask_func
        
    Returns
    -------
    ee.ImageCollection
        Collection of prepared Landsat images intersecting study area within the 
        given time frame.   
        
    Authors
    -------
    Steven Filippelli
    
    """
    # TODO: Add ability to keep a list of given bands and use them in masking
    
    sat_dict = {
        'L4':{'id':'LANDSAT/LT04/C01/T1_SR',
              'bands':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']},
        'L5':{'id':'LANDSAT/LT05/C01/T1_SR',
              'bands':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']},
        'L7':{'id':'LANDSAT/LE07/C01/T1_SR',
              'bands':['blue', 'green', 'red', 'nir', 'swir1','therm1', 'swir2', 'sr_atmos_opacity', 'sr_cloud_qa', 'pixel_qa', 'radsat_qa']},
        'L8':{'id':'LANDSAT/LC08/C01/T1_SR',
              'bands':['cb', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'therm1', 'therm2', 'sr_aerosol', 'pixel_qa', 'radsat_qa']}
               }
    
    lx_imgs = ee.ImageCollection([])
    for sat in sats:
        sat_bands = sat_dict[sat]['bands']
        imgs = (ee.ImageCollection(sat_dict[sat]['id'])
                .filterBounds(aoi)
                .filterDate(start, end)
                .filter(ee.Filter.calendarRange(startdoy,enddoy))
                .select(list(range(len(sat_bands))), sat_bands)
               )
        if (sat=='L7') and (slc_on):
            imgs = imgs.filterDate(ee.Date.fromYMD(1998,1,1),ee.Date.fromYMD(2003,5,31))
        if cloud_cover<100:
            imgs = imgs.filterMetadata('CLOUD_COVER', 'not_greater_than', cloud_cover)
        if exclude:
            imgs = imgs.filter(ee.Filter.inList('system:index', exclude).Not())
        if (sat=='L8') and (harmonize):
            imgs = imgs.map(l8sr_harmonize)
        if mask_func:
            imgs = imgs.map(lambda i: i.updateMask(mask_func(i, **mask_kwargs)))
        if rescale:
            imgs = imgs.map(sr_rescale)  # TODO: try to move rescale below band selection, need to allow for .select on bands that don't exist in the image in sr_rescale
        if bands:
            imgs = imgs.select(bands) 
        
        lx_imgs = lx_imgs.merge(imgs)

    return lx_imgs    
    
    
    
###############################################################################
#-------------Landsat spectral index functions --------------------------------
# TODO: TSAVI (Baret 1989), TVI (Mcdaniel and Haas 1982), 
# TODO: Move to separate file.
# TODO: allow specification of bands instead of requiring certain band names
# TODO: consider casting all to float

def sri(i):
    """
    Simple Ratio Index (Jordan 1969): https://doi.org/10.2307/1936256
    """
    return (i.expression("nir / red", 
                        {'nir' : i.select('nir'),
                         'red' : i.select('red')
                        })
            .select([0], ['sri']))

def b5b4r(i):
    """
    SWIR1 / NIR ratio. B5 and B4 of TM. (Vogelmann et al., 2009) https://doi.org/10.1016/j.rse.2009.04.014
    """
    return (i.expression("swir1 / nir", 
                        {'nir' : i.select('nir'),
                         'swir1' : i.select('swir1')
                        })
            .select([0], ['b5b4r']))

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


def ndmi(i):
    """
    Normalized Difference Moisture Index (Gao 1996) https://doi.org/10.1016/S0034-4257(96)00067-3
    """
    return (i.normalizedDifference(["nir", "swir1"])
           .select([0], ["ndmi"]))


def nbr(i):
    """
    Normalized Burn Ratio (Key and Benson 2006)  No. RMRS-GTR-164-CD
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
    Carotenoid Reflectance Index 1 (Gitelson 2002). https://doi.org/10.1016/S0034-4257(01)00289-9 
    Hill 2013 RSE found it discriminated treed landscape from others in texas savannahs
    """
    return (i.expression("(1 / blue) - (1 / green)",
                      {'green':i.select('green'),
                       'blue':i.select('blue')
                      })
          .select([0], ["cri1"]))


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


def gcvi(i):
    """
    Green Chlorophyll Vegetation Index (Gitelson et al.,2005) https://doi.org/10.1029/2005GL022688
    Deines et al. 2019 used for mapping irrigation. Works well in most cropland applications.
    """
    return (i.expression("nir / green - 1", 
                        {'nir' : i.select('nir'),
                         'green' : i.select('green')
                        })
            .select([0], ['gcvi']))

def ndsi(i):
    """
    Normalized Difference Snow Index
    """
    return (i.normalizedDifference(["green", "swir1"])
           .select([0], ["ndsi"]))

# ---Tasseled cap band coefficients--------
# TODO: add TOA coefficients for TM TOA

# Baig 2014 coeffs - TOA refl (http://www.tandfonline.com/doi/pdf/10.1080/2150704X.2014.915434)
l8_tc_coeffs = [ee.Image([0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872]),
              ee.Image([-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608]),
              ee.Image([ 0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559]),
              ee.Image([-0.8239, 0.0849, 0.4396, -0.058, 0.2013, -0.2773]),
              ee.Image([-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085]),
              ee.Image([0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252])]

# Huang 2002 coeffs - TOA refl (http://landcover.usgs.gov/pdf/tasseled.pdf)
l7_tc_coeffs = [ee.Image([0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596]),
              ee.Image([-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630]),
              ee.Image([0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388]),
              ee.Image([0.0805, -0.0498, 0.1950, -0.1327, 0.5752, -0.7775]),
              ee.Image([-0.7252, -0.0202, 0.6683, 0.0631, -0.1494, -0.0274]),
              ee.Image([0.4000, -0.8172, 0.3832, 0.0602, -0.1095, 0.0985])]
  
# Crist 1985 coeffs - SR refl (https://doi.org/10.1016/0034-4257(85)90102-6)
# This can be applied to all TM/ETM+ SR data and OLI after harmonization
sr_tc_coeffs = [ee.Image([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]),
                ee.Image([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]),
                ee.Image([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]),
                ee.Image([-0.2117, -0.0284, 0.1302, -0.1007, 0.6529, -0.7078]),
                ee.Image([-0.8669, -0.1835, 0.3856, 0.0408, -0.1132, 0.2272]),
                ee.Image([0.3677, -0.8200, 0.4354, 0.0518, -0.0066, -0.0104])]
            
toa_coeff_dict =ee.Dictionary({"LANDSAT_4": l7_tc_coeffs,
                               "LANDSAT_5": l7_tc_coeffs,
                               "LANDSAT_7": l7_tc_coeffs,
                               "LANDSAT_8": l8_tc_coeffs})
tc_names = ['tcb','tcg', 'tcw', 'tc4','tc5','tc6']


def tasseled_cap(image, sr=True):
    """ Function for computing tasseled cap transform

    Parameters
    ----------
    image: ee.Image
        Image to perform tasseled cap on. Bands must be named by color not num.
        
    sr: bool
        Whether to use surface reflectance coefficients. OLI must be harmonized to TM/ETM+ first.

    Returns
    -------
    ee.Image
        Image with tasseled cap bands.
        
    Authors
    -------
    Ian Housman, adapted by Steven Filippelli
    """
    # TODO: tcap bands lose dimensions. Fix this.
    # TODO: allow computing individual indices rather than full set and discarding rest
    
    image = ee.Image(image).select(["blue", "green", "red", "nir", "swir1", "swir2"])
    # Nested function to do the multiplication and addition
    def mult_sum(matrix):
        return image.multiply(matrix).reduce(ee.call("Reducer.sum")).float()

    # Get the correct coefficients
    if sr:
        coeffs = ee.List(sr_tc_coeffs)
    else:
        coeffs = ee.List(tc_coeff_dict.get(image.get('SATELLITE')))
                  
    # Compute the tc transformation and give it useful names
    tco = coeffs.map(mult_sum)
    empty = ee.Image().select()
    tco = ee.Image(tco.iterate(lambda image, result: 
                              ee.Image(result).addBands(image),
                              empty))
  
    tco = tco.select([0,1,2,3,4,5], tc_names)
    
    return tco


def greenness(i, sr=True):
    return tasseled_cap(i, sr).select("tcg")


def wetness(i, sr=True):
    return tasseled_cap(i, sr).select("tcw")


def brightness(i, sr=True):
    return tasseled_cap(i, sr).select("tcb")

def angle(i, sr=True):
    tcb = tasseled_cap(i, sr).select("tcb")
    tcg = tasseled_cap(i, sr).select("tcg")
    return (tcg.divide(tcb)).atan().multiply(180/np.pi).multiply(100).rename(["tca"])

def specixs(image, ixlist='all', sr=True):
    """ Get landsat spectral indices

    Parameters
    ----------
    image: ee.Image
        Landsat 5-8 image with bands labeled by color rather than band number.

    indices: list
        List of desired band indices to return in a new image.

     Returns
     -------
     ee.Image
         Image with given band indices.

    Author
    ------
    Steven Filippelli
    """
    
    ixbands = {"sri":sri(image),
               "b5b4r":b5b4r(image),
               "ndvi":ndvi(image),
               "evi":evi(image),
               "savi":savi(image),
               "msavi":msavi(image),
               "ndmi":ndmi(image),
               "nbr":nbr(image),
               "nbr2":nbr2(image),
               "cri1":cri1(image),
               "satvi":satvi(image),
               "gcvi":gcvi(image),
               "ndsi":ndsi(image),
               "tasseled_cap":tasseled_cap(image, sr),
               "tcg":greenness(image, sr),
               "tcw":wetness(image, sr),
               "tcb":brightness(image, sr),
               "tca":angle(image, sr)
            }
    
    if ixlist == 'all':
        to_drop = ['tcg', 'tcw', 'tcb', 'tca']
        ixlist = [b for b in ixbands.keys() if b not in to_drop]
    if ixlist == 'no_tcap':
        to_drop= ['tasseled_cap', 'tcg', 'tcw', 'tcb', 'tca']
        ixlist = [b for b in ixbands.keys() if b not in to_drop]
    
    ixbands = [img for b, img in ixbands.items() if b in ixlist]
    iximg = ee.Image.cat(ixbands)
                              
    return iximg
