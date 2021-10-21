# -*- coding: utf-8 -*- 
# Functions for preparing Landsat data in GEE. 

# TODO: Add function to get a scene's snow and cloud cover from pqa band (see redcedar_landtrendr_v4)

import numpy as np
import ee
import time_series
ee.Initialize()


__sr_dict = {
    'LANDSAT_4':{'id':'LANDSAT/LT04/C02/T1_L2',
                 'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'st_other':['ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD'],
                 'atmos':['SR_ATMOS_OPACITY'],
                 'qa':['SR_CLOUD_QA', 'QA_PIXEL', 'QA_RADSAT'],
                 'orig':['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'ST_B6', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'st', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT']
                },
    'LANDSAT_5':{'id':'LANDSAT/LT05/C02/T1_L2',
                 'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'st_other':['ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD'],
                 'atmos':['SR_ATMOS_OPACITY'],
                 'qa':['SR_CLOUD_QA', 'QA_PIXEL', 'QA_RADSAT'],
                 'orig':['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'ST_B6', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'st', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT']
                },
    'LANDSAT_7':{'id':'LANDSAT/LE07/C02/T1_L2',
                 'opt':['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'st_other':['ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD'],
                 'atmos':['SR_ATMOS_OPACITY'],
                 'qa':['SR_CLOUD_QA', 'QA_PIXEL', 'QA_RADSAT'],
                 'orig':['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'ST_B6', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT'],
                 'all':['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'SR_ATMOS_OPACITY', 'SR_CLOUD_QA', 'st', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT']
                },
    'LANDSAT_8':{'id':'LANDSAT/LC08/C02/T1_L2',
                 'opt':['cb','blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
                 'st_other':['ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD'],
                 'atmos':[],
                 'qa':['SR_QA_AEROSOL', 'QA_PIXEL', 'QA_RADSAT'],
                 'orig':['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'SR_QA_AEROSOL', 'ST_B10', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT'],
                 'all':['cb', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'SR_QA_AEROSOL', 'st', 'ST_ATRAN', 'ST_CDIST', 'ST_DRAD', 'ST_EMIS', 'ST_EMSD', 'ST_QA', 'ST_TRAD', 'ST_URAD', 'QA_PIXEL', 'QA_RADSAT']
                }   
    }
__sr_eedict = ee.Dictionary(__sr_dict)


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
    'LANDSAT_8': {'SR_B1':'cb', 'SR_B2':'blue', 'SR_B3':'green', 'SR_B4':'red', 'SR_B5':'nir', 
                  'SR_B6':'swir1', 'SR_B7':'swir2', 'ST_B10':'st'},
    'LANDSAT_7': {'SR_B1':'blue', 'SR_B2':'green', 'SR_B3':'red', 'SR_B4':'nir', 
                  'SR_B5':'swir1', 'SR_B7':'swir2', 'ST_B6':'st'},
    'LANDSAT_5': {'SR_B1':'blue', 'SR_B2':'green', 'SR_B3':'red', 'SR_B4':'nir', 
                  'SR_B5':'swir1', 'SR_B7':'swir2', 'ST_B6':'st'},
    'LANDSAT_4': {'SR_B1':'blue', 'SR_B2':'green', 'SR_B3':'red', 'SR_B4':'nir', 
                  'SR_B5':'swir1', 'SR_B7':'swir2', 'ST_B6':'st'}
    })
    band_dict = ee.Dictionary(band_dicts.get(image.get('SPACECRAFT_ID')))

    def swap_name(b):
        # TODO: try to use filter instead of if
        return ee.Algorithms.If(band_dict.contains(b), band_dict.get(b), b)
  
    new_bands = image.bandNames().map(swap_name)
    return image.rename(new_bands)


def sr_rename(img):
    """ Rename bands based on SPACECRAFT_ID, assuming all bands present in original order.
    """
    bnames = ee.Dictionary(__sr_eedict.get(img.get('SPACECRAFT_ID'))).get('all')
    return img.rename(bnames)


def sr_rescale(img):
    """ Rescale Landsat SR images from any satellite. Ancillary surface temperature bands are not scaled.
    
    TODO: Add option to scale other ST bands - unlikely to use these, so unneccessary compute time for most cases; other fix would be to apply scale only to present bands and perform select beforehand.
    """
    sat = img.get('SPACECRAFT_ID')
    bdict = ee.Dictionary(__sr_eedict.get(sat))
    opt = img.select(bdict.get('opt'))
    therm = img.select(['st'])
    atmos = img.select(bdict.get('atmos'))
    qa = img.select(bdict.get('qa'))
    st_other = img.select(bdict.get('st_other'))
    
    # scale and cast to float (may reduce memory bottlenecks?)
    opt = opt.multiply(0.0000275).add(-0.2).toFloat()
    atmos = atmos.multiply(__atmos_scale.get(sat)).toFloat() #0.001 for L57, empty for L8
    therm = therm.multiply(0.00341802).add(149.0).toFloat()

    scaled = ee.Image(ee.Image.cat([opt, therm, qa, atmos, st_other])).copyProperties(img)
    scaled = scaled.set('system:time_start', img.get('system:time_start'))
    return scaled

# TODO: Need to check coefficients and update to use reflectance 0-1, but this function likely no longer necessary anyway.
# def l8sr_harmonize(img):
#     """Harmonize L8 to L5/7 following Roy 2016 for surface reflectance images.
#        Expects images in the original scale (i.e.0-10000) with all renamed bands.
#        This correction likely makes differences worse for Collection 2: https://groups.google.com/g/google-earth-engine-developers/c/mx5LmM9SI2g/m/nsg7fqv5AQAJ
#     """
#     # RMA slopes and intercepts per band
#     slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949])
# #     itcp_coefs = [-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029] # img scaled 0-1
#     itcp_coefs = [-95.0, -16.0, -22.0, -21.0, -30.0,  29.0] # img scaled 0-10,000
#     itcp = ee.Image.constant(itcp_coefs)
# #     # Least squares slopes and intercepts per band
# #     slopes = ee.Image.constant([0.885, 0.9317, 0.9372, 0.8339, 0.8639, 0.9165])
# #     itcp = ee.Image.constant([0.0183, 0.0123, 0.0123, 0.0448, 0.0306, 0.0116])

#     # apply reverse of regression (Roy gives ETM to OLI)
#     y = (img.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
#             .subtract(itcp).divide(slopes)
#             .toInt16()
#         )
#     new = ee.Image.cat([img.select(['cb']),
#                         y,
#                         img.select('st', 'sr_aerosol', 'QA_PIXEL', 'QA_RADSAT')])
    
#     img = (new.copyProperties(img)
#              .set('system:time_start', img.get('system:time_start'))
#           )
# #     # SF: now treat as L5/7 for calculating TCAP
# #     img = img.set('SPACECRAFT_ID', 'LANDSAT_7')
#     return img


def pqa_mask(img, fill=0, dilated_cloud=0, cirrus=0, cloud=0, shadow=0, snow=0, clear=None, water=0, cloud_conf=1, shadow_conf=1, snow_conf=1, cirrus_conf=1):
    """Keep pixels which meet all the given critera, and mask out others.
    Values should match the Landsat 8 and Landsat 4-7 product guide QA table with 
    0=='No' and 1=='Yes', None='accept either' and conf values: 0,1,2,3=='None', 'low', 'med', 'high'
    https://www.usgs.gov/media/files/landsat-8-collection-2-level-2-science-product-guide
    https://www.usgs.gov/media/files/landsat-4-7-collection-2-level-2-science-product-guide
   
    'Less than or equal to' operation applied to conf values. For example, 
    cloud_conf=2 means keep medium and low confidence clouds. 
    
    Defaults retain non-fill pixels with no water, snow, clouds, and only low or no confidence for all the confidence bits. This is conservative as medium confidence clouds are masked out.
    
    L4-7 does not have cirrus or cirrus confidence bits. Setting cirrus=1 results in a completely masked out image for L4-7.
    
    returns: ee.Image
        Mask which can be applied using ee.Image.updateMask()
    """    
    pixqa = img.select('QA_PIXEL')
    
    mask = img.mask() # use existing mask or ee.Image(1)?
    if fill is not None:
        mask = mask.And(pixqa.bitwiseAnd(1).eq(fill))
    if dilated_cloud is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<1).eq(dilated_cloud<<1))
    if cirrus is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<2).eq(cirrus<<2))
    if cloud is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<3).eq(cloud<<3))                    
    if shadow is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<4).eq(shadow<<4))
    if snow is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<5).eq(snow<<5))
    if clear is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<6).eq(clear<<6))
    if water is not None:
        mask = mask.And(pixqa.bitwiseAnd(1<<7).eq(water<<7))
    if cloud_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<8).lte(cloud_conf<<8))
    if shadow_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<10).lte(shadow_conf<<10))
    if snow_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<12).lte(snow_conf<<12))
    if cirrus_conf is not None:
        mask = mask.And(pixqa.bitwiseAnd(3<<14).lte(cirrus_conf<<14))
        
    return mask
    

def sr_mask(img, opt_lo=7273, opt_hi=43636, edges=True, aerosol=True, **kwargs):
    """ Return pixel mask based on surface reflectance QA bands and valid ranges.
    Assumes optical bands have been renamed as "cb", "blue", etc.
    
    opt_low: float or list of floats
        The lowest valid pixel value for all optical bands (float) or each band (list)
    
    opt_hi: float or list of floats
        The highest valid pixel value for all optical (float) or each band (list)
        
    aerosol: bool
        Remove high aerosols in L8. Ignored for L4-7 which do not have sr_qa_aerosol.
        
    **kwargs:
        Options passed to pqa_mask function
        
    returns: ee.Image
        Mask
    """
    # TODO: Allow selection of bands to return and to consider when masking
    
    sat = img.get('SPACECRAFT_ID')
    
    # Apply pqa mask
    pqa = pqa_mask(img, **kwargs)  
    
    # Keep only valid values in the optical range
    # TODO: Mask for thermal too?
    opt = img.select(ee.Dictionary(__sr_eedict.get(sat)).get('opt'))
    valid = opt.lt(opt_lo).Or(opt.gt(opt_hi)).reduce(ee.Reducer.max()).eq(0)
    
    # Remove edge pixels that don't occur in all bands (common in L5)
    if edges:
        edge = opt.mask().reduce(ee.Reducer.min())
    else:
        edge = ee.Image(1)
    
    # Keep pixels without radiometric saturation in any band. 
    # Note: oversaturation can rollover into valid range, and pixels outside 
    #       valid range don't necessarily register as oversaturated
    # Note 2: Terrain occlusion included in radsat band for Landsat 8,
    #         so enforcing 0 also filters for no terrain occlusion.
    # TODO: blue and coastal often oversatured, consider just clamping 0-1
    # TODO: QA_RADSAT bitmask could be used to only filter for certain bands.
    rad = img.select('QA_RADSAT').eq(0)
    
    # Mask out high aerosols in Landsat 8. Not possible in L4-7.
    if aerosol:
        # TODO: avoid using ee.Algorithm.If
        has_aero = img.bandNames().contains('SR_QA_AEROSOL')
        aero = ee.Algorithms.If(has_aero,
                                img.select('SR_QA_AEROSOL').lt(192),
                                ee.Image(1))
    else:
        aero = ee.Image(1)
    
    # combine masks
    mask = pqa.And(edge).And(valid).And(rad).And(aero)
    
    return mask


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


def sr_collection(aoi, start, end, startdoy=1, enddoy=366, bands=None, 
                  sats=["LANDSAT_4", "LANDSAT_5", "LANDSAT_7", "LANDSAT_8"],
                  rescale=True, cloud_cover=70, mask_func=sr_mask, slc_on=False,
                  tdom=False, exclude=None, mask_kwargs={}, tdom_kwargs={}
                 ):
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
        
    cloud_cover: int
        Maximum of percent of clouds allowed in an image
        
    mask_func: function
        a pixel masking function to apply to all images.
        
    slc_on: bool
        Keep Landsat 7 images only when SLC is on, January 1998 - May 2003.
        
    tdom: bool
        Apply Temporal Dark Outlier Mask to remove cloud shadows
        
    exclude: list
        List of system:index to exclude from the collection.
        
    mask_kwargs: dict
        keyword arguments passed to mask_func
        
    tdom_kwargs: dict
        keyword arguments passed to tdom
        
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
    
    lx_imgs = ee.ImageCollection([])
    for sat in sats:
        sat_bands = __sr_dict[sat]['all']
        imgs = (ee.ImageCollection(__sr_dict[sat]['id'])
                .filterBounds(aoi)
                .filterDate(start, end)
                .filter(ee.Filter.calendarRange(startdoy,enddoy))
                .select(list(range(19)), sat_bands)
               )
        if (sat=='LANDSAT_7') and (slc_on):
            imgs = imgs.filterMetadata('SENSOR_MODE_SLC', 'equals', 'ON')
        if cloud_cover<100:
            imgs = imgs.filterMetadata('CLOUD_COVER', 'not_greater_than', cloud_cover)
        if exclude:
            imgs = imgs.filter(ee.Filter.inList('system:index', exclude).Not())
        if mask_func:
            imgs = imgs.map(lambda i: i.updateMask(mask_func(i, **mask_kwargs)))
        if rescale:
            imgs = imgs.map(sr_rescale)  # TODO: try to move rescale below band selection, need to allow for .select on bands that don't exist in the image in sr_rescale
        if tdom:
            # tdom should be after rescale otherwise defaults for sum_thresh will need to be changed for typical runs
            imgs = tdom2(imgs, **tdom_kwargs)
        if bands:
            # TODO: allow selection of bands that don't exist while keeping those in 'bands' that do? Or keep error?
            imgs = imgs.select(bands)
        
        lx_imgs = lx_imgs.merge(imgs)

    return lx_imgs
    
    
###############################################################################
#-------------Landsat spectral index functions --------------------------------
# TODO: TSAVI (Baret 1989), TVI (Mcdaniel and Haas 1982), 
# TODO: Move to separate file since most are sensor agnostic, except TCAP.
# TODO: allow specification of bands instead of requiring certain band names
# TODO: consider casting all to float
# TODO: create params for rescaled 0-1 true/false since some depend on rescaled

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
    return (i.expression("2.5 * ((nir - red) / (nir + 6.0 * red - 7.5 * blue + 1.))",
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
# TODO: add crist and cicone 1984 coefficients

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
        Whether to use surface reflectance coefficients. For Collection 2 OLI harmonization with TM/ETM+ may not be necessary.

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
        coeffs = ee.List(toa_coeff_dict.get(image.get('SPACECRAFT_ID')))
                  
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
        
    sr: bool
        Image is surface reflectance

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
