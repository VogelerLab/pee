# -*- coding: utf-8 -*-
"""

Functions and scripts related to processing synthetic aperture radar (SAR) data.

Created on Wed Feb 7 2024

@author: Steven Filippelli
"""

def dn_to_db(image, dataset):
    """
    Calibrate from DNs to decibels

    Parameters
    ----------
    image : ee.Image
        Image to calibrated

    dataset: str
        The dataset the image belongs to for determining the calibration equation

    Returns
    -------
    ee.Image
        Calibrated image in decibels

    TODO: could select calibration automatically from id, but need to do without ee.Algorithm.If which is slow
    
    """
    bands = image.bandNames().filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    img = image.select(bands)
    
    if dataset=='PALSAR_Yearly':
        out = img.pow(2).log10().multiply(10).subtract(83)

    return image.addBands(out, None, True)


def dn_to_pow(image, dataset):
    """
    Calibrate from DNs directly to linear power

    Parameters
    ----------
    image : ee.Image
        Image to calibrated

    dataset: str
        The dataset the image belongs to for determining the calibration equation

    Returns
    -------
    ee.Image
        Calibrated image in linear power, which is used in speckle filtering

    TODO: could select calibration automatically from id, but need to do without ee.Algorithm.If which is slow
    
    """

    bands = image.bandNames().filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    img = image.select(bands)

    if dataset=='PALSAR_Yearly':
        out = ee.Image(img).pow(2).divide(ee.Number(10.0).pow(8.3))
    
    return image.addBands(out, None, True)


def pow_to_db(image):
    """ Convert from linear power to decibels
    """
    bands = image.bandNames().filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    out = ee.Image(image.select(bands)).log10().multiply(10.0)
    return image.addBands(out, None, True)


def db_to_pow(image):
    """ Convert from decibels to linear power
    """
    bands = image.bandNames().filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    out = ee.Image(10.0).pow(image.select(bands).divide(10.0))
    return image.addBands(out, None, True)


def epoch_to_doy(image, dateband):
    """ Convert date band of image from millis to day of year
    """
    ibands = image.bandNames().filter(ee.Filter.neq('item', dateband))
    ymil = ee.Image.constant(ee.Date.fromYMD(image.date().get('year'), 1, 1).millis())
    imil = image.select(dateband)
    doy = imil.subtract(ymil).divide(1000*60*60*24).rename('doy')
    return image.select(ibands).addBands(doy)
    

"""
Speckle filtering functions below are a copy and modification of
https://github.com/adugnag/gee_s1_ard/blob/main/python-api/speckle_filter.py
to make them sensor agnostic.

Version: v1.0
Date: 2021-03-12
Authors: Mullissa A., Vollrath A., Braun, C., Slagter B., Balling J., Gou Y., Gorelick N.,  Reiche J.
Description: A collection of functions to perform mono-temporal and multi-temporal speckle filtering
"""
import ee
import math

# ---------------------------------------------------------------------------//
# 1.SPECKLE FILTERS
# ---------------------------------------------------------------------------//

def boxcar(image, ksize=3):
    """
    Apply boxcar filter on every image in the collection.

    Parameters
    ----------
    image : ee.Image
        Image to be filtered
    ksize : positive odd integer
        Neighbourhood window size

    Returns
    -------
    ee.Image
        Filtered Image

    """
    
    bands = image.bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    
      #Define a boxcar kernel
    kernel = ee.Kernel.square((ksize/2), units='pixels', normalize=True)
     #Apply boxcar
    output = image.select(ibands).convolve(kernel).rename(ibands)
    return image.addBands(output, None, True)
    

def leefilter(image, ksize=3, enl=5):
    """
    Lee Filter applied to one image.
    It is implemented as described in
    J. S. Lee, “Digital image enhancement and noise filtering by use of local statistics,”
    IEEE Pattern Anal. Machine Intell., vol. PAMI-2, pp. 165–168, Mar. 1980.

    Parameters
    ----------
    image : ee.Image
        Image to be filtered
    ksize : positive odd integer
        Neighbourhood window size
    enl: integer
        Equivalent number of looks of the image pixel. Default 5 is for S1-GRD. PALSAR-1 mosaics should be 16.


    Returns
    -------
    ee.Image
        Filtered Image

    """
    
    bands = image.bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))

    # Compute the speckle standard deviation
    eta = 1.0/math.sqrt(enl)
    eta = ee.Image.constant(eta)

    # MMSE estimator
    # Neighbourhood mean and variance
    oneImg = ee.Image.constant(1)
    # Estimate stats
    reducers = ee.Reducer.mean().combine(
                      reducer2= ee.Reducer.variance()
                      ,sharedInputs= True
                      )
    stats = (image.select(ibands).reduceNeighborhood(
                      reducer= reducers
                          ,kernel= ee.Kernel.square(ksize/2, 'pixels')
                              ,optimization= 'window'))
    meanBand = ibands.map(lambda bname: ee.String(bname).cat('_mean'))
    varBand = ibands.map(lambda bname:  ee.String(bname).cat('_variance'))

    z_bar = stats.select(meanBand)
    varz = stats.select(varBand)
    # Estimate weight 
    varx = (varz.subtract(z_bar.pow(2).multiply(eta.pow(2)))).divide(oneImg.add(eta.pow(2)))
    b = varx.divide(varz)
  
    # if b is negative set it to zero
    new_b = b.where(b.lt(0), 0)
    output = oneImg.subtract(new_b).multiply(z_bar.abs()).add(new_b.multiply(image.select(ibands))).rename(ibands)
    return image.addBands(output, None, True)


def gammamap(image, ksize, enl=5): 
    
    """
    Gamma Maximum a-posterior Filter applied to one image. It is implemented as described in 
    Lopes A., Nezry, E., Touzi, R., and Laur, H., 1990.  
    Maximum A Posteriori Speckle Filtering and First Order texture Models in SAR Images.  
    International  Geoscience  and  Remote  Sensing  Symposium (IGARSS).
    Parameters
    ----------
    image : ee.Image
        Image to be filtered
    ksize : positive odd integer
        Neighbourhood window size
    Returns
    -------
    ee.Image
        Filtered Image
    """
    
    bands = image.bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    
    #local mean
    reducers = ee.Reducer.mean().combine( \
                      reducer2= ee.Reducer.stdDev(), \
                      sharedInputs= True
                      )
    stats = (image.select(ibands).reduceNeighborhood( \
                      reducer= reducers, \
                          kernel= ee.Kernel.square(ksize/2,'pixels'), \
                              optimization= 'window'))
    meanBand = ibands.map(lambda bandName: ee.String(bandName).cat('_mean'))
    stdDevBand = ibands.map(lambda bandName:  ee.String(bandName).cat('_stdDev'))
        
    z = stats.select(meanBand)
    sigz = stats.select(stdDevBand)
    
    #local observed coefficient of variation
    ci = sigz.divide(z)
    #noise coefficient of variation (or noise sigma)
    cu = 1.0/math.sqrt(enl)
    #threshold for the observed coefficient of variation
    cmax = math.sqrt(2.0) * cu
    cu = ee.Image.constant(cu)
    cmax = ee.Image.constant(cmax)
    enlImg = ee.Image.constant(enl)
    oneImg = ee.Image.constant(1)
    twoImg = ee.Image.constant(2)

    alpha = oneImg.add(cu.pow(2)).divide(ci.pow(2).subtract(cu.pow(2)))

    #Implements the Gamma MAP filter described in equation 11 in Lopez et al. 1990
    q = image.select(ibands).expression('z**2 * (z * alpha - enl - 1)**2 + 4 * alpha * enl * b() * z', { 'z': z,  'alpha':alpha,'enl': enl})
    rHat = z.multiply(alpha.subtract(enlImg).subtract(oneImg)).add(q.sqrt()).divide(twoImg.multiply(alpha))
  
    #if ci <= cu then its a homogenous region ->> boxcar filter
    zHat = (z.updateMask(ci.lte(cu))).rename(ibands)
    #if cmax > ci > cu then its a textured medium ->> apply Gamma MAP filter
    rHat = (rHat.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))).rename(ibands)
    #ci>cmax then its strong signal ->> retain
    x = image.select(ibands).updateMask(ci.gte(cmax)).rename(ibands)  
    #Merge
    output = ee.ImageCollection([zHat,rHat,x]).sum()
    return image.addBands(output, None, True)

def RefinedLee(image):
    """
    This filter is modified from the implementation by Guido Lemoine 
    Source: Lemoine et al. https://code.earthengine.google.com/5d1ed0a0f0417f098fdfd2fa137c3d0c

    Parameters
    ----------
    image: ee.Image
        Image to be filtered. Must be in linear units, not dB!.

    Returns
    -------
    result: ee.Image
        Filtered Image

    """

    bands = image.bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))

    def inner(b):

        img = image.select([b]);
    
        # img must be linear, i.e. not in dB!
        # Set up 3x3 kernels 
        weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
        kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False);
  
        mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
        variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
  
        # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
        sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
  
        sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False);
  
        # Calculate mean and variance for the sampled windows and store as 9 bands
        sample_mean = mean3.neighborhoodToBands(sample_kernel); 
        sample_var = variance3.neighborhoodToBands(sample_kernel);
  
        # Determine the 4 gradients for the sampled windows
        gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
        gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
        gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
        gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
  
        # And find the maximum gradient amongst gradient bands
        max_gradient = gradients.reduce(ee.Reducer.max());
  
        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(max_gradient);
  
        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask);
  
        # Determine the 8 directions
        directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
        directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
        directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
        directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).Not().multiply(5));
        directions = directions.addBands(directions.select(1).Not().multiply(6));
        directions = directions.addBands(directions.select(2).Not().multiply(7));
        directions = directions.addBands(directions.select(3).Not().multiply(8));
  
        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask);
  
        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum());  
  
        sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
  
        #Calculate localNoiseVariance
        sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
  
        # Set up the 7*7 kernels for directional statistics
        rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
  
        diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);
  
        rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False);
        diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False);
  
        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
        dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
  
        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
        dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
  
        # and add the bands for rotated kernels
        for i in range(1, 4):
            dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
            dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))

  
        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
        dir_mean = dir_mean.reduce(ee.Reducer.sum());
        dir_var = dir_var.reduce(ee.Reducer.sum());
  
        # A finally generate the filtered value
        varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))
  
        b = varX.divide(dir_var)
        result = dir_mean.add(b.multiply(img.subtract(dir_mean)))
  
        return result.arrayProject([0]).arrayFlatten([['sum']]).float()
    
    result = ee.ImageCollection(ibands.map(inner)).toBands().rename(ibands).copyProperties(image)
    
    return image.addBands(result, None, True) 



def leesigma(image, ksize, enl=5):
    """
    Implements the improved lee sigma filter to one image. 
    It is implemented as described in, Lee, J.-S. Wen, J.-H. Ainsworth, T.L. Chen, K.-S. Chen, A.J. 
    Improved sigma filter for speckle filtering of SAR imagery. 
    IEEE Trans. Geosci. Remote Sens. 2009, 47, 202–213.

    Parameters
    ----------
    image : ee.Image
        Image to be filtered
    ksize : positive odd integer
        Neighbourhood window size

    Returns
    -------
    ee.Image
        Filtered Image

    """

    #parameters
    Tk = ee.Image.constant(7) #number of bright pixels in a 3x3 window
    sigma = 0.9
    
    target_kernel = 3
    bands = image.bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
  
    #compute the 98 percentile intensity 
    z98 = ee.Dictionary(image.select(ibands).reduceRegion(
                reducer= ee.Reducer.percentile([98]),
                geometry= image.geometry(),
                scale=10,
                maxPixels=1e13
            )).toImage()
  

    #select the strong scatterers to retain
    brightPixel = image.select(ibands).gte(z98)
    K = brightPixel.reduceNeighborhood(ee.Reducer.countDistinctNonNull()
            ,ee.Kernel.square(target_kernel/2)) 
    retainPixel = K.gte(Tk)

  
    #compute the a-priori mean within a 3x3 local window
    #original noise standard deviation since the data is 5 look
    eta = 1.0/math.sqrt(enl) 
    eta = ee.Image.constant(eta)
    #MMSE applied to estimate the apriori mean
    reducers = ee.Reducer.mean().combine( \
                      reducer2= ee.Reducer.variance(), \
                      sharedInputs= True
                      )
    stats = image.select(ibands).reduceNeighborhood( \
                      reducer= reducers, \
                          kernel= ee.Kernel.square(target_kernel/2,'pixels'), \
                              optimization= 'window')
    meanBand = ibands.map(lambda bandName: ee.String(bandName).cat('_mean'))
    varBand = ibands.map(lambda bandName:  ee.String(bandName).cat('_variance'))
        
    z_bar = stats.select(meanBand)
    varz = stats.select(varBand)
    
    oneImg = ee.Image.constant(1)
    varx = (varz.subtract(z_bar.abs().pow(2).multiply(eta.pow(2)))).divide(oneImg.add(eta.pow(2)))
    b = varx.divide(varz)
    xTilde = oneImg.subtract(b).multiply(z_bar.abs()).add(b.multiply(image.select(ibands)))
  
    #step 3: compute the sigma range
    #Lookup table (J.S.Lee et al 2009) for range and eta values for intensity (only 4 look is shown here)
    LUT = ee.Dictionary({0.5: ee.Dictionary({'I1': 0.694,'I2': 1.385,'eta': 0.1921}),
                                 0.6: ee.Dictionary({'I1': 0.630,'I2': 1.495,'eta': 0.2348}),
                                 0.7: ee.Dictionary({'I1': 0.560,'I2': 1.627,'eta': 0.2825}),
                                 0.8: ee.Dictionary({'I1': 0.480,'I2': 1.804,'eta': 0.3354}),
                                 0.9: ee.Dictionary({'I1': 0.378,'I2': 2.094,'eta': 0.3991}),
                                 0.95: ee.Dictionary({'I1': 0.302,'I2': 2.360,'eta': 0.4391})});
  
    #extract data from lookup
    sigmaImage = ee.Dictionary(LUT.get(str(sigma))).toImage()
    I1 = sigmaImage.select('I1')
    I2 = sigmaImage.select('I2')
    #new speckle sigma
    nEta = sigmaImage.select('eta')
    #establish the sigma ranges
    I1 = I1.multiply(xTilde)
    I2 = I2.multiply(xTilde)
  
    #step 3: apply MMSE filter for pixels in the sigma range
    #MMSE estimator
    mask = image.select(ibands).gte(I1).Or(image.select(ibands).lte(I2))
    z = image.select(ibands).updateMask(mask)
  
    stats = z.reduceNeighborhood( \
                      reducer= reducers, \
                          kernel= ee.Kernel.square(ksize/2,'pixels'), \
                              optimization= 'window')
        
    z_bar = stats.select(meanBand)
    varz = stats.select(varBand)
    
    
    varx = (varz.subtract(z_bar.abs().pow(2).multiply(nEta.pow(2)))).divide(oneImg.add(nEta.pow(2)))
    b = varx.divide(varz)
    #if b is negative set it to zero
    new_b = b.where(b.lt(0), 0)
    xHat = oneImg.subtract(new_b).multiply(z_bar.abs()).add(new_b.multiply(z))
  
    #remove the applied masks and merge the retained pixels and the filtered pixels
    xHat = image.select(ibands).updateMask(retainPixel).unmask(xHat)
    output = ee.Image(xHat).rename(ibands)
    return image.addBands(output, None, True)


#---------------------------------------------------------------------------//
# 2. MONO-TEMPORAL SPECKLE FILTER (WRAPPER)
#---------------------------------------------------------------------------//


def MonoTemporal_Filter(coll, sfilter='BOXCAR', ksize=3, enl=5):
    """
    A wrapper function for monotemporal filter

    Parameters
    ----------
    coll : ee Image collection
        the image collection to be filtered
    ksize : odd integer
        Spatial Neighbourhood window
    sfilter : String
        Type of speckle filter

    Returns
    -------
    ee.ImageCollection
        An image collection where a mono-temporal filter is applied to each 
        image individually

    """
    def _filter(image):    
        if (sfilter=='BOXCAR'):
            _filtered = boxcar(image, ksize)
        elif (sfilter=='LEE'):
            _filtered = leefilter(image, ksize, enl)
        elif (sfilter=='GAMMA MAP'):
            _filtered = gammamap(image, ksize, enl)
        elif (sfilter=='REFINED LEE'):
            _filtered = RefinedLee(image)
        elif (sfilter=='LEE SIGMA'):
            _filtered = leesigma(image, ksize, enl)
        return _filtered
    return coll.map(_filter)

# ---------------------------------------------------------------------------//
# 3. MULTI-TEMPORAL SPECKLE FILTER
# ---------------------------------------------------------------------------//

def MultiTemporal_Filter(coll, sfilter='BOXCAR', ksize=3, enl=5):
    """

    A wrapper function for multi-temporal filter

    Parameters
    ----------
    coll : ee Image collection
        the image collection to be filtered
    sfilter : String
        Type of speckle filter
    ksize : odd integer
        Spatial Neighbourhood window

    Returns
    -------
    ee.ImageCollection
        An image collection where a multi-temporal filter is applied to each
        image individually

    """
    bands = coll.first().bandNames()
    ibands = bands.filter(ee.Filter.inList('item', ['HH', 'HV', 'VV', 'VH']))
    
    icoll = coll.select(ibands)
    meanBands = ibands.map(lambda bandName: ee.String(bandName).cat('_mean'))
    ratioBands = ibands.map(lambda bandName: ee.String(bandName).cat('_ratio'))
    count_img = icoll.reduce(ee.Reducer.count())

    def inner(image):
        """
        Creats an image whose bands are the filtered image and image ratio

        Parameters
        ----------
        image : ee.Image
            Image to be filtered

        Returns
        -------
        ee.Image
            Filtered image and image ratio

        """
        
        if (sfilter=='BOXCAR'):
            _filtered = boxcar(image, ksize).select(ibands).rename(meanBands) 
        elif (sfilter=='LEE'):
            _filtered = leefilter(image, ksize, enl).select(ibands).rename(meanBands)
        elif (sfilter=='GAMMA MAP'):
            _filtered = gammamap(image, ksize, enl).select(ibands).rename(meanBands)
        elif (sfilter=='REFINED LEE'):
            _filtered = RefinedLee(image).select(ibands).rename(meanBands)
        elif (sfilter=='LEE SIGMA'):
            _filtered = leesigma(image, ksize, enl).select(ibands).rename(meanBands)

        _ratio = image.select(ibands).divide(_filtered).rename(ratioBands) 
        return _filtered.addBands(_ratio)
        
    ifilt = icoll.map(inner)
    isum = ifilt.select(ratioBands).reduce(ee.Reducer.sum())
    ocoll = ifilt.select(meanBands).map(lambda i: i.divide(count_img).multiply(isum).rename(ibands))
    out =  coll.combine(ocoll, overwrite=True)
    out = out.map(lambda i: i.set('system:time_start', i.date().millis())) # for some reason timestart property getting made unusable (maybe temporary bug in GEE datasets)
    return out