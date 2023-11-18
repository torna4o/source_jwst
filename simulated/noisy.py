import numpy as np

# Generate Gaussian noise

def add_gauss(img, mean=0, std=1, speckle = False):

    """
    This function adds a gaussian noise to a numpy.ndarray.
    It takes a numpy.ndarray as input and returns gaussian added array
    and gaussian noise itself separately.
    
    "img" = a numpy.ndarray format image data in 2D
    "mean" = mean of the gaussian noise array to be generated.
    "std" = standard deviation of the gaussian noise array to be generated.
    """

    gauss = np.random.normal(loc=mean,scale=std, size=(img.shape[0],img.shape[1]))

    if speckle==False:
        img_gauss = img + gauss
    elif speckle==True:
        img_gauss = img + img * gauss
    return img_gauss, gauss

def add_poisson(img, lam):
    
    """
    
    This function adds a Poisson noise to a numpy.ndarray.
    It takes a numpy.ndarray as input and returns Poisson added array
    and Poisson noise itself separately.

    "img" = a numpy.ndarray format image data in 2D   
    "lam" = lambda of Poisson distribution.

    """

    pois = np.random.poisson(lam, (img.shape[0], img.shape[1]))
    img_pois = img + pois
    
    return img_pois, pois

def saltpep(img, multi, threshold=0.99):
    """
    Adds salt & pepper noise using NumPy package
    
    img = image data in numpy.ndarray that the noise will be ingested
    threshold = cutoff for making salt and pepper noise, 1 will result
        in full salt, 0 full pepper noise
    multi = this is a multiplier, in case one needs different numbers than 
        1 salt noise generation
    """
    base = np.random.rand(img.shape[0], img.shape[1])
    base[base >= threshold] = 1
    base[base < threshold] = 0
    base *= multi
    return np.add(img, base), base


def rstrip(img, dire, s=1):
    """
    Creates single row (dire='h') or column (dire='v') of random numbers scaled
    between 0 and s (default is s=1). Then, it broadcast the strips in the same shape as img
    
    img is a 2D array
    dire is the direction of the random values strips to be generated
    s is the right border of random number generator, np.random.rand
    """
    if dire == 'h':
        return img + np.zeros((img.shape[0],img.shape[1]))+np.random.rand(img.shape[0],1)*s
    elif dire == 'v':
        return img + np.zeros((img.shape[0],img.shape[1]))+np.random.rand(1,img.shape[1])*s


def gradn(img, dire, max_value):
    """
    img is a 2D array
    dire is the direction of the gradient from 0 to max_value
    max_value is a close value to the maximum of this noise component
    """
    
    ny, nx = img.shape
    y, x = np.mgrid[:ny, :nx]
    
    if dire == 'ld': # left diagonal start
        a = img.shape[0]*img.shape[1]/max_value
        gradient = x * y / a
    elif dire == 'rd': # left diagonal start
        a = img.shape[0]*img.shape[1]/max_value
        gradient = (img.shape[0] - x  - 1)* y / a
        
    elif dire == 'rc': # right column  start
        a = img.shape[0]/max_value
        gradient = x / a
    elif dire == 'lc': # left column start
        a = img.shape[0]/max_value
        gradient = (img.shape[0] - x - 1)/ a
        
    elif dire == 'mb': # middle bottom row start
        a = img.shape[1]/max_value
        gradient = y / a
    elif dire == 'mt': # middle top row start
        a = img.shape[1]/max_value
        gradient = (img.shape[1] - y  -1) / a
    return img + gradient