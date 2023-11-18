from simu import testing as tt
from rna import phottest as pt
import noisy as ns
from matplotlib import pyplot as plt
import numpy as np

def rsq(truth, model):
    """
    truth is the real value of the simulated image
    model is what PSF Photometry estimates for x, y, flux
    plt below is PyPlot sublibrary of matplotlib
    """
    
    # The following prints R2 (coefficient of correlation) values of x,y,flux
    # between actual values and truth
    
    corx = np.corrcoef(truth['x'],model['x_fit'])
    cory = np.corrcoef(truth['y'],model['y_fit'])
    corf = np.corrcoef(truth['flux'],model['flux_fit'])
    print("Correlation between x and x_fit", corx[0,1]**2)
    print("Correlation between y and y_fit", cory[0,1]**2)
    print("Correlation between flux and flux_fit", corf[0,1]**2)
    
    fig, ax = plt.subplots(1, 3)
    ax[0].scatter(truth['x'], model['x_fit'])
    ax[0].set_title('Truth X value vs Model X Value')
    ax[0].set_ylabel('Model values')
    ax[1].scatter(truth['y'], model['y_fit'])
    ax[1].set_title('Truth Y value vs Model Y Value')
    ax[2].scatter(truth['flux'], model['flux_fit'])
    ax[2].set_title('Truth Flux value vs Model Flux Value')    
    plt.show()

def imgq(img_s, fitim_s, sp):
    if sp==True:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(img_s, cmap='Greys_r')
        ax[0].set_title('Sparse Truth')
        ax[1].imshow(fitim_s, cmap='Greys_r')
        ax[1].set_title('Sparse Model')
        plt.show()
    elif sp==False:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(img_s, cmap='Greys_r')
        ax[0].set_title('Crowded Truth')
        ax[1].imshow(fitim_s, cmap='Greys_r')
        ax[1].set_title('Crowded Model')
        plt.show()        


### Relatively crowded field simulation

img_c, truth_c = tt(75, (256,256),(110,150))
xyflux_c, fitim_c, resim_c = pt(img_c)
print("Truth for crowded:", truth_c)

# Sanity checks to see if PSFPhotometry retrieves the expected sources

rsq(truth_c, xyflux_c)
rsq(truth_s, xyflux_s)
imgq(img_c, fitim_c, sp=False)
imgq(img_s, fitim_s, sp=True)

"""
NOISING

Now we add various types of noises to both sparse and crowded images
and see if an algorithm robustly recover sources from either (or both)
images.

"""

# Gaussian, Salt and Pepper, Random vertical strips, top down gradient noises
# After some trials, the hard-coded numbers below chosen so as to have both
# false positives and negatives on estimated point sources

def noiser(img):
    img2, _ = ns.add_gauss(img, mean=1, std=1, speckle=False)
    img3, _ = ns.saltpep(img2, multi=4, threshold=0.95)
    img4 = ns.rstrip(img3, dire='v', s=10) 
    img5 = ns.gradn(img4, dire='mt', max_value=15)
    return img5

img_n = noiser(img_c)
xyflux_n, fitim_n, resim_n = pt(img_n)
imgq(img_c, fitim_n, sp=False)

"""
DENOISING

Primarily SSA, several ways to denoise this image will be tried. Starting from 
the available background2D removal techniques

"""

from astropy.stats import SigmaClip
from photutils.background import (MeanBackground, MedianBackground, ModeEstimatorBackground, 
                                MMMBackground, SExtractorBackground, BiweightLocationBackground)
from photutils.background import Background2D


sigma_clip = SigmaClip(sigma=3.0)

bkgm = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=MeanBackground())
bkgmd = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=MedianBackground())
bkgmo = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=ModeEstimatorBackground())
bkgmmm = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=MMMBackground())
bkgse = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=SExtractorBackground())
bkgbi = Background2D(img_n, (6, 6), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=BiweightLocationBackground())


# None were successful in retrieving everything and nothig more, hence the results were not reported. 

xyflux_nr, fitim_nr, resim_nr = pt(img_n - bkgmd.background - (img_n - bkgmd.background).min())



# Now going for SSA way

from pyts.decomposition import SingularSpectrumAnalysis

flat_data = (img_n - bkgmd.background).flatten(order='F') # 'F' better since the noise is vertical
flatt = np.array(range(img_n.shape[0]*img_n.shape[1]))
tossa = np.vstack((flatt, flat_data))
wind = 5 # 5 after several trials and errors
ssa = SingularSpectrumAnalysis(window_size=wind)
X_ssa = ssa.fit_transform(tossa)
xx = np.reshape(X_ssa, (2,wind,img_n.shape[0],img_n.shape[1]))

# This part is to create a source mask from the Second SSA component

from photutils.segmentation import detect_sources
dat = xx[1,1,:,:].T
bkg_estimator = BiweightLocationBackground()
bkg = Background2D(dat, (128, 128), filter_size=(3, 3), bkg_estimator=bkg_estimator)
threshold =  2 * bkg.background_rms

from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.utils import circular_footprint
kernel = make_2dgaussian_kernel(3.0, size=7) 
convolved_data = convolve(xx[1,1,:,:].T, kernel)
segment_map = detect_sources(convolved_data, threshold, npixels=4)
footprint = circular_footprint(radius=1)
mask2 = segment_map.make_source_mask(footprint=footprint)

# We now apply this mask to the noisy image to mask estimated source positions and have only background.

back3 = img_n*abs(mask2-1)
back3=back3.astype('float')
back3[back3 == 0.0] = np.NaN

# The following interfiller function fills the source-masked regions in the background with
# linear interpolation

def interpolate_column(col):
    nan_indices = np.isnan(col)
    valid_indices = np.where(~nan_indices)[0]
    if len(valid_indices) > 1:
        interp_func = np.interp(np.arange(len(col)), valid_indices, col[valid_indices])
        col[nan_indices] = interp_func[nan_indices]
    return col

def interfiller(img):
    rowsh, colsh = img.shape
    for i in range(colsh):
        img[:, i] = interpolate_column(img[:, i])
    return img

## Now we will have our estimated background.

dd2 = interfiller(back3)

# Checking the performance of the PSFPhotometry after removing this background.

xyflux_end4, fitim_end4, resim_end4 = pt(img_n - dd2)
rsq(truth_c, xyflux_end4)

fig, ax = plt.subplots(2,3)
ax[0,0].imshow(img_n, cmap='Greys_r')
ax[0,0].set_title('Noisy Image')
ax[0,1].imshow(xx[1,0,:,:].T, cmap='Greys_r')
ax[0,1].set_title('SSA 1st Component')
ax[0,2].imshow(xx[1,1,:,:].T, cmap='Greys_r')
ax[0,2].set_title('SSA 2nd Component')
ax[1,0].imshow(dd2, cmap='Greys_r')
ax[1,0].set_title('Estimated Background')
ax[1,1].imshow(fitim_end4, cmap='Greys_r')
ax[1,1].set_title('SSA Source Model')
ax[1,2].imshow(img_c, cmap='Greys_r')
ax[1,2].set_title('Original Image')
plt.tight_layout()
plt.show()