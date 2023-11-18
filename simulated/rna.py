from photutils.detection import DAOStarFinder
from photutils.psf import IntegratedGaussianPRF, PSFPhotometry

def phottest(img):
    psf_model = IntegratedGaussianPRF(flux=1, sigma=2.7 / 2.35)
    fit_shape = (3, 3)
    finder = DAOStarFinder(4.0, 2.0)
    psfphot = PSFPhotometry(psf_model, fit_shape, finder=finder,
                            aperture_radius=4)
    phot = psfphot(img)
    xyflux = phot[('x_fit', 'y_fit', 'flux_fit')]
    xyflux.sort(['x_fit','y_fit'])
    print(xyflux)
    fitim = psfphot.make_model_image(img.shape,  (25, 25))
    resim = psfphot.make_residual_image(img, (25, 25))
    return xyflux, fitim, resim