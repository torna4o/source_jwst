from photutils.datasets import make_test_psf_data
from photutils.psf import IntegratedGaussianPRF

def testing(n, shape, flu):

    psf_model = IntegratedGaussianPRF(flux=1, sigma=2.7 / 2.35)
    psf_shape = (25, 25)
    img, truth = make_test_psf_data(shape, psf_model, psf_shape,
                                           n, flux_range=flu,
                                           min_separation=10)
    truth.sort(['x', 'y'])
    return img, truth