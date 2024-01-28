from astropy.io import fits # FITS file opener, knowledge extractor and manipulator
from astropy.wcs import WCS # World Coordinate System wrapper and manipulator

import os # Operating System library, locating files etc.
import matplotlib.pyplot as plt # MATlab like PLOTting LIBrary as PYthon PLoTs

from matplotlib.patches import Rectangle as re # Patching Matplotlib PLOTS with Rectangles
import numpy as np #Numerical Python, basic matrix operations etc





def opener(i2d):
    with fits.open(i2d) as hdu_fits: 
        image_data = hdu_fits['SCI'].data # this is a numpy nd array
        wcs_helix = WCS(hdu_fits[1].header) # WCS object from FITS
    # The following plots according to the WCS information from the FITS
    ax = plt.subplot(projection=wcs_helix) # subplot with specific projection
    ax.imshow(image_data, origin='lower', cmap='Greys_r')
    rect = re((430, 20), 590, 980, linewidth=2, edgecolor='y', facecolor='none')
    ax.add_patch(rect)
    overlay = ax.get_coords_overlay('fk5') #Overlaying ICRS coordinates on counter axes
    overlay.grid(color='red', ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')
    plt.show()
    headheader = {"PI:": hdu_fits[0].header['PI_NAME'],
                  "OBS ID:": hdu_fits[0].header['OBS_ID'],
                  "OBS DATE:": hdu_fits[0].header['DATE-OBS'],
                  "Effect. Exp:": hdu_fits[0].header['EFFEXPTM'],
                  "Guide Star ID:": hdu_fits[0].header['GDSTARID'],
                  "Guide Star RA:": hdu_fits[0].header['GS_RA'],
                  "Guide Star DEC:": hdu_fits[0].header['GS_DEC']}
    return image_data, wcs_helix, headheader


i2d = [i for i in os.listdir(os.getcwd()) if i[-8:] == "i2d.fits"]

img1, wcs1, head1 = opener(i2d[0]) #vmin 5.9 ~ 6.0 for graph in img
img2, wcs2, head2 = opener(i2d[1]) #vmin 5.9 ~ 6.3 for graph in img
img3, wcs3, head3 = opener(i2d[2]) #vmin 10.7 ~ 11 for graph in img
img4, wcs4, head4 = opener(i2d[3]) #vmin 10 ~ 11 for graph in img
img5, wcs5, head5 = opener(i2d[4]) #vmin 4.3 ~ 4.4 for graph in img
img6, wcs6, head6 = opener(i2d[5]) #vmin 10 ~ 11 for graph in img
img7, wcs7, head7 = opener(i2d[6]) #vmin 4.3 ~ 4.4 for graph in img


# This thing will print dictionaries beautifully
def pp(head):
    for i in head:
        print(i, head[i])


from astroquery.vizier import Vizier # To query catalogs in VizieR
from astropy.coordinates import SkyCoord # Celestial coordinate repr. interface
from astropy.coordinates import ICRS # A common choice of coordinate ref. system
import astropy.units as u # ensuring standard unit to represent custom lengths

def unquery(wcs_helix):
    # Parameters, coordinated retrieved here
    coord = SkyCoord(ra=wcs_helix.wcs.crval[0], dec=wcs_helix.wcs.crval[1], unit="deg", frame="icrs")
    exp_coeff = 0.3 # A multiplier of width and height for querying
    width = u.Quantity(0.1*exp_coeff, u.deg)
    height = u.Quantity(0.1*exp_coeff, u.deg)
    # Querying here
    r2 = Vizier.query_region(coord, width=width, catalog = 'II/363/unwise')
    Vizier.ROW_LIMIT = -1 # alters the 50 rows default limit
    print(r2)
    ralist = r2['II/363/unwise']["RAJ2000"].tolist()
    declist = r2['II/363/unwise']["DEJ2000"].tolist()
    names = r2['II/363/unwise']["objID"].tolist()
    return  ralist, declist, names

def plotter(image_data, wcs_helix):
    # Retrieving the RA, DEC, name of the objects appeared in the query of UNWISE catalog
    ralist, declist, names = unquery(wcs_helix)
    # Plotting the query results overlayed on the image
    ax = plt.subplot(projection=wcs_helix) # subplot with specific projection
    ax.imshow(image_data, origin='lower', cmap='Greys_r', vmin = np.median(image_data)-0.5, vmax = np.median(image_data)+1)
    overlay = ax.get_coords_overlay('fk5') #Overlaying ICRS coordinates on counter axes
    overlay.grid(color='red', ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')
    ax.scatter(ralist,declist,transform=ax.get_transform('fk5'), s=15, edgecolor='red')
    for (xi, yi, ni) in zip(ralist, declist, names):
        plt.text(xi, yi, ni, va='bottom', ha='center', fontsize=5, color= 'yellow', transform=ax.get_transform('fk5'))
    plt.show()
    return ralist, declist, names


ra1, dec1, nam1 = plotter(img1, wcs1)
ra2, dec2, nam2 = plotter(img2, wcs2)
ra3, dec3, nam3 = plotter(img3, wcs3)
ra4, dec4, nam4 = plotter(img4, wcs4)
ra5, dec5, nam5 = plotter(img5, wcs5)
ra6, dec6, nam6 = plotter(img6, wcs6)
ra7, dec7, nam7 = plotter(img7, wcs7)


seg = [i for i in os.listdir(os.getcwd()) if i[-8:] == "egm.fits"]

seg1, _, _ = opener(seg[0])
seg2, _, _ = opener(seg[1])
seg3, _, _ = opener(seg[2])
seg4, _, _ = opener(seg[3])
seg5, _, _ = opener(seg[4])
seg6, _, _ = opener(seg[5])
seg7, _, _ = opener(seg[6])

# The following code snippet creates the histogram image I draw in the manuscript


fig, ax = plt.subplots(2,4)
ax[0,0].hist(img1[20:20+980,430:430+590].flatten(), range=(4,15),bins=100)
ax[0,0].axvline(np.mean(img1[20:20+980,430:430+590]) + np.std(img1[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[0,0].axvline(np.percentile(img1[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[0,1].hist(img2[20:20+980,430:430+590].flatten(), range=(5,7),bins=100)
ax[0,1].axvline(np.mean(img2[20:20+980,430:430+590]) + np.std(img2[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[0,1].axvline(np.percentile(img2[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[0,2].hist(img3[20:20+980,430:430+590].flatten(), range=(9.8,12),bins=100)
ax[0,2].axvline(np.mean(img3[20:20+980,430:430+590]) + np.std(img3[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[0,2].axvline(np.percentile(img3[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[0,3].hist(img4[20:20+980,430:430+590].flatten(), range=(9,11.5),bins=100)
ax[0,3].axvline(np.mean(img4[20:20+980,430:430+590]) + np.std(img4[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[0,3].axvline(np.percentile(img4[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[1,0].hist(img5[20:20+980,430:430+590].flatten(), range=(5,7),bins=100)
ax[1,0].axvline(np.mean(img5[20:20+980,430:430+590]) + np.std(img5[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[1,0].axvline(np.percentile(img5[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[1,1].hist(img6[20:20+980,430:430+590].flatten(), range=(4.2,5),bins=100)
ax[1,1].axvline(np.mean(img6[20:20+980,430:430+590]) + np.std(img6[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[1,1].axvline(np.percentile(img6[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[1,2].hist(img7[20:20+980,430:430+590].flatten(), range=(6, 7.5),bins=100)
ax[1,2].axvline(np.mean(img7[20:20+980,430:430+590]) + np.std(img7[20:20+980,430:430+590])*2, linestyle='dashed', color='r')
ax[1,2].axvline(np.percentile(img7[20:20+980,430:430+590].flatten(), 95.4), linestyle='dashed', color='k')
ax[0,0].set_title("Img1 Histogram")
ax[0,0].set_ylabel("Count")
ax[1,0].set_ylabel("Count")
ax[0,1].set_ylabel("Count")
ax[1,1].set_ylabel("Count")
ax[0,2].set_ylabel("Count")
ax[1,2].set_ylabel("Count")
ax[0,3].set_ylabel("Count")
ax[1,0].set_xlabel("Digital Number")
ax[1,1].set_xlabel("Digital Number")
ax[1,2].set_xlabel("Digital Number")
ax[0,0].set_xlabel("Digital Number")
ax[0,1].set_xlabel("Digital Number")
ax[0,2].set_xlabel("Digital Number")
ax[0,3].set_xlabel("Digital Number")
ax[0,1].set_title("Img2 Histogram")
ax[0,2].set_title("Img3 Histogram")
ax[0,3].set_title("Img4 Histogram")
ax[1,0].set_title("Img5 Histogram")
ax[1,1].set_title("Img6 Histogram")
ax[1,2].set_title("Img7 Histogram")
fig.delaxes(ax[1,3])
plt.tight_layout()
plt.show()


from matplotlib.collections import EllipseCollection # For drawing circles around unWISE source locations


# Creating derivative images

diff1 = np.gradient(img1[20:20+980,430:430+590], axis=0)
diff2 = np.gradient(img2[20:20+980,430:430+590], axis=0)
diff3 = np.gradient(img3[20:20+980,430:430+590], axis=0)
diff4 = np.gradient(img4[20:20+980,430:430+590], axis=0)
diff5 = np.gradient(img5[20:20+980,430:430+590], axis=0)
diff6 = np.gradient(img6[20:20+980,430:430+590], axis=0)
diff7 = np.gradient(img7[20:20+980,430:430+590], axis=0)




from scipy.stats import kurtosis
from scipy.stats import skew

# To retrieve image statistics related parameters
def summ(kk):
    print("Min-Max: ", kk.min(), kk.max())
    print("Mean: ", kk.mean())
    print("Median: ", np.median(kk))
    print("variance  : ", np.var(kk))
    print("skew : ",skew(kk.flatten()))
    print("kurt : ",kurtosis(kk.flatten()))


# The following is the 4-image figure plotter, you can change the ratios and retrieve other figures

def dplot(img, seg, diff, wcs, ra, dec):
    ssr1 = np.zeros([seg2.shape[0], seg2.shape[1]])
    ssr1[20:20+980,430:430+590] = diff
    offsets = list(zip(ra, dec))
    fig = plt.figure()
    ax = fig.add_subplot(1, 4, 1, projection=wcs)
    rec = re((430, 20), 590, 980, linewidth=2, edgecolor='y', facecolor='none', transform=ax.transData)
    ax.imshow(img, origin='lower', cmap='Greys_r', vmin = np.percentile(img, 80.0), vmax = np.percentile(img, 99.0)).set_clip_path(rec)
    ax.add_collection(EllipseCollection(widths=np.ones([len(ra)])*36, heights=np.ones([len(ra)])*36, angles=0, units='xy',
                                       facecolors = 'None',edgecolors='red', offsets=offsets, transOffset = ax.get_transform('fk5')))        
    ax.set_title("Original Image")
    ax0 = fig.add_subplot(1, 4, 2, projection=wcs)
    rec = re((430, 20), 590, 980, linewidth=2, edgecolor='y', facecolor='none', transform=ax0.transData)
    ax0.imshow(img, origin='lower', cmap='Greys_r', vmin = np.percentile(img, 99.7), vmax = np.percentile(img, 99.71)).set_clip_path(rec)
    ax0.add_collection(EllipseCollection(widths=np.ones([len(ra)])*36, heights=np.ones([len(ra)])*36, angles=0, units='xy',
                                       facecolors = 'None',edgecolors='red', offsets=offsets, transOffset = ax0.get_transform('fk5')))      
    ax0.set_title("Original Image, top 0.3%")
    ax1 = fig.add_subplot(1, 4, 3, projection=wcs)
    rec = re((430, 20), 590, 980, linewidth=2, edgecolor='y', facecolor='none', transform=ax1.transData)
    ax1.imshow(ssr1, origin='lower', cmap='Greys_r', vmin = np.percentile(ssr1, 99.7), vmax = np.percentile(ssr1, 99.71)).set_clip_path(rec)    
    ax1.add_collection(EllipseCollection(widths=np.ones([len(ra)])*36, heights=np.ones([len(ra)])*36, angles=0, units='xy',
                                       facecolors = 'None',edgecolors='red', offsets=offsets, transOffset = ax1.get_transform('fk5')))    
    ax1.set_title("Gradient Image, top 0.3%")
    ax2 = fig.add_subplot(1, 4, 4, projection=wcs)
    rec = re((430, 20), 590, 980, linewidth=2, edgecolor='y', facecolor='none', transform=ax2.transData)
    ax2.imshow(seg, origin='lower', cmap='Greys_r', vmin = 0, vmax = 0.01).set_clip_path(rec)
    ax2.add_collection(EllipseCollection(widths=np.ones([len(ra)])*36, heights=np.ones([len(ra)])*36, angles=0, units='xy',
                                       facecolors = 'None',edgecolors='red', offsets=offsets, transOffset = ax2.get_transform('fk5'))) 
    ax2.set_title("Original Segmentation Image")
    ax2.set_xlabel("RA")
    ax1.set_xlabel("RA")
    ax0.set_xlabel("RA")
    ax.set_xlabel("RA")
    ax2.set_ylabel("DEC")
    ax1.set_ylabel("DEC")
    ax.set_ylabel("DEC")
    plt.show()



dplot(img1, seg1, diff1, wcs1, ra1, dec1)
dplot(img2, seg2, diff2, wcs2, ra2, dec2)
dplot(img3, seg3, diff3, wcs3, ra3, dec3)
dplot(img4, seg4, diff4, wcs4, ra4, dec4)
dplot(img5, seg5, diff5, wcs5, ra5, dec5)
dplot(img6, seg6, diff6, wcs6, ra6, dec6)
dplot(img7, seg7, diff7, wcs7, ra7, dec7)

summ(img1[20:20+980,430:430+590])
summ(img2[20:20+980,430:430+590])
summ(img3[20:20+980,430:430+590])
summ(img4[20:20+980,430:430+590])
summ(img5[20:20+980,430:430+590])
summ(img6[20:20+980,430:430+590])
summ(img7[20:20+980,430:430+590])


