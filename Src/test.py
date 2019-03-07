from astropy.io import fits
file = fits.open("../Skymaps/skymap_221_46.fits")
matrix = file[0].data
for v in matrix:
    print(v)

