from astropy.io import fits
file = fits.open("../Skymaps/skymap_221_47_1.fits")
matrix = file[0].data
for i in matrix:
    for j in i:
        if j >= 4:
            print(j)
