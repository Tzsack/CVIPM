import cv2 as cv
import numpy as np
from util import Util
from astropy.io import fits
from astropy.wcs import WCS


class SourceFinder:

    def __init__(self, file_name):
        """Constructor"""
        self.file_name = file_name
        self.matrix = None

    def load_fits_file(self):
        """Save the fits image in the matrix attribute"""
        file = fits.open(self.file_name)
        self.matrix = np.array(file[0].data, np.float64)
        file.close()

    def source_pixels(self):
        """Search the source pixels position in the image"""
        sigma = 3.0
        temp = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma), sigma)
        return 0, 0

        # Thresholding

    def source_world(self, source_pixels):
        """Return world coordinates for the given pixel position source_pixels"""
        w = WCS(self.file_name)
        return w.all_pix2world(source_pixels[0], source_pixels[1], 0)

# Test main


fits_file = "../Skymaps/skymap_221-0_47-0_0-8_no_source.fits"
sf = SourceFinder(fits_file)
sf.load_fits_file()
pixels_coord = sf.source_pixels()
lon, lat = sf.source_world(pixels_coord)
print(lon, lat)

cv.waitKey(0)
cv.destroyAllWindows()
