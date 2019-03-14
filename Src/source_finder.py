import cv2 as cv
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


class SourceFinder:

    def __init__(self, file_name):
        """Constructor"""
        self.file_name = file_name
        self.matrix = None

    @staticmethod
    def kernel_size(sigma):
        """Determine kernel size according to sigma"""
        size = int(sigma * 6)
        if size >= 5:
            if size % 2 == 0:
                size += 1
            return size, size
        else:
            return 5, 5

    @staticmethod
    def best_candidate(candidates):
        """Return most voted candidate in candidates"""
        result = candidates[0]
        max_votes = 0

        for source in candidates:
            if candidates.count(source) > max_votes:
                max_votes = candidates.count(source)
                result = source

        return result

    def load_fits_file(self):
        """Save the fits image in the matrix attribute"""
        file = fits.open(self.file_name)
        self.matrix = np.array(file[0].data, np.float64)
        file.close()

    def source_pixels(self, min_sigma=1.0, max_sigma=3.0, steps=10):
        """Search the source pixels position in the image"""
        candidates = []

        for i in range(0, steps):
            sigma = min_sigma + ((max_sigma - min_sigma)/steps)*i
            out = cv.GaussianBlur(self.matrix, self.kernel_size(sigma), sigma)
            # cv.imshow("a"+str(i), out)
            candidates.append(np.unravel_index(np.argmax(out, axis=None), out.shape))

        best_candidate = self.best_candidate(candidates)
        return best_candidate[1], best_candidate[0]

    def source_world(self, source_pixels):
        """Return world coordinates for the given pixel position source_pixels"""
        w = WCS(self.file_name)
        return w.all_pix2world(source_pixels[0], source_pixels[1], 0)


fits_file = "../Skymaps/skymap_221_47.fits"
sf = SourceFinder(fits_file)
sf.load_fits_file()
pixels_coord = sf.source_pixels(1.5, 3.5, 10)
lon, lat = sf.source_world(pixels_coord)
print(lon, lat)

cv.waitKey(0)
cv.destroyAllWindows()
