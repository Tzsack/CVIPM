import cv2 as cv
import numpy as np
from astropy.io import fits


class SourceFinder:

    def __init__(self):
        """"""
        self.matrix = None

    @staticmethod
    def kernel_size(sigma):
        size = int(sigma * 6)
        if size % 2 != 1:
            size += 1
        if size >= 5:
            return size, size
        else:
            return 5, 5

    def load_fits_file(self, filename):
        """"""
        file = fits.open(filename)
        self.matrix = np.array(file[0].data, np.float64)

    def search(self, min_sigma=1, max_sigma=3, steps=10):
        """"""
        candidates = []
        for i in range(0, steps):
            sigma = min_sigma + (max_sigma - min_sigma)/steps*i
            out = cv.GaussianBlur(self.matrix, self.kernel_size(sigma), sigma)
            cv.imshow("a"+str(i), out)
            candidates.append(np.unravel_index(np.argmax(out, axis=None), out.shape))

        result = candidates[0]
        max_votes = 0

        for source in candidates:
            if candidates.count(source) > max_votes:
                max_votes = candidates.count(source)
                result = source

        return result


sf = SourceFinder()
sf.load_fits_file("../Skymaps/skymap_221_47_1.fits")
print(sf.search(0.1, 6, 20))

cv.waitKey(0)
cv.destroyAllWindows()
