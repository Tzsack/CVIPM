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

    @staticmethod
    def best_candidate(candidates):
        """Return most voted candidate in candidates"""

        labels = [0]
        threshold = 1

        for i in range(1, len(candidates)):
            if Util.distance4(candidates[i], candidates[i - 1]) <= threshold:
                labels.append(labels[i-1])
            else:
                labels.append(labels[i-1]+1)

        winner = labels[0]
        max_votes = 0

        for label in labels:
            if labels.count(label) > max_votes:
                max_votes = labels.count(label)
                winner = label

        average_x = 0
        average_y = 0

        for i in range(0, len(candidates)):
            print(candidates[i], labels[i])
            if labels[i] == winner:
                average_x += candidates[i][0]
                average_y += candidates[i][1]

        average_x /= labels.count(winner)
        average_y /= labels.count(winner)

        return average_x, average_y

    def load_fits_file(self):
        """Save the fits image in the matrix attribute"""
        file = fits.open(self.file_name)
        self.matrix = np.array(file[0].data, np.float64)
        file.close()

    # def source_pixels(self, min_sigma=1.0, max_sigma=3.0, steps=10):
    #     """Search the source pixels position in the image."""
    #     candidates = []
    #
    #     for i in range(0, steps):
    #         sigma = min_sigma + ((max_sigma - min_sigma)/steps)*i
    #         out = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma), sigma)
    #         # cv.imshow("a"+str(i), out)
    #         candidates.append(np.unravel_index(np.argmax(out, axis=None), out.shape))
    #
    #     best_candidate = self.best_candidate(candidates)
    #     return best_candidate[1], best_candidate[0]

    def source_pixels(self, min_sigma=1.0, step_sigma=0.5, steps=10):
        """Search the source pixels position in the image"""
        candidates = []

        out = cv.GaussianBlur(self.matrix, Util.kernel_size(min_sigma), min_sigma)
        candidates.append(np.unravel_index(np.argmax(out, axis=None), out.shape))
        cv.imshow("a"+str(0), 1/(out[candidates[0][0]][candidates[0][1]])*out)

        for i in range(0, steps):
            out = cv.GaussianBlur(out, Util.kernel_size(step_sigma), step_sigma)
            candidates.append(np.unravel_index(np.argmax(out, axis=None), out.shape))
            cv.imshow("a"+str(i+1), 1/(out[candidates[i][0]][candidates[i][1]])*out)

        best_candidate = self.best_candidate(candidates)
        return best_candidate[1], best_candidate[0]

    def source_world(self, source_pixels):
        """Return world coordinates for the given pixel position source_pixels"""
        w = WCS(self.file_name)
        return w.all_pix2world(source_pixels[0], source_pixels[1], 0)


def main():
    fits_file = "../Skymaps/skymap_221-0_47-0_0-8_no_source.fits"
    sf = SourceFinder(fits_file)
    sf.load_fits_file()
    pixels_coord = sf.source_pixels(0.0001, 1.5, 10)
    lon, lat = sf.source_world(pixels_coord)
    print(lon, lat)

    cv.waitKey(0)
    cv.destroyAllWindows()


# Test main
# main()