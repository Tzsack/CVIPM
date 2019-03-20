import cv2 as cv
import numpy as np
from util import Util
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    def source_pixels_routine(self, image, step):
        maxima = Util.local_maxima(image)
        maxima = Util.sort_list(maxima, 1, True)

        threshold = 0.1
        if maxima[0][1] - maxima[1][1] > threshold:
            x, y = Util.from_pix_to_wcs((maxima[0][0][1], maxima[0][0][0]), self.wcs)
            print("Step:", step, "-", maxima[0][0], x, y)
        else:
            print("Step:", step, "- No source found.")

    def source_pixels(self, sigma_min, sigma_step, steps):
        """Search the source pixels position in the image"""
        smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma_min), sigma_min)
        self.source_pixels_routine(smoothed, 0)

        for i in range(1, steps):
            smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_step), sigma_step)
            self.source_pixels_routine(smoothed, i)

        return 0, 0


def main(index):
    index -= 1
    sf = SourceFinder(skymaps[index])
    print(skymaps[index])
    pixels_coord = sf.source_pixels(1, 2, 5)
    # print(pixels_coord)
    # lon, lat = Util.from_pix_to_wcs(pixels_coord, sf.wcs)
    # print(lon, lat)
    plt.show()
    cv.waitKey(0)
    cv.destroyAllWindows()

def main2(indexes, sigma):
    for index in indexes:
        sf = SourceFinder(skymaps[index - 1])
        mat = cv.GaussianBlur(sf.matrix, Util.kernel_size(sigma), sigma)
        maxima = Util.local_maxima(mat)
        Util.sort_list(maxima, 1, True)
        print(maxima)
        plt.plot(Util.project_list(maxima, 1), label=''+str(index))
        plt.legend()
    plt.show()


skymaps = ["../Skymaps/01-skymap_220-0_46-0.fits",
           "../Skymaps/02-skymap_221-0_46-0.fits",
           "../Skymaps/03-skymap_221-0_46-0_4-0.fits",
           "../Skymaps/04-skymap_221-0_47-0.fits",
           "../Skymaps/05-skymap_221-0_47-0_0-0.fits",
           "../Skymaps/06-skymap_221-0_47-0_0-5.fits",
           "../Skymaps/07-skymap_221-0_47-0_0-8.fits",
           "../Skymaps/08-skymap_221-0_47-0_1-0.fits",
           "../Skymaps/09-skymap_221-5_45-5.fits",
           "../Skymaps/10-skymap_bgonly.fits",
           "../Skymaps/11-skymap_221-0_47-0_1-2.fits",
           "../Skymaps/12-skymap_221-0_47-0_2-0.fits",
           ]

# test main
main2([2, 6], 0.0001)
