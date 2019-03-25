import cv2 as cv
from util import Util
from gaussian_interpolation import GInterpolation as gi
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
import os
import math


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    def source_pixels(self, sigma_min=1.0, sigma_step=0.5, steps=10):
        """Search the source pixels position in the image"""
        sigma = 0
        smoothed = None

        for i in range(0, steps):
            if i == 0:
                smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma_min), sigma_min)
                sigma = sigma_min
            else:
                smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_step), sigma_step)
                sigma = math.sqrt(pow(sigma, 2) + pow(sigma_step, 2))
            maxima = Util.local_maxima(smoothed)
            if sigma >= 0.5:
                maxima = gi.optimize_maxima(smoothed, maxima, sigma)
            maxima = Util.sort_list(maxima, 1, True)

            print(sigma)

            q = 20
            last_d = 5
            score = []

            for maximum_a in maxima[:q]:
                distances = []
                for maximum_b in maxima:
                    distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
                distances.sort(reverse=False)
                s = 0
                for d in distances[1:last_d]:
                    s += d
                score.append(s)
            print(maxima[:q])
            print(score[:q])
            plt.plot(score[:q], label=str(sigma))
            plt.legend()
        plt.figure()

        return 0, 0


def main(index):
    curdir = os.getcwd()
    os.chdir("../Skymaps")
    skymaps = os.listdir(".")
    skymaps.sort()
    index -= 1
    sf = SourceFinder(skymaps[index])
    print(skymaps[index])
    pixels_coord = sf.source_pixels(1, 1, 10)
    # print(pixels_coord)
    #lon, lat = Util.from_pix_to_wcs(pixels_coord, sf.wcs)
    #print(lon, lat)
    #plt.show()
    os.chdir(curdir)
    # cv.waitKey(0)
    # cv.destroyAllWindows()


def main2(indexes, sigma):
    curdir = os.getcwd()
    os.chdir("../Skymaps")
    skymaps = os.listdir(".")
    skymaps.sort()

    for index in indexes:
        sf = SourceFinder(skymaps[index - 1])
        mat = cv.GaussianBlur(sf.matrix, Util.kernel_size(sigma), sigma)
        maxima = Util.local_maxima(mat)
        # print(maxima)
        maxima = gi.optimize_maxima(mat, maxima, sigma)
        # print(maxima)
        # cv.imshow(skymaps[index - 1] + str(index), Util.flip_mat(mat))
        Util.sort_list(maxima, 1, True)
        print(skymaps[index - 1], maxima)

        maxima = Util.project_list(maxima, 1)
        for i in range(0, len(maxima) - 1):
            maxima[i] = maxima[i] - maxima[i + 1]

        plt.plot(maxima, label=''+str(index))
        plt.legend()
        # cv.waitKey(0)
        # cv.destroyAllWindows()

    os.chdir(curdir)


# main
main(7)
main(8)
main(9)
main(10)
main(11)
main(12)
main(13)
# main(8)
# main(13)
plt.show()
# test main
# for s in range(1, 4):
#     main2([4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 1*s + 0.5)
#     plt.figure()
# plt.show()

