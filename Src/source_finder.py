import math
import os

import cv2 as cv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.wcs import WCS

from gaussian_interpolation import GInterpolation as gi
from util import Util


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    def source_pixels(self, sigma_min=1.0, sigma_step=0.5, steps=10):
        """Search the source pixels position in the image"""
        sigma = 0
        smoothed = None

        data = []

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

            last_max = 3
            last_d = 2
            dist_thresh = 4

            for maximum_a in maxima[:last_max]:
                distances = []
                for maximum_b in maxima:
                    distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
                distances.sort(reverse=False)
                isolatedness = 0
                for d in distances[1:last_d]:
                    isolatedness += d

                # data matrix construction

                found = False
                for max_id in range(0, len(data)):
                    if Util.distance_eu(maximum_a[0], data[max_id][0]) <= dist_thresh:
                        data[max_id][i + 1] = round(isolatedness, 2)
                        data[max_id][0] = maximum_a[0]
                        found = True
                if not found:
                    data.append([maximum_a[0]])
                    for j in range(0, steps):
                        data[-1].append(float(0))
                    data[-1][i + 1] = round(isolatedness, 2)

        k = 1 + 1/3
        candidates = []
        votes = []

        for col in range(1, steps + 1):
            values = Util.project_list(data, col)
            first = values.index(max(values))
            values.remove(values[first])
            second = values.index(max(values))
            if first >= second * k:
                if candidates.count(data[first][0]) == 0:
                    candidates.append(data[first][0])
                    votes.append(1)
                else:
                    votes[candidates.index(data[first][0])] += 1

        if candidates:
            print(str(candidates[votes.index(max(votes))]))
        else:
            print("No source found.")

        for row in data:
            print(*row, sep="\t")
            plt.plot(row[1:-1], label=str(row[0]))
            plt.legend()

        plt.figure()

        return 0, 0


"""
                    SIGMA
            1   2   3   4   5   6 ...
 M (x, y)   10  11
 A (x, y) 
 X ...
 I 
 M 
 A
 
 I
 D
"""


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

main(1)
# main(2)
# main(3)
# main(4)
# main(5)
# main(6)
# main(7)
main(8)
main(9)
main(10)
main(11)
main(12)
main(13)
plt.show()
# test main
# for s in range(1, 4):
#     main2([4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 1*s + 0.5)
#     plt.figure()
# plt.show()

