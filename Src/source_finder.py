import math
import os

import cv2 as cv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.wcs import WCS

from gaussian_interpolation import GInterpolation as gi
from util import Util
import numpy as np


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    def source_pixels(self, sigmas):
        """Search the source pixels position in the image"""

        smoothed = None
        last_sigma = 0
        data = []

        for step in range(0, len(sigmas)):
            sigma = sigmas[step]
            if last_sigma == 0 or sigma <= last_sigma:
                smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma), sigma)
                last_sigma = sigma
            else:
                sigma_to_add = math.sqrt(pow(sigma, 2) - pow(last_sigma, 2))
                smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_to_add), sigma_to_add)
                last_sigma = sigma
            maxima = Util.local_maxima(smoothed)
            if sigma >= 0.5:
                maxima = gi.optimize_maxima(smoothed, maxima, sigma)
            maxima = Util.sort_list(maxima, 1, True)

            # Trying to create a synthetic version of the image with less background noise

            # cv.imshow(str(round(sigma, 2)), smoothed * 1/np.max(smoothed))

            # dist_threshold = 5

            # out = np.copy(smoothed)
            # for x in range(0, len(out)):
            #     print(x)
            #     for y in range(0, len(out[x])):
            #         out[x][y] = 0
            #         for m in maxima:
            #             out[x][y] += m[1] * gi.gaussian((x, y), m[0], sigma)
            #
            # cv.imshow("max " + str(round(sigma, 2)), out * 1/np.max(out))


            # Construction of a data matrix containing for each local maximum the corrispondig values of isolatedness
            # for increasing values of sigma

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
                        data[max_id][step + 2] = round(isolatedness, 2)
                        data[max_id][0] = maximum_a[0]
                        found = True
                if not found:
                    data.append([maximum_a[0]])
                    data[-1].append(maximum_a[0])
                    for j in range(0, len(sigmas)):
                        data[-1].append(float(0))
                    data[-1][step + 2] = round(isolatedness, 2)

        candidates = []
        votes = []
        k = 1.5
        quorum = 3/10
        total_votes = 0

        for step in range(0, len(sigmas)):
            values = Util.project_list(data, step + 2)
            sorted_values = sorted(values, reverse=True)
            if sorted_values[0] >= sorted_values[1] * k:
                if candidates.count(data[values.index(sorted_values[0])][1]) == 0:
                    candidates.append(data[values.index(sorted_values[0])][1])
                    votes.append(1)
                else:
                    votes[candidates.index(data[values.index(sorted_values[0])][1])] += 1
                total_votes += 1

        sorted_votes = sorted(votes, reverse=True)

        print(candidates)
        print(votes)

        if sorted_votes[0] >= len(sigmas) * quorum:
            print(candidates[votes.index(sorted_votes[0])])
        else:
            print("No source")

        for row in data:
            print(*row, sep="\t")
            plt.plot(sigmas, row[2:], label=str(int(row[0][0]))+", "+str(int(row[0][1])))
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
    pixels_coord = sf.source_pixels([2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3])
    # print(pixels_coord)
    #lon, lat = Util.from_pix_to_wcs(pixels_coord, sf.wcs)
    #print(lon, lat)
    #plt.show()
    os.chdir(curdir)


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

# main(1)   # 145, 76
# main(2)   # 150, 117
# main(3)   # 86,  62
# main(4)   # 116, 73
# main(5)   # 121, 112
# main(6)   # 155, 84
# main(7)   # 142, 97
main(8)   # 91, 131
main(9)   # 97, 93
main(10)  # 86, 121 (?)
main(11)  # 115, 80 (?)
main(12)  # 70, 107 (?)
main(13)  # Nowhere
plt.show()
# cv.waitKey(0)
# cv.destroyAllWindows()
# test main
# for s in range(1, 4):
#     main2([4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 1*s + 0.5)
#     plt.figure()
# plt.show()

