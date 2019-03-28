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

    @staticmethod
    def best_candidate(isolatedness_values, k=1.5, quorum=0.3):
        candidates = []
        votes = []
        total_votes = 0

        for step in range(0, len(isolatedness_values[0]) - 2):
            values = Util.project_list(isolatedness_values, step + 2)
            sorted_values = sorted(values, reverse=True)
            if sorted_values[0] >= sorted_values[1] * k:
                if candidates.count(isolatedness_values[values.index(sorted_values[0])][1]) == 0:
                    candidates.append(isolatedness_values[values.index(sorted_values[0])][1])
                    votes.append(1)
                else:
                    votes[candidates.index(isolatedness_values[values.index(sorted_values[0])][1])] += 1
                total_votes += 1

        if not votes:
            return None

        sorted_votes = sorted(votes, reverse=True)
        if sorted_votes[0] >= (len(isolatedness_values[0]) - 2) * quorum:
            return candidates[votes.index(sorted_votes[0])][1], candidates[votes.index(sorted_votes[0])][0]
        else:
            return None

    def compute_isolatedness(self, sigmas, maxima_per_iter=3, furthest_index=2, dist_thresh=4):
        """Search the source pixels position in the image"""

        smoothed = None
        last_sigma = 0
        data = []

        for step in range(0, len(sigmas)):
            sigma = sigmas[step]
            if last_sigma == 0 or sigma <= last_sigma:
                smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma), sigma)
            else:
                sigma_to_add = math.sqrt(pow(sigma, 2) - pow(last_sigma, 2))
                smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_to_add), sigma_to_add)
            last_sigma = sigma
            maxima = Util.local_maxima(smoothed)
            if sigma >= 0.5:
                maxima = gi.optimize_maxima(smoothed, maxima, sigma)
            maxima = Util.sort_list(maxima, 1, True)

            for maximum_a in maxima[:maxima_per_iter]:
                distances = []
                for maximum_b in maxima:
                    distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
                distances.sort()
                isolatedness = 0
                for d in distances[1:furthest_index]:
                    isolatedness += d

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

        """
        for row in data:
            print(*row, sep="\t")
            plt.plot(sigmas, row[2:], label=str(int(row[0][0]))+", "+str(int(row[0][1])))
            plt.legend()

        plt.figure()
        """

        return data

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
    isolatedness_values = sf.compute_isolatedness([2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3])
    pixels_coord = sf.best_candidate(isolatedness_values)
    # print(pixels_coord)
    if pixels_coord:
        src_cor = np.round(Util.from_pix_to_wcs(pixels_coord, sf.wcs), 2)
        print(src_cor[0], src_cor[1])
    else:
        print("No source found")
    # plt.show()
    os.chdir(curdir)


# main
main(1)   # 145, 76
main(2)   # 150, 117
main(3)   # 86,  62
main(4)   # 116, 73
main(5)   # 121, 112
main(6)   # 155, 84
main(7)   # 142, 97
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
