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

    def __init__(self, conf_file_name):
        """Constructor"""
        self.parameters = Util.read_list_from_json_file(conf_file_name)[0]
        self.matrix = None

    def compute_coords(self):
        skymaps_dir = self.parameters['dir']
        cur_dir = os.getcwd()
        os.chdir(skymaps_dir)
        coords = []
        for skymap in sorted(os.listdir('.')):
            if skymap == "computed_coordinates.json":
                continue
            # print(skymap)
            self.matrix = Util.from_fits_to_mat(skymap)
            isolatedness_values = self.compute_isolatedness()
            # self.visualize_data(isolatedness_values)  # UNCOMMENT TO SHOW GRAPHS
            src_pix_pos = self.best_candidate(isolatedness_values)
            if src_pix_pos:
                src_eq_pos = Util.from_pix_to_wcs(src_pix_pos, WCS(skymap))
                coords.append((float(src_eq_pos[0]), float(src_eq_pos[1])))
            else:
                coords.append(None)
        Util.write_list_to_json_file(coords, "computed_coordinates.json")
        # plt.show()  # UNCOMMENT TO SHOW GRAPHS
        os.chdir(cur_dir)
        return coords

    def best_candidate(self, isolatedness_values):
        candidates = []
        votes = []
        total_votes = 0

        for step in range(0, len(isolatedness_values[0]) - 2):
            values = Util.project_list(isolatedness_values, step + 2)
            if len(values) < 2:
                return None
            sorted_values = sorted(values, reverse=True)
            if sorted_values[0] >= sorted_values[1] * self.parameters['k']:
                if candidates.count(isolatedness_values[values.index(sorted_values[0])][1]) == 0:
                    candidates.append(isolatedness_values[values.index(sorted_values[0])][1])
                    votes.append(1)
                else:
                    votes[candidates.index(isolatedness_values[values.index(sorted_values[0])][1])] += 1
                total_votes += 1

        if not votes:
            return None

        # print(votes)
        sorted_votes = sorted(votes, reverse=True)
        if sorted_votes[0] >= (len(isolatedness_values[0]) - 2) * self.parameters['quorum']:
            return candidates[votes.index(sorted_votes[0])][1], candidates[votes.index(sorted_votes[0])][0]
        else:
            return None

    def compute_isolatedness(self):
        """Search the source pixels position in the image"""

        smoothed = None
        last_sigma = 0
        data = []

        for step in range(0, len(self.parameters['sigma_array'])):
            sigma = self.parameters['sigma_array'][step]
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

            for maximum_a in maxima[:self.parameters['maxima_per_iter']]:
                distances = []
                for maximum_b in maxima:
                    distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
                distances.sort()
                isolatedness = 0
                for d in distances[1:self.parameters['furthest_index']]:
                    isolatedness += d

                found = False
                for max_id in range(0, len(data)):
                    if Util.distance_eu(maximum_a[0], data[max_id][0]) <= self.parameters['dist_thresh']:
                        data[max_id][step + 2] = round(isolatedness, 2)
                        data[max_id][0] = maximum_a[0]
                        found = True
                if not found:
                    data.append([maximum_a[0]])
                    data[-1].append(maximum_a[0])
                    for j in range(0, len(self.parameters['sigma_array'])):
                        data[-1].append(float(0))
                    data[-1][step + 2] = round(isolatedness, 2)

        return data

    def visualize_data(self, data):
        """Visualize data"""
        for row in data:
            print(*row, sep="\t")
            plt.plot(self.parameters['sigma_array'], row[2:], label=str(int(row[0][0]))+", "+str(int(row[0][1])))
            plt.legend()
        plt.figure()


def main():
    sf = SourceFinder("conf.json")
    coords = sf.compute_coords()
    print(coords)


# main()
