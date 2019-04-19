import math
import os

import cv2 as cv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.wcs import WCS

from gaussian_interpolation import GInterpolation as gi
from util import Util
import numpy as np
from point_data import PointData
from point_data import PointDataList


class SourceFinder:

    def __init__(self, conf_file_name):
        """Constructor"""
        self.parameters = Util.read_list_from_json_file(conf_file_name)
        self.matrix = None

    def compute_coords(self):
        skymaps_dir = self.parameters['dir']
        cur_dir = os.getcwd()
        os.chdir(skymaps_dir)
        coords = []
        os.mkdir('Measures')
        for skymap in sorted(os.listdir('.')):
            if skymap.endswith('.fits'):
                # print(skymap)
                self.matrix = Util.from_fits_to_mat(skymap)
                measures = self.compute_measures(skymap.split('.')[0])
                self.visualize_data(measures, WCS(skymap))  # UNCOMMENT TO SHOW GRAPHS
                src_pix_pos = self.best_candidate(measures)
                if src_pix_pos:
                    src_eq_pos = Util.from_pix_to_wcs(src_pix_pos, WCS(skymap))
                    coords.append((float(src_eq_pos[0]), float(src_eq_pos[1])))
                else:
                    coords.append(None)
        Util.write_list_to_json_file(coords, "computed_coordinates.json")
        os.chdir(cur_dir)
        return coords

    def best_candidate(self, measures):
        candidates = []
        votes = []
        total_votes = 0

        intensity_pdl = measures['intensity']
        for step in range(0, len(intensity_pdl.points)):
            values = intensity_pdl.project_list(step)
            if len(values) < 2:
                return None
            sorted_values = sorted(values, reverse=True)
            if sorted_values[0] >= sorted_values[1] * self.parameters['k']:
                if candidates.count(intensity_pdl.points[values.index(sorted_values[0])].original) == 0:
                    candidates.append(intensity_pdl.points[values.index(sorted_values[0])].original)
                    votes.append(1)
                else:
                    votes[candidates.index(intensity_pdl.points[values.index(sorted_values[0])].original)] += 1
                total_votes += 1

        if not votes:
            return None

        # print(votes)
        sorted_votes = sorted(votes, reverse=True)
        if sorted_votes[0] >= len(intensity_pdl.points[0].values) * self.parameters['quorum']:
            return candidates[votes.index(sorted_votes[0])][1], candidates[votes.index(sorted_votes[0])][0]
        else:
            return None

    def visualize_data(self, measures, wcs):
        """Visualize data"""
        for key in measures.keys():
            measure = measures[key]
            print(measure.to_string())
            for point in measure.points:
                eq = Util.from_pix_to_wcs((point.original[1], point.original[0]), wcs)
                plt.plot(self.parameters['sigma_array'], point.values, label=str(round(float(eq[0]), 2))+", "+str(round(float(eq[1]), 2)))
                plt.title(key)
                plt.legend()
            plt.figure()

    def update_isolatedness(self, maxima, step, point_data_list):
        for maximum_a in maxima[:self.parameters['maxima_per_iter']]:
            distances = []
            for maximum_b in maxima:
                distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
            distances.sort()
            isolatedness = 0
            for d in distances[1:self.parameters['furthest_index']]:
                isolatedness += d
            point_data_list.set(maximum_a[0], isolatedness, step)

    def update_intensity(self, maxima, step, point_data_list):
        for maximum_a in maxima[:self.parameters['maxima_per_iter']]:
            point_data_list.set(maximum_a[0], maximum_a[1], step)

    @staticmethod
    def sorted_maxima(smoothed, sigma):
        """"""
        maxima = Util.local_maxima(smoothed)
        if sigma >= 0.5:
            maxima = gi.optimize_maxima(smoothed, maxima, sigma)
        maxima = Util.sort_list(maxima, 1, True)
        return maxima

    def smooth_image(self, sigma, last_sigma, smoothed):
        if last_sigma == 0 or sigma <= last_sigma:
            smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma), sigma)
        else:
            sigma_to_add = math.sqrt(pow(sigma, 2) - pow(last_sigma, 2))
            smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_to_add), sigma_to_add)
        return smoothed

    def compute_measures(self, skymap):
        """"""
        smoothed = None
        last_sigma = 0
        measures = {'isolatedness': PointDataList(self.parameters['dist_thresh']),
                    'intensity': PointDataList(self.parameters['dist_thresh'])}

        for step in range(0, len(self.parameters['sigma_array'])):
            sigma = self.parameters['sigma_array'][step]
            smoothed = self.smooth_image(sigma, last_sigma, smoothed)
            last_sigma = sigma
            maxima = self.sorted_maxima(smoothed, sigma)
            self.update_isolatedness(maxima, step, measures['isolatedness'])
            self.update_intensity(maxima, step, measures['intensity'])

        for key in measures.keys():
            Util.write_list_to_json_file(measures[key].to_list(), 'Measures/' + skymap + '_' + key + '.json')

        return measures


def main():
    sf = SourceFinder("conf.json")
    coords = sf.compute_coords()
    print(coords)
    plt.show()


# main()
