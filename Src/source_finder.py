import math
import os
import shutil

import cv2 as cv
import matplotlib.pyplot as plt
from astropy.wcs import WCS

from gaussian_interpolation import GInterpolation as gi
from util import Util
from point_data import PointDataList
import numpy as np


class SourceFinder:

    def __init__(self, conf_file_name):
        """Constructor"""
        self.parameters = Util.read_list_from_json_file(conf_file_name)
        self.original = None
        self.images = []

    def compute_coords(self):
        """"""
        skymaps_dir = self.parameters['dir']
        cur_dir = os.getcwd()
        os.chdir(skymaps_dir)
        coords = []
        if os.listdir('.').count('Measures') != 0:
            shutil.rmtree('Measures')
        os.mkdir('Measures')

        for skymap in sorted(os.listdir('.')):
            if skymap.endswith('.fits'):
                print(skymap)
                self.original = Util.from_fits_to_mat(skymap)
                measures = self.compute_measures(skymap.split('.')[0])
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
        candidates = PointDataList(self.parameters['election_dist_thresh'])
        no_source = 0

        for key in measures.keys():
            pdl = measures[key]
            for step in range(0, len(pdl.points[0].values)):
                values = pdl.project_list(step)
                sorted_values = sorted(values, reverse=True)
                votes = self.parameters[key]['vote_weight']
                # print(key, step, sorted_values[0], '/', sorted_values[1], '=', sorted_values[0] / sorted_values[1])
                if len(sorted_values) < 2 or \
                        sorted_values[0] >= sorted_values[1] * self.parameters[key]['vote_threshold']:
                    coords = pdl.points[values.index(sorted_values[0])].original
                    point = candidates.get_point(coords)
                    if point:
                        votes += point.values[0]
                    candidates.set(coords, votes, 0)
                else:
                    no_source += votes

        print(candidates.to_string())
        print('no_source', no_source)
        print('')

        winner = None
        winner_votes = no_source
        for point in candidates.points:
            if point.values[0] > winner_votes:
                winner_votes = point.values[0]
                winner = point.original

        return winner

    def update_isolatedness(self, maxima, step, point_data_list):
        for maximum_a in maxima[:self.parameters['isolatedness']['maxima_per_iter']]:
            distances = []
            for maximum_b in maxima:
                distances.append(Util.distance_eu(maximum_a[0], maximum_b[0]))
            distances.sort()
            isolatedness = 0
            for d in distances[1:self.parameters['isolatedness']['furthest_index']]:
                isolatedness += d
            point_data_list.set(maximum_a[0], isolatedness, step)

    def update_intensity(self, maxima, step, point_data_list):
        for maximum_a in maxima[:self.parameters['intensity']['maxima_per_iter']]:
            val = maximum_a[1]
            point_data_list.set(maximum_a[0], val, step)

    @staticmethod
    def sorted_maxima(smoothed, sigma):
        """"""
        maxima = Util.local_maxima(smoothed)
        # if sigma >= 0.5:
        #    maxima = gi.optimize_maxima(smoothed, maxima, sigma)
        maxima = Util.sort_list(maxima, 1, True)
        return maxima

    def smooth_image(self):
        self.images = []
        last_sigma = 0
        for sigma in self.parameters['sigma_array']:
            if last_sigma == 0 or sigma <= last_sigma:
                self.images.append(cv.GaussianBlur(self.original, Util.kernel_size(sigma), sigma))
            else:
                sigma_to_add = math.sqrt(pow(sigma, 2) - pow(last_sigma, 2))
                self.images.append(cv.GaussianBlur(self.images[-1], Util.kernel_size(sigma_to_add), sigma_to_add))
            last_sigma = sigma

    def compute_measures(self, skymap):
        """"""
        measures = {}

        for key in self.parameters['active_measures']:
            measures[key] = PointDataList(self.parameters[key]['dist_thresh'])

        self.smooth_image()
        for step in range(0, len(self.images)):
            sigma = self.parameters['sigma_array'][step]
            maxima = self.sorted_maxima(self.images[step], sigma)
            for key in self.parameters['active_measures']:
                getattr(self, self.parameters[key]['method'])(maxima, step, measures[key])

        for key in measures.keys():
            # print(measures[key].to_string())
            Util.write_list_to_json_file(measures[key].to_list(), 'Measures/' + skymap + '_' + key + '.json')

        return measures


def main():
    sf = SourceFinder("conf.json")
    coords = sf.compute_coords()
    print(coords)
    plt.show()


# main()
