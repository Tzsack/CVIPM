import numpy as np
from astropy.io import fits
import math
import json


class Util:

    @staticmethod
    def from_fits_to_mat(fits_file_name):
        """Return numpy matrix from fits file fits_file_name"""
        file = fits.open(fits_file_name)
        matrix = np.array(file[0].data, np.float64)
        file.close()
        return matrix

    @staticmethod
    def flip_mat(mat):
        return np.flip(mat, 0)

    @staticmethod
    def from_pix_to_wcs(pixel_coord, wcs):
        """Return world coordinates for the given pixel coordinates pixel_coord"""
        return wcs.all_pix2world(pixel_coord[1], pixel_coord[0], 0)

    @staticmethod
    def kernel_size(sigma):
        """Determine kernel size according to sigma"""
        size = int(sigma * 6)
        if size >= 5:
            if size % 2 == 0:
                size += 1
            return size, size
        else:
            return 5, 5

    @staticmethod
    def distance_eu(a, b):
        return math.sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2))

    @staticmethod
    def distance4(a, b):
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    @staticmethod
    def neigh4(p):
        return [
            (p[0] - 1, p[1]),
            (p[0], p[1] - 1),
            (p[0] + 1, p[1]),
            (p[0], p[1] + 1)
        ]

    @staticmethod
    def distance8(a, b):
        return max(abs(a[0] - b[0]), abs(a[1] - b[1]))

    @staticmethod
    def neigh8(p):
        return [
            (p[0] - 1, p[1] - 1),
            (p[0] - 1, p[1]),
            (p[0] - 1, p[1] + 1),
            (p[0], p[1] - 1),
            (p[0], p[1] + 1),
            (p[0] + 1, p[1] - 1),
            (p[0] + 1, p[1]),
            (p[0] + 1, p[1] + 1)
        ]

    @staticmethod
    def neigh8_radius(pixel, radius):
        starters = [pixel]
        derived = []
        for i in range(1, radius + 1):
            starters.append((pixel[0] - i, pixel[1]))
            starters.append((pixel[0] + i, pixel[1]))
        for i in range(1, radius + 1):
            for starter in starters:
                derived.append((starter[0], starter[1] - i))
                derived.append((starter[0], starter[1] + i))
        starters.remove(pixel)
        neighbourhood = starters + derived
        # neighbourhood.sort(key=lambda x: (x[0], x[1]))
        return neighbourhood

    @staticmethod
    def local_maxima(img):
        maxima = []
        maximum = ((0, 0), img[0][0])
        for i in range(0, len(img)):
            for j in range(0, len(img[i])):
                is_max = True
                if img[i][j] != 0:
                    if img[i][j] > maximum[1]:
                        maximum = ((i, j), img[i][j])
                    for p in Util.neigh8((i, j)):
                        if 0 <= p[0] < len(img) and 0 <= p[1] < len(img[i]) and img[p[0]][p[1]] > img[i][j]:
                            is_max = False
                    if is_max:
                        maxima.append(((i, j), img[i][j]))
        if maxima.count(maximum) == 0:
            maxima.append(maximum)
        return maxima

    @staticmethod
    def project_list(l, pos):
        return [x[pos] for x in l]

    @staticmethod
    def sort_list(l, pos, reverse=False):
        l.sort(key=lambda x: x[pos], reverse=reverse)
        return l

    @staticmethod
    def read_list_from_json_file(file_name):
        """Return list of dictionaries from json file file_name"""
        with open(file_name) as f:
            return json.load(f)

    @staticmethod
    def json_beautifier(json_data, sort_keys):
        """Return json_data in a beautified output"""
        return json.dumps(json.loads(json_data), indent=4, sort_keys=sort_keys)

    @staticmethod
    def from_list_to_json(a_list):
        """Return json of a_list"""
        return json.dumps(a_list)

    @staticmethod
    def write_list_to_json_file(a_list, file_name, sort_keys=False):
        """Save a_list in json beautified format in the file file_name"""
        beauty_output = Util.json_beautifier(json_data=Util.from_list_to_json(a_list), sort_keys=sort_keys)
        text_file = open(file_name, "w")
        text_file.write("%s" % beauty_output)
        text_file.close()

    @staticmethod
    def blob_test(mat):
        res = np.copy(mat)
        for i in range(1, len(mat) - 1):
            for j in range(1, len(mat[i]) - 1):
                u = 0
                d = 1
                if mat[i+1][j] >= mat[i][j]:
                    u += 1
                    d *= -1
                if mat[i][j+1] >= mat[i][j]:
                    u += 1
                if mat[i-1][j] >= mat[i][j]:
                    u += 1
                    d *= -1
                if mat[i][j-1] >= mat[i][j]:
                    u += 1

                if u == 4:
                    res[i][j] = 0  # minimum
                if u == 3:
                    res[i][j] = 0  # valley
                if u == 1:
                    res[i][j] = 1  # crest
                if u == 2:
                    if d > 0:
                        res[i][j] = 0  # saddle
                    else:
                        res[i][j] = 1  # slope
                if u == 0:
                    res[i][j] = 0.5  # maximum
        return res

    @staticmethod
    def minmax_normalize(mat):
        return (mat - np.amin(mat)) * 1/(np.amax(mat) - np.amin(mat))

    @staticmethod
    def binarize(mat, t):
        res = np.copy(mat)
        for i in range(0, len(res)):
            for j in range(0, len(res[i])):
                if res[i][j] < t:
                    res[i][j] = 0
                else:
                    res[i][j] = 1
        return res
