import numpy as np
from astropy.io import fits


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
        return wcs.all_pix2world(pixel_coord[0], pixel_coord[1], 0)

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
