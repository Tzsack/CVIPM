import cv2 as cv
import numpy as np
from util import Util
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt


skymaps = ["../Skymaps/01-skymap_220-0_46-0.fits",
           "../Skymaps/02-skymap_221-0_46-0.fits",
           "../Skymaps/03-skymap_221-0_46-0_4-0.fits",
           "../Skymaps/04-skymap_221-0_47-0.fits",
           "../Skymaps/05-skymap_221-0_47-0_0-0.fits",
           "../Skymaps/06-skymap_221-0_47-0_0-5.fits",
           "../Skymaps/07-skymap_221-0_47-0_0-8.fits",
           "../Skymaps/08-skymap_221-0_47-0_1-0.fits",
           "../Skymaps/09-skymap_221-5_45-5.fits",
           "../Skymaps/10-skymap_bgonly.fits"
           ]


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    @staticmethod
    def display_data(img, maxima):
        cv.imshow("original", img)

        for i in range(0, len(img)):
            for j in range(0, len(img[i])):
                img[i][j] = 0

        for i in range(0, len(maxima)):
            img[maxima[i][0][0]][maxima[i][0][1]] = 1

        cv.imshow("maxima", img)

        # print(maxima)
        maxima.sort(key=lambda x: x[1], reverse=True)
        print(maxima)

        plot_values = []
        for i in range(1, len(maxima)):
            plot_values.append(maxima[i][1] - maxima[i - 1][1])

        plt.subplot(2, 1, 1)
        plt.plot(plot_values)
        plt.subplot(2, 1, 2)
        plt.plot([x[1] for x in maxima])

    def source_pixels_routine(self, image):
        maximum = np.max(image, axis=None)
        thresholded = np.copy(image)
        # print(maximum)
        for i in range(0, len(image)):
            for j in range(0, len(image)):
                if image[i][j] < maximum * 0:
                    thresholded[i][j] = 0
                else:
                    thresholded[i][j] = image[i][j]

        maxima = Util.local_maxima(thresholded)
        self.display_data(thresholded, maxima)

    def source_pixels(self, sigma_min, sigma_step, steps):
        """Search the source pixels position in the image"""
        smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma_min), sigma_min)
        self.source_pixels_routine(smoothed)

        for i in range(0, steps):
            smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_step), sigma_step)
            plt.figure()
            self.source_pixels_routine(smoothed)

        return 0, 0

    @staticmethod
    def cut(arr, threshold):
        result = []
        for element in arr:
            if element >= threshold:
                result.append(element)
        return result

    def source_pixels2(self, sigma_min=1.0, sigma_step=0.5, steps=10):
        candidates = []
        smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma_min), sigma_min)
        maxima = Util.local_maxima(smoothed)
        print(Util.sort_list(maxima, 1, True))
        maxima_val = Util.project_list(maxima, 1)
        print(maxima_val[0] - maxima_val[1])

        for i in range(0, steps):
            smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_step), sigma_step)
            maxima = Util.local_maxima(smoothed)
            print(Util.sort_list(maxima, 1, True))
            maxima_val = Util.project_list(maxima, 1)
            print(maxima_val[0] - maxima_val[1])

        return 0, 0


def main(index):
    sf = SourceFinder(skymaps[index])
    print(skymaps[index])
    # cv.imshow("a", sf.matrix)
    """
    maxima = Util.local_maxima(sf.matrix)
    Util.sort_list(maxima, 1, True)
    maxima_values = Util.project_list(maxima, 1)
    # print(maxima_values)

    smoothed = cv.GaussianBlur(sf.matrix, (3, 3), 1.5)
    # cv.imshow("b", smoothed)
    s_maxima = Util.local_maxima(smoothed)
    s_maxima = Util.sort_list(s_maxima, 1, True)

    s_maxima_values = Util.project_list(s_maxima, 1)
    s_maxima_pos = Util.project_list(s_maxima, 0)
    print(s_maxima_values)
    print(s_maxima_pos)
    x, y = Util.from_pix_to_wcs((s_maxima_pos[0][1], s_maxima_pos[0][0]), sf.wcs)
    print(x, y)
    """
    pixels_coord = sf.source_pixels2(1, 2, 5)
    print(pixels_coord)
    lon, lat = Util.from_pix_to_wcs(pixels_coord, sf.wcs)
    print(lon, lat)
    plt.show()
    cv.waitKey(0)
    cv.destroyAllWindows()


# test main
main(2)
