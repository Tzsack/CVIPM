import cv2 as cv
import numpy as np
from util import Util
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt


class SourceFinder:

    def __init__(self, file_name):
        """Constructor"""
        self.file_name = file_name
        self.matrix = None

    def load_fits_file(self):
        """Save the fits image in the matrix attribute"""
        file = fits.open(self.file_name)
        self.matrix = np.array(file[0].data, np.float64)
        file.close()

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
        maxima.sort(key=lambda x: x[1], reverse=False)
        # print(maxima)

        plot_values = []
        for i in range(1, len(maxima)):
            plot_values.append(maxima[i][1] - maxima[i - 1][1])

        plt.subplot(2, 1, 1)
        plt.plot(plot_values)
        plt.subplot(2, 1, 2)
        plt.plot([x[1] for x in maxima])

    def source_pixels_routine(self, image):
        maximum = np.max(image, axis=None)
        thresholded = image
        print(maximum)
        for i in range(0, len(image)):
            for j in range(0, len(image)):
                if image[i][j] < maximum * 0.5:
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

    def source_world(self, source_pixels):
        """Return world coordinates for the given pixel position source_pixels"""
        w = WCS(self.file_name)
        return w.all_pix2world(source_pixels[0], source_pixels[1], 0)

# Test main


fits_file = "../Skymaps/skymap_221-0_47-0_0-8.fits"
sf = SourceFinder(fits_file)
sf.load_fits_file()
pixels_coord = sf.source_pixels(1, 2, 5)
# lon, lat = sf.source_world(pixels_coord)
# print(lon, lat)
plt.show()
cv.waitKey(0)
cv.destroyAllWindows()
