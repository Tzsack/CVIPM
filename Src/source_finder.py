import cv2 as cv
from util import Util
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os


class SourceFinder:

    def __init__(self, fits_file_name):
        """Constructor"""
        self.matrix = Util.from_fits_to_mat(fits_file_name)
        self.wcs = WCS(fits_file_name)

    def source_pixels_routine(self, image, step, maxima_0):
        maxima = Util.sort_list(Util.local_maxima(image), 1, True)
        for maximum in Util.project_list(maxima, 0):
            Util.neigh8(maximum)
        # print()
        # print(maxima)
        # print("pos", Util.project_list(maxima, 0))
        # print("val", Util.project_list(maxima, 1))

        threshold = 0

        if maxima[0][1] - maxima[1][1] > threshold:
            x, y = Util.from_pix_to_wcs((maxima[0][0][1], maxima[0][0][0]), self.wcs)
            print("Step:", step, "-", maxima[0][0], x, y)
        else:
            print("Step:", step, "- No source found.")

    def source_pixels(self, sigma_min=1.0, sigma_step=0.5, steps=10):
        """Search the source pixels position in the image"""
        maxima_0 = Util.sort_list(Util.local_maxima(self.matrix), 1, True)
        # print(maxima_0)
        # print("pos", Util.project_list(maxima_0, 0))
        # print("val", Util.project_list(maxima_0, 1))
        smoothed = cv.GaussianBlur(self.matrix, Util.kernel_size(sigma_min), sigma_min)
        self.source_pixels_routine(smoothed, 0, maxima_0)

        for i in range(1, steps):
            smoothed = cv.GaussianBlur(smoothed, Util.kernel_size(sigma_step), sigma_step)
            self.source_pixels_routine(smoothed, i, maxima_0)

        return 0, 0


def main(index):
    os.chdir("../Skymaps")
    skymaps = os.listdir(".")
    skymaps.sort()
    index -= 1
    sf = SourceFinder(skymaps[index])
    print(skymaps[index])
    #pixels_coord = sf.source_pixels(1, 2, 5)
    # print(pixels_coord)
    #lon, lat = Util.from_pix_to_wcs(pixels_coord, sf.wcs)
    #print(lon, lat)
    #plt.show()
    cv.waitKey(0)
    cv.destroyAllWindows()


def main2(indexes, sigma):
    os.chdir("../Skymaps")
    skymaps = os.listdir(".")
    skymaps.sort()

    index = 1
    sf = SourceFinder(skymaps[index - 1])
    mat = cv.GaussianBlur(sf.matrix, Util.kernel_size(sigma), sigma)
    maxima = Util.local_maxima(mat)
    # cv.imshow(skymaps[index - 1] + str(index), Util.flip_mat(mat))
    Util.sort_list(maxima, 1, True)
    # print(skymaps[index - 1], maxima)
    # plt.plot(Util.project_list(maxima, 1), label=''+str(index))
    # plt.legend()
    # cv.waitKey(0)
    # cv.destroyAllWindows()


main(1)

"""
# test main
for i in range(0, 10):
    main2([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 0.5*i + 0.01)
    #plt.figure()
#plt.show()
"""
