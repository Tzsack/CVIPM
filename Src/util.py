import cv2 as cv
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt


class Util:

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
    def neigh4(p):
        return [
            (p[0] - 1, p[1]),
            (p[0], p[1] - 1),
            (p[0] + 1, p[1]),
            (p[0], p[1] + 1)
        ]

    @staticmethod
    def local_maxima(img):
        maxima = []
        for i in range(1, len(img) - 1):
            for j in range(1, len(img[i]) - 1):
                is_max = True
                for p in Util.neigh8((i, j)):
                    if img[p[0]][p[1]] >= img[i][j]:
                        is_max = False
                if is_max:
                    maxima.append(((i, j), img[i][j]))
        return maxima


file = fits.open("../Skymaps/skymap_221-0_47-0_0-8.fits")
im = np.array(file[0].data, np.float64)
file.close()
cv.imshow("1", im)

sigma = 2
im = cv.GaussianBlur(im, Util.kernel_size(sigma), sigma)
cv.imshow("2", im)
maxima_ = Util.local_maxima(im)
for ii in range(0, len(im)):
    for jj in range(0, len(im[ii])):
        im[ii][jj] = 0

for ii in range(0, len(maxima_)):
    im[maxima_[ii][0][0]][maxima_[ii][0][1]] = 1

print(maxima_)

maxima_.sort(key=lambda x: x[1], reverse=True)

print(maxima_)

cv.imshow("3", im)

plt.plot([1,2,3,4])
plt.ylabel('some numbers')
plt.show()


cv.waitKey(0)
cv.destroyAllWindows()
