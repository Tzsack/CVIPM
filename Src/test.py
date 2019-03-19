import cv2 as cv
import numpy as np
from astropy.io import fits

# file = fits.open("../Skymaps/skymap_221_47_1.fits")
# img_input = np.array(file[0].data, np.float64)

img_input = cv.imread("../../cat.jpg")

kernel = 1 * np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]], np.float64)

out_1 = cv.GaussianBlur(img_input, (3, 3), 2)

out_2 = cv.filter2D(out_1, -1, kernel)

cv.imshow("b", img_input)
cv.imshow("3", out_1)
cv.imshow("laplace", out_2)

cv.waitKey(0)
cv.destroyAllWindows()
