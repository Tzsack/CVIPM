import cv2 as cv
import numpy as np
from astropy.io import fits

file = fits.open("../Skymaps/skymap_221_47_1.fits")
img_input = np.array(file[0].data, np.float64)

# kernel = 1 * np.array([[1/16, 1/8, 1/16], [1/8, 1/4, 1/8], [1/16, 1/8, 1/16]], np.float64)
# out = cv.filter2D(img_input, -1, kernel)

out_1 = cv.GaussianBlur(img_input, (3, 3), 2)
out_2 = cv.GaussianBlur(img_input, (5, 5), 2)
out_3 = cv.GaussianBlur(img_input, (7, 7), 2)
out_4 = cv.GaussianBlur(img_input, (9, 9), 2)
out_5 = cv.GaussianBlur(img_input, (11, 11), 2)
out_6 = cv.GaussianBlur(img_input, (13, 13), 2)

cv.imshow("b", img_input)

cv.imshow("3", out_1)
cv.imshow("5", out_2)
cv.imshow("7", out_3)
cv.imshow("9", out_4)
cv.imshow("11", out_5)
cv.imshow("13", out_6)

cv.waitKey(0)
cv.destroyAllWindows()
