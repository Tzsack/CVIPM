from gaussian_interpolation import GInterpolation as gi
import os
from source_finder import SourceFinder
import cv2 as cv
from util import Util
import numpy as np


def test_gi_on_synth(mu, sigma, prefactor):
    curdir = os.getcwd()
    os.chdir("../../Skymaps/")
    skymaps = os.listdir(".")
    print(skymaps)
    skymaps.sort()
    index = 1
    sf = SourceFinder(skymaps[index - 1])
    os.chdir(curdir)

    for i in range(0, len(sf.matrix)):
        for j in range(0, len(sf.matrix[0])):
            sf.matrix[i][j] = prefactor * gi.gaussian((i, j), mu, sigma)

    cv.imshow("img", sf.matrix * 1 / (np.amax(sf.matrix)))

    maxima = Util.sort_list(Util.local_maxima(sf.matrix), 1, True)
    maximum = maxima[0]

    estimated_mu = gi.estimate_mu(sf.matrix, maximum[0], sigma)
    estimated_prefactor = gi.estimate_prefactor(sf.matrix, estimated_mu, sigma)
    print("Apex =", mu, prefactor * gi.gaussian(mu, mu, sigma))
    print("Maximum =", maximum[0], maximum[1])
    print("Estimated Apex =", estimated_mu, estimated_prefactor * gi.gaussian(estimated_mu, estimated_mu, sigma))
    print("Estimated Prefactor =", estimated_prefactor)


def test_gi_on_image():
    curdir = os.getcwd()
    os.chdir("../../Skymaps/")
    skymaps = os.listdir(".")
    skymaps.sort()
    index = 1
    sf = SourceFinder(skymaps[index - 1])
    os.chdir(curdir)

    sigma = 5.681

    sf.matrix = cv.GaussianBlur(sf.matrix, Util.kernel_size(sigma), sigma)
    cv.imshow("img", sf.matrix * 1/(np.amax(sf.matrix)))

    maxima = Util.sort_list(Util.local_maxima(sf.matrix), 1, True)
    maximum = maxima[0]

    estimated_mu = gi.estimate_mu(sf.matrix, maximum[0], sigma)
    estimated_prefactor = gi.estimate_prefactor(sf.matrix, estimated_mu, sigma)

    print("Maximum =", maximum[0], maximum[1])
    print("Estimated Apex =", estimated_mu, estimated_prefactor * gi.gaussian(estimated_mu, estimated_mu, sigma))
    print("Estimated Prefactor =", estimated_prefactor)


test_gi_on_image()
test_gi_on_synth((201.5, -3.45), 3.52, 5.52)

cv.waitKey(0)
cv.destroyAllWindows()
