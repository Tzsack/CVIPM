from gaussian_interpolation import GInterpolation as gi
import os
from source_finder import SourceFinder
import cv2 as cv
from util import Util
import numpy as np


def test_gi_on_synth(mu, sigma, prefactor):
    curdir = os.getcwd()
    os.chdir("../../Tests/gi")
    skymaps = os.listdir(".")
    print(skymaps)
    skymaps.sort()
    index = 1
    mat = Util.from_fits_to_mat(skymaps[index - 1])
    os.chdir(curdir)

    for i in range(0, len(mat)):
        for j in range(0, len(mat[0])):
            mat[i][j] = prefactor * gi.gaussian((i, j), mu, sigma)

    cv.imshow("img", mat * 1 / (np.amax(mat)))

    maxima = Util.sort_list(Util.local_maxima(mat), 1, True)
    maximum = maxima[0]

    estimated_mu = gi.estimate_mu(mat, maximum[0], sigma)
    estimated_prefactor = gi.estimate_prefactor(mat, estimated_mu, sigma)
    print("Apex =", mu, prefactor * gi.gaussian(mu, mu, sigma))
    print("Maximum =", maximum[0], maximum[1])
    print("Estimated Apex =", estimated_mu, estimated_prefactor * gi.gaussian(estimated_mu, estimated_mu, sigma))
    print("Estimated Prefactor =", estimated_prefactor)


def test_gi_on_image():
    curdir = os.getcwd()
    os.chdir("../../Tests/gi")
    skymaps = os.listdir(".")
    skymaps.sort()
    index = 1
    mat = Util.from_fits_to_mat(skymaps[index - 1])
    os.chdir(curdir)

    sigma = 4.0

    mat = cv.GaussianBlur(mat, Util.kernel_size(sigma), sigma)
    cv.imshow("img", mat * 1/(np.amax(mat)))

    maxima = Util.sort_list(Util.local_maxima(mat), 1, True)
    maximum = maxima[0]

    estimated_mu = gi.estimate_mu(mat, maximum[0], sigma)
    estimated_prefactor = gi.estimate_prefactor(mat, estimated_mu, sigma)

    print("Maximum =", maximum[0], maximum[1])
    print("Estimated Apex =", estimated_mu, estimated_prefactor * gi.gaussian(estimated_mu, estimated_mu, sigma))
    print("Estimated Prefactor =", estimated_prefactor)

    synth = np.copy(mat)

    synth = Util.blob_test(mat)

    # for i in range(0, len(synth)):
    #     for j in range(0, len(synth[0])):
    #         synth[i][j] = 0
    #
    # n = 0
    # for maximum in maxima:
    #     estimated_mu = gi.estimate_mu(mat, maximum[0], sigma)
    #     estimated_prefactor = gi.estimate_prefactor(mat, estimated_mu, sigma)
    #     n += 1
    #     print(n, len(maxima))
    #     for i in range(0, len(synth)):
    #         for j in range(0, len(synth[0])):
    #             synth[i][j] += gi.gaussian((i, j), estimated_mu, sigma) * estimated_prefactor
    #
    # cv.imshow("synth", synth * 1 / np.amax(synth))
    # for i in range(0, len(synth)):
    #     for j in range(0, len(synth[0])):
    #         synth[i][j] = 0
    # n = 0
    # for maximum in maxima:
    #     es = gi.optimize_maximum_z(mat, maximum[0])
    #     estimated_mu = (es[0], es[1])
    #     estimated_sigma = es[2]
    #     estimated_prefactor = es[3]
    #     n += 1
    #     print(n, len(maxima))
    #     for i in range(0, len(synth)):
    #         for j in range(0, len(synth[0])):
    #             synth[i][j] += gi.gaussian((i, j), estimated_mu, estimated_sigma) * estimated_prefactor

    cv.imshow("synth_z", synth * 1/np.amax(synth))


def test():
    curdir = os.getcwd()
    os.chdir("../../Tests/gi")
    skymaps = os.listdir(".")
    skymaps.sort()
    index = 1
    mat = Util.from_fits_to_mat(skymaps[index - 1])
    os.chdir(curdir)

    sigma = 4.0

    mat = cv.GaussianBlur(mat, Util.kernel_size(sigma), sigma)
    cv.imshow("img", (mat - np.amin(mat)) * 1/(np.amax(mat) - np.amin(mat)))

    synth = Util.blob_test(mat)
    cv.imshow("synth_z", (synth - np.amin(synth)) * 1/(np.amax(synth) - np.amin(synth)))

    kernel = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]])
    # kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    synth = cv.filter2D(mat, -1, kernel=kernel)
    synth = Util.minmax_normalize(synth)
    cv.imshow("laplacian 0", synth)
    mask = Util.binarize(synth, 0.6)
    cv.imshow("binarized laplacian 0", mask)

    for i in range(0, len(mat)):
        for j in range(0, len(mat[i])):
            synth[i][j] = mask[i][j] * mat[i][j]

    synth = Util.minmax_normalize(synth)

    cv.imshow("masked", synth)

test()
# test_gi_on_synth((201.5, -3.45), 3.52, 5.52)

cv.waitKey(0)
cv.destroyAllWindows()
