import numpy as np
import math
import cv2 as cv


class GInterpolation:

    @staticmethod
    def mu_c_est(c1, z1, c2, z2, sigma):
        return (pow(c1, 2) - pow(c2, 2) + 2*pow(sigma, 2)*math.log(z1/z2))/(2*(c1 - c2))

    @staticmethod
    def gaussian(p, mu, sigma):
        return 1/(2*pow(sigma, 2)*math.pi)*math.exp((-pow(p[0] - mu[0], 2) - pow(p[1] - mu[1], 2))/(2*pow(sigma, 2)))

    @staticmethod
    def estimate_mu(mat, maximum, sigma, radius):
        if radius <= 0:
            return maximum
        x_a = maximum[0]
        y_a = maximum[1]
        x_b = maximum[0]
        y_b = maximum[1]
        mu_x = 0
        mu_y = 0
        for i in range(0, radius):
            if x_a + 1 < len(mat):
                x_a += 1
            if x_b > 0:
                x_b -= 1
            if y_a + 1 < len(mat[0]):
                y_a += 1
            if y_b > 0:
                y_b -= 1

            dmu_x = GInterpolation.mu_c_est(x_a, mat[x_a][maximum[1]], x_b, mat[x_b][maximum[1]], sigma)
            dmu_y = GInterpolation.mu_c_est(y_a, mat[maximum[0]][y_a], y_b, mat[maximum[0]][y_b], sigma)
            # print(dmu_x, dmu_y)
            mu_x += dmu_x
            mu_y += dmu_y
        mu_x /= radius
        mu_y /= radius
        return mu_x, mu_y
