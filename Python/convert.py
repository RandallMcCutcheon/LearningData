""" Imports """"
import numpy as np

def unbiased_pol2cart(range, bearing, range_sd, bearing_sd):
    """
    Date: 14 Dec 2020
    Author: Dr Sanjeev Arulampalam
    Modified: 9 Jun 2024 by Randall McCutcheon
    
    This is the unbiased polar to Cartesian conversion based on Longbin Mo's
    paper IEEE Transactions on Aerospace and Electronic Systems,
    Vol 34. No. 3, July 1998, p1023-1027.

    Inputs:
        r - range
        beta - Bearing
        var_r - Standard Deviation of Range
        var_beta - Standard Deviation of Bearing

    Outputs:
        z -
        R -

    9 Jun - modifications include ability to input arrays
    """
  
    # Equation 12
    lambda_beta = np.exp(-bearing_sd / 2)
    lambda_beta2 = np.power(lambda_beta1, 4)

    inv_lambda_beta1 = 1 / lambda_beta1
    inv_lambda_beta2 = 1 / lambda_beta2

    xm_u = inv_lambda_beta1 * r * np.cos(beta)
    ym_u = inv_lambda_beta1 * r * np.sin(beta)

    R11 = (np.power(inv_lambda_beta1, 2) - 2) * np.power(r, 2) * np.power(
        np.cos(beta), 2
    ) + 0.5 * (np.power(r, 2) + var_r) * (1 + lambda_beta2 * np.cos(2 * beta))

    R22 = (np.power(inv_lambda_beta1, 2) - 2) * np.power(r, 2) * np.power(
        np.sin(beta), 2
    ) + 0.5 * (np.power(r, 2) + var_r) * (1 - lambda_beta2 * np.cos(2 * beta))

    R12 = (np.power(inv_lambda_beta1, 2) - 2) * np.power(r, 2) * np.cos(beta) * np.sin(
        beta
    ) + 0.5 * (np.power(r, 2) + var_r) * lambda_beta2 * np.sin(2 * beta)

    z = np.transpose([xm_u, ym_u])
    R = np.array([[R11, R12], [R12, R22]])

    return z, R
