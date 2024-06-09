""" Imports """"
import numpy as np

def unbiased_pol2cart(r, beta, sd_r, sd_beta):
    """
    Date: 14 Dec 2020
    Author: Dr Sanjeev Arulampalam
    Modified: 9 Jun 2024 by Randall McCutcheon
    
    This is the unbiased polar to Cartesian conversion based on Longbin Mo's
    paper IEEE Transactions on Aerospace and Electronic Systems,
    Vol 34. No. 3, July 1998, p1023-1027

    RM note: Title of the paper is: Unbiased converted measurements for tracking
    
    Inputs:
        r - range
        beta - Bearing
        sd_r - Standard Deviation of Range
        sd_beta - Standard Deviation of Bearing

    Outputs:
        z -
        R -

    9 Jun - modifications include the ability to input arrays
    """
  
    # Equation 12
    lambda_beta = np.exp(-np.power(sd_beta,2) / 2)
    lambda_beta_prime = np.power(lambda_beta, 4)

    inv_lambda_beta = 1 / lambda_beta
    inv_lambda_beta_prime = 1 / lambda_beta_prime

    # Equation 5
    xm_u = inv_lambda_beta * r * np.cos(bearing)
    ym_u = inv_lambda_beta * r * np.sin(bearing)

    # Equation 7a
    R11 = (np.power(lambda_beta, -2) - 2) * np.power(r, 2) * np.power(np.cos(beta), 2) + 0.5 * (np.power(r, 2) + np.power(sd_r, 2)) * (1 + lambda_beta_prime * np.cos(2 * beta))

    # Equation 7b
    R22 = (np.power(lambda_beta, -2) - 2) * np.power(r, 2) * np.power(np.sin(beta), 2) + 0.5 * (np.power(r, 2) + np.power(sd_r, 2)) * (1 - lambda_beta_prime * np.cos(2 * beta))

    # Equation 7c
    R12 = (np.power(lambda_beta, -2) - 2) * np.power(r, 2) * np.cos(beta) * np.sin(beta) + 0.5 * (np.power(r, 2) + np.power(sd_r, 2)) * lambda_beta_prime * np.sin(2 * beta)

    z = np.transpose([xm_u, ym_u])
    R = np.array([[R11, R12], [R12, R22]])

    return z, R
