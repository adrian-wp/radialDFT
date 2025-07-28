import numpy as np
import scipy


def slater_x(x, r, rho):
    rho = np.maximum(rho, 1e-12)
    e_x = - 3 * np.pi * (3 / np.pi) ** (1 / 3) * scipy.integrate.simpson(rho ** (4 / 3) * r ** 3, x)
    v_x = -(3 / np.pi) ** (1 / 3) * rho ** (1 / 3)
    return e_x, v_x


def epsilon_vwn(rs):
    a = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    # helper functions
    fx = rs + b * np.sqrt(rs) + c
    q = np.sqrt(4 * c - b ** 2)
    f1 = 2 * b / q
    f2 = b * x0 / (x0 ** 2 + b * x0 + c)
    f3 = 2 * (2 * x0 + b) / q

    # calculate epsilon_xc
    epsilon_x = - (3 / 4) * (3 / (2 * np.pi)) ** (2 / 3) / rs
    epsilon_c = a * (np.log(rs / fx) + (f1 - f2 * f3) * np.arctan(q / (2 * np.sqrt(rs) + b))
                     - f2 * np.log((np.sqrt(rs) - x0) ** 2 / fx))
    epsilon_xc = epsilon_x + epsilon_c

    # calculate grad_epsilon_xc
    fx_prime = 1 + b / (2 * np.sqrt(rs))
    epsilon_x_prime = (3 / 4) * (3 / (2 * np.pi)) ** (2 / 3) / rs ** 2
    epsilon_c_prime = a * ((1 / rs) + (f2 - 1) * (fx_prime / fx) - f2 / (rs - x0 * np.sqrt(rs))
                           - (f1 - f2 * f3) * q / (np.sqrt(rs) * (q ** 2 + (2 * np.sqrt(rs) + b) ** 2)))
    grad_epsilon_xc = epsilon_x_prime + epsilon_c_prime
    return epsilon_xc, grad_epsilon_xc


def vwn_xc(x, r, rho):
    rho = np.maximum(rho, 1e-12)
    rs = (3 / (4 * np.pi * rho)) ** (1 / 3)

    epsilon, grad_epsilon = epsilon_vwn(rs)

    e_xc = 4 * np.pi * scipy.integrate.simpson(rho * epsilon * r ** 3, x)
    v_xc = epsilon - (rs / 3) * grad_epsilon
    return e_xc, v_xc
