import numpy as np
import scipy
import jax
import jax.numpy as jnp


def slater_x(x, r, rho):
    e_x = - 3 * np.pi * (3 / np.pi) ** (1 / 3) * scipy.integrate.simpson(rho ** (4 / 3) * r ** 3, x)
    v_x = -(3 / np.pi) ** (1 / 3) * rho ** (1 / 3)
    return e_x, v_x


def epsilon_vwn(rs):
    # epsilon_c is taken from libxc (maple code for VWN)
    # epsilon_x is taken from the original paper (and divided by 2)
    a = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    # helper functions
    fx = rs + b * jnp.sqrt(rs) + c
    q = jnp.sqrt(4 * c - b ** 2)
    f1 = 2 * b / q
    f2 = b * x0 / (x0 ** 2 + b * x0 + c)
    f3 = 2 * (2 * x0 + b) / q

    epsilon_c = a * (jnp.log(rs / fx) + (f1 - f2 * f3) * jnp.arctan(q / (2 * jnp.sqrt(rs) + b))
                     - f2 * jnp.log((jnp.sqrt(rs) - x0) ** 2 / fx))
    epsilon_x = - (3 / 4) * (3 / (2 * jnp.pi)) ** (2 / 3) * (1 / rs)
    return epsilon_x + epsilon_c


grad_epsilon_vwn = jax.vmap(jax.value_and_grad(epsilon_vwn))


def vwn_xc(x, r, rho):
    rho = jnp.where(rho < 1e-12, 1e-12, rho)
    rs = (3 / (4 * jnp.pi * rho)) ** (1 / 3)

    epsilon, grad_epsilon = grad_epsilon_vwn(rs)
    e_xc = 4 * np.pi * scipy.integrate.simpson(rho * epsilon * r ** 3, x)
    v_xc = epsilon - (rs / 3) * grad_epsilon
    return np.array(e_xc), np.array(v_xc)
