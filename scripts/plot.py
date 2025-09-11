import os
import sys
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import utils


def slater_orbital(r, n, zeta):
    N = (2 * zeta) ** (n + 0.5) / np.sqrt(scipy.special.factorial(2 * n))
    return N * r ** (n - 1) * np.exp(-zeta * r)


def clementi_roetti_orbital(r, coeffs, ns, zetas):
    R = np.zeros_like(r)
    for c, n, z in zip(coeffs, ns, zetas):
        R += c * slater_orbital(r, n, z)
    return R


if __name__ == "__main__":
    z = 10 if len(sys.argv) < 2 else int(sys.argv[1])
    if os.path.exists(f"../densities/{z}.csv"):
        df = pd.read_csv(f"../densities/{z}.csv")
    elif os.path.exists(f"densities/{z}.csv"):
        df = pd.read_csv(f"densities/{z}.csv")
    else:
        print("Could not find densities folder.")
        exit(1)
    r = df["r"]

    fig, axs = plt.subplots(2, 1, figsize=(4, 6))

    # electron density
    axs[0].sharex(axs[1])
    axs[0].plot(r, 4 * np.pi * df["density"] * r ** 2, color="C0")
    axs[0].set(ylabel=r"Radial distribution function $ P(r) = 4 \pi r^2 \rho(r) $", title=f"{utils.chemical_symbols[z]} (Z={z})")
    axs[0].tick_params('x', labelbottom=False)
    axs[0].grid()

    # orbitals from DFT
    line_styles = ["-", "--", "-.", ":"]
    for i, c in enumerate(filter(lambda s: s[0].isnumeric(), df.columns)):
        axs[1].plot(r, df[c], label=f"orbital {c}", color="C3", linestyle=line_styles[i % len(line_styles)])

    # reference plots for neon
    if z == 10:
        neon_exponents = [9.48486, 15.56590, 1.96184, 2.86423, 4.82530, 7.79242]
        neon_n = [1, 1, 2, 2, 2, 2]
        neon_1s = clementi_roetti_orbital(r, [0.93717, 0.04899, 0.00058, -0.00064, 0.00551, 0.01999], neon_n, neon_exponents)
        neon_2s = clementi_roetti_orbital(r, [-0.23093, -0.00635, 0.18620, 0.66899, 0.30910, -0.13871], neon_n, neon_exponents)
        neon_2p = clementi_roetti_orbital(r, [0.21799, 0.53338, 0.32933, 0.01872], [2, 2, 2, 2], [1.45208, 2.38168, 4.48489, 9.13464])
        axs[1].plot(r, neon_1s * r, label="reference 1s", color="0.2", linestyle=line_styles[0])
        axs[1].plot(r, neon_2s * r, label="reference 2s", color="0.2", linestyle=line_styles[1])
        axs[1].plot(r, neon_2p * r, label="reference 2p", color="0.2", linestyle=line_styles[2])


    axs[1].set(xlim=[-0.1, 5.1], ylabel=r"Radial orbitals $ u(r) = r \cdot R(r) $", xlabel=r"r [Bohr]")
    axs[1].grid()
    axs[1].legend()
    plt.tight_layout()
    fig.align_ylabels()
    plt.show()
