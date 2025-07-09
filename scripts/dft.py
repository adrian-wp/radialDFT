import numpy as np
import scipy
import matplotlib.pyplot as plt

from utils import generate_occupations
from functionals import slater_x, vwn_xc


class SingleAtomDFT:
    def __init__(self, z, occupations, r_max, n_grid, xc_functional=None):
        self.z = z
        self.occupations = occupations
        # simpson rule works better with uneven number of samples
        if n_grid % 2 == 0:
            n_grid = n_grid + 1
        self.r = np.linspace(1e-7, r_max, n_grid)
        self.h = r_max / (n_grid - 1)
        self.rho = np.zeros_like(self.r)
        self.e_total = None
        self.xc_functional = xc_functional
        if xc_functional is None:
            self.xc_functional = slater_x

    def get_v_external(self):
        e_ext = -4 * np.pi * self.z * scipy.integrate.simpson(self.rho * self.r, self.r)
        v_ext = -self.z / self.r
        return e_ext, v_ext

    def get_v_hartree(self):
        u_diagonal = np.full(self.r.shape[0], -2 / self.h ** 2)
        u_off = np.full(self.r.shape[0] - 1, 1 / self.h ** 2)
        rhs = -4 * np.pi * self.rho * self.r

        u_diagonal[0] = 1
        u_off[0] = 0
        rhs[0] = 0

        u_diagonal[-1] = 1
        u_off[-1] = 0
        rhs[-1] = sum(map(sum, self.occupations))

        # solve banded matrix
        ab = np.zeros(shape=(3, self.r.shape[0]))
        ab[0, 1:] = u_off
        ab[1] = u_diagonal
        ab[2, :-1] = u_off
        u = scipy.linalg.solve_banded((1, 1), ab, rhs)
        v_hartree = u / self.r
        v_hartree[0] = v_hartree[1]

        e_hartree = 2 * np.pi * scipy.integrate.simpson(v_hartree * self.rho * self.r ** 2, self.r)
        return e_hartree, v_hartree

    def solve_ks(self, l, v_eff):
        u_diagonal = 1 / self.h ** 2 + l * (l + 1) / (2 * self.r ** 2) + v_eff
        u_off = np.full(self.r.shape[0] - 1, - 1 / (2 * self.h ** 2))
        u_diagonal[0] = 1
        u_diagonal[-1] = 1
        u_off[0] = 0
        u_off[-1] = 0

        values, vectors = scipy.linalg.eigh_tridiagonal(u_diagonal, u_off)

        n_vectors = len(self.occupations[l])
        eigen_values = values[:n_vectors]
        eigen_vectors = np.empty(shape=(n_vectors, self.r.shape[0]))
        for n in range(n_vectors):
            psi = vectors[:, n]
            psi[0] = 0
            psi[-1] = 0
            norm = scipy.integrate.simpson(psi ** 2, self.r)
            eigen_vectors[n] = psi / np.sqrt(norm)
        return eigen_values, eigen_vectors

    def construct_rho(self, vectors, mixing=0.7):
        new_rho = np.zeros_like(self.rho)
        for l in range(len(self.occupations)):
            for n in range(len(self.occupations[l])):
                new_rho += self.occupations[l][n] * vectors[l][n] ** 2
        new_rho /= (4 * np.pi * self.r ** 2)
        # linear extrapolation
        new_rho[0] = 2 * new_rho[1] - new_rho[2]
        self.rho = mixing * self.rho + (1 - mixing) * new_rho

    def kinetic_energy(self, eigen_values, v_eff):
        eigen_sum = 0
        for l in range(len(self.occupations)):
            for n in range(len(self.occupations[l])):
                eigen_sum += self.occupations[l][n] * eigen_values[l][n]
        return eigen_sum - 4 * np.pi * scipy.integrate.simpson(v_eff * self.rho * self.r ** 2, self.r)

    def electron_count(self):
        return 4 * np.pi * scipy.integrate.simpson(self.rho * self.r ** 2, self.r)

    def run_scf(self, convergence=1e-5, max_iterations=100):
        v_eff = np.zeros_like(self.r)
        for i in range(max_iterations):
            print(f"Iteration {i}")
            eigen_values, eigen_vectors = list(), list()
            for l in range(len(self.occupations)):
                values_l, vectors_l = self.solve_ks(l, v_eff)
                eigen_values.append(values_l)
                eigen_vectors.append(vectors_l)
            if i == 0:
                self.construct_rho(eigen_vectors, mixing=0)
            else:
                self.construct_rho(eigen_vectors)

            e_ext, v_ext = self.get_v_external()
            e_hartree, v_hartree = self.get_v_hartree()
            e_xc, v_xc = self.xc_functional(self.r, self.rho)
            v_eff = v_ext + v_xc + v_hartree
            e_kin = self.kinetic_energy(eigen_values, v_eff)
            e_total = e_kin + e_ext + e_xc + e_hartree
            print(f"\tElectrons: {self.electron_count(): .5g}\n"
                  f"\tE_kinetic: {e_kin: .5g}\n"
                  f"\tE_hartree: {e_hartree: .5g}\n"
                  f"\tE_nucleus: {e_ext: .5g}\n"
                  f"\tE_xc     : {e_xc: .5g}\n"
                  f"\tE_total  : {e_total: .5g}")
            if self.e_total is not None:
                print(f"\tDelta    : {e_total - self.e_total: .5g}")

            if self.e_total is not None and np.abs(e_total - self.e_total) < convergence:
                print("Reached desired convergence.")
                return e_total, e_kin, e_hartree, e_ext, e_xc, i
            self.e_total = e_total
        print("Maximum number of iterations reached.")


if __name__ == "__main__":
    dft = SingleAtomDFT(8, generate_occupations(8), 30, 3001, vwn_xc)
    dft.run_scf()

    fig, ax = plt.subplots()
    ax.plot(dft.r, dft.rho)
    ax.set(xlabel=r"r", ylabel=r"rho", xlim=[-0.1, 5.1])
    ax.grid()
    plt.show()
