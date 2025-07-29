import time

import numpy as np
import scipy

from utils import generate_occupations, occupation_string
from functionals import slater_x, vwn_xc


class SingleAtomDFT:
    def __init__(self, z, r_min, r_max, n_grid, xc_functional=None, occupations=None):
        self.z = z
        # occupation string
        if occupations is None:
            self.occupations = generate_occupations(z)
        else:
            self.occupations = occupation_string(occupations)
        # simpson rule works better with uneven number of samples
        if n_grid % 2 == 0:
            n_grid = n_grid + 1
        self.x = np.linspace(np.log(r_min), np.log(r_max), n_grid)
        self.dx = (self.x[-1] - self.x[0]) / (n_grid - 1)
        self.r = np.exp(self.x)
        self.rho = np.zeros_like(self.r)
        self.xc_functional = xc_functional
        if xc_functional is None:
            self.xc_functional = slater_x
        # used to access energies after completing SCF loop
        self.e_kin = None
        self.e_ext = None
        self.e_xc = None
        self.e_hartree = None
        self.e_total = None
        # histories used for pulay mixing
        self.densities = list()
        self.residuals = list()

    def get_v_external(self):
        e_ext = -4 * np.pi * self.z * scipy.integrate.simpson(self.rho * self.r ** 2, self.x)
        v_ext = -self.z / self.r
        return e_ext, v_ext

    def get_v_hartree(self):
        integral_inner = scipy.integrate.cumulative_simpson(self.rho * self.r ** 3, dx=self.dx, initial=0) / self.r
        integral_outer = scipy.integrate.cumulative_simpson((self.rho * self.r ** 2)[::-1], dx=self.dx, initial=0)[::-1]
        v_hartree = 4 * np.pi * (integral_inner + integral_outer)
        e_hartree = 2 * np.pi * scipy.integrate.simpson(v_hartree * self.rho * self.r ** 3, self.x)
        return e_hartree, v_hartree

    def solve_ks(self, l, v_eff):
        # solve generalized hermitian eigenvalue problem A * chi = epsilon * B * chi
        # with diagonal matrix B
        a_diagonal = 1 / self.dx ** 2 + 1 / 8 + l * (l + 1) / 2 + self.r ** 2 * v_eff
        a_off = np.full(self.r.shape[0] - 1, - 1 / (2 * self.dx ** 2))
        b = self.r ** 2

        # Neumann BC: chi'(x_min) = (l+1/2) * chi(x_min)
        a_diagonal[0] += (l + 0.5) / self.dx
        a_diagonal[0] /= 2
        b[0] /= 2

        # C = B^(-1/2) * A * B^(-1/2)
        c_diagonal = a_diagonal / b
        c_off = a_off / np.sqrt(b[:-1] * b[1:])

        # solve eigenvalue problem C * y = epsilon * y for first n eigenvalues
        # remove last values for Dirichlet boundary conditions
        occupation = len(self.occupations[l])

        epsilon, y = scipy.linalg.eigh_tridiagonal(c_diagonal[:-1], c_off[:-1], select="i",
                                                   select_range=(0, occupation - 1), lapack_driver="stebz")

        # transform back to chi and append zeros for boundary
        chi = np.zeros(shape=(self.r.shape[0], occupation))
        chi[:-1] = y / np.sqrt(b[:-1, np.newaxis])
        # transform back to u
        u = chi * np.sqrt(self.r)[:, np.newaxis]

        # normalize u
        for n in range(occupation):
            norm = scipy.integrate.simpson(u[:, n] ** 2 * self.r, self.x)
            u[:, n] /= np.sqrt(norm)
        return epsilon, u

    def construct_rho(self, eigen_vectors):
        new_rho = np.zeros_like(self.rho)
        for l in range(len(self.occupations)):
            for n in range(len(self.occupations[l])):
                new_rho += self.occupations[l][n] * eigen_vectors[l][:, n] ** 2
        new_rho /= (4 * np.pi * self.r ** 2)
        return new_rho

    def pulay_mixing(self, new_rho, steps=5):
        self.densities.append(new_rho)
        self.residuals.append(new_rho - self.rho)
        if len(self.residuals) >= steps:
            last_residuals = np.array(self.residuals[-steps:])
            weights = self.r ** 3
            weighted_residuals = last_residuals * weights[np.newaxis, :]
            overlap = np.full(shape=(steps + 1, steps + 1), fill_value=-1.0)
            overlap[:steps, :steps] = weighted_residuals @ last_residuals.T
            overlap[-1, -1] = 0
            rhs = np.zeros(steps + 1)
            rhs[-1] = -1
            coefficients = scipy.linalg.solve(overlap, rhs)[:-1]
            last_densities = np.array(self.densities[-steps:])
            return np.tensordot(last_densities, coefficients, (0, 0))
        elif len(self.residuals) > 1:
            return self.rho * 0.8 + new_rho * 0.2
        else:
            return new_rho

    def kinetic_energy(self, eigen_values, v_eff):
        eigen_sum = 0
        for l in range(len(self.occupations)):
            for n in range(len(self.occupations[l])):
                eigen_sum += self.occupations[l][n] * eigen_values[l][n]
        return eigen_sum - 4 * np.pi * scipy.integrate.simpson(v_eff * self.rho * self.r ** 3, self.x)

    def electron_count(self):
        return 4 * np.pi * scipy.integrate.simpson(self.rho * self.r ** 3, self.x)

    def run_scf(self, e_tol=1e-5, rho_tol=1e-3, max_iter=100, min_iter=5):
        v_eff = np.zeros_like(self.r)
        for i in range(max_iter):
            print(f"Iteration {i}")

            # get eigen values and vectors for each l
            eigen_values, eigen_vectors = list(), list()
            for l in range(len(self.occupations)):
                values_l, vectors_l = self.solve_ks(l, v_eff)
                eigen_values.append(values_l)
                eigen_vectors.append(vectors_l)

            # new rho and mixing
            new_rho = self.construct_rho(eigen_vectors)
            rho_diff = np.linalg.norm(new_rho - self.rho)
            self.rho = self.pulay_mixing(new_rho)

            # get energies and v_eff with new rho
            self.e_ext, v_ext = self.get_v_external()
            self.e_hartree, v_hartree = self.get_v_hartree()
            self.e_xc, v_xc = self.xc_functional(self.x, self.r, self.rho)
            v_eff = v_ext + v_xc + v_hartree
            self.e_kin = self.kinetic_energy(eigen_values, v_eff)
            new_e_total = self.e_kin + self.e_ext + self.e_xc + self.e_hartree
            if self.e_total is not None:
                e_total_diff = np.abs(self.e_total - new_e_total)
            else:
                e_total_diff = None
            self.e_total = new_e_total

            print(f"\tElectrons: {self.electron_count(): .5g}\n"
                  f"\tε        : {[a.tolist() for a in eigen_values]}\n"
                  f"\tE_xc     : {self.e_xc: .5g}\n"
                  f"\tE_nucleus: {self.e_ext: .5g}\n"
                  f"\tE_coulomb: {self.e_hartree: .5g}\n"
                  f"\tE_kinetic: {self.e_kin: .5g}\n"
                  f"\tE_total  : {self.e_total: .5g}")
            if e_total_diff is not None:
                print(f"\tDelta E  : {e_total_diff: .5g}\n"
                      f"\tDelta rho: {rho_diff: .5g}")

            # check convergence criteria
            if i >= min_iter and rho_diff < rho_tol and e_total_diff < e_tol:
                print(f"Reached desired convergence. (rho = {rho_diff:.5g}, E = {e_total_diff:.5g})")
                return True, i
        print("Maximum number of iterations reached.")
        return False, max_iter


if __name__ == "__main__":
    start = time.time()
    dft = SingleAtomDFT(1, 1e-4, 100, 3001, xc_functional=vwn_xc)
    dft.run_scf()
    print(f"{time.time() - start:.3f}s elapsed.")
