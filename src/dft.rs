//! Module containing the main code of the DFT implementation.
//! single_atom_dft() is the main function that should be called to perform
//! the SCF loop. All other functions are called during the SCF loop.

use std::f64::consts;
use crate::integrate;
use crate::linalg;

/// Used to pass the integer occupations to the DFT function.
#[derive(Clone)]
pub struct Occupations {
    // integer occupation number for each (l, n) e.g. [[2, 2], [4], [], []] for Oxygen
    pub f: [Vec<i32>; 4],
    // for how many l is f not empty; e.g. 2 for Oxygen
    pub max_l: usize,
}

/// Contains the grid points of the log grid with n points, where the points in
/// x are equally spaced by h_x and the real space points are r=e^x.
pub struct LogGrid {
    pub n: usize,
    pub x: Vec<f64>,
    pub r: Vec<f64>,
    pub h_x: f64,
}

/// Stores the eigenvalues and eigenvectors for each l.
pub struct Orbitals {
    pub eigenvalues: [Vec<f64>; 4],
    pub eigenvectors: [Vec<f64>; 4],
}

/// All the components of the total energy for a single iteration.
pub struct Energy {
    pub external: f64,
    pub hartree: f64,
    pub xc: f64,
    pub kinetic: f64,
    pub total: f64,
}

/// The return value after the SCF loop converged. Contains the energies,
/// the density, the orbitals etc.
pub struct DFTResult {
    pub energy: Energy,
    pub density: Vec<f64>,
    pub orbitals: Orbitals,
    pub grid: LogGrid,
    pub occupations: Occupations,
    pub iterations: usize,
    pub success: bool,
}

/// Alias for the function pointer used to pass the XC functionals.
pub type XCFunctional = fn(&LogGrid, &Vec<f64>, &mut Vec<f64>) -> f64;

/// Adds the external potential to v_eff and returns the potential energy.
fn add_v_external(z: i32, grid: &LogGrid, rho: &Vec<f64>, v_eff: &mut Vec<f64>) -> f64 {
    // external potential
    for i in 0..grid.n {
        v_eff[i] += -z as f64 / grid.r[i];
    }
    // external energy
    let integrand = (0..grid.n).map(|i| rho[i] * grid.r[i].powi(2)).collect();
    -4.0 * consts::PI * z as f64 * integrate::simpson(&integrand, grid.h_x)
}

/// Adds the hartree potential to v_eff and returns the hartree energy. Uses a
/// Greens function approach to solve the poisson equation.
fn add_v_hartree(grid: &LogGrid, rho: &Vec<f64>, v_eff: &mut Vec<f64>) -> f64 {
    // each point is determined from two integrals
    // using the cumulative integration the integrals are determined for all grid points at once
    let mut inner_integral = (0..grid.n).map(|i| rho[i] * grid.r[i].powi(3)).collect();
    let mut outer_integral = (0..grid.n).map(|i| rho[i] * grid.r[i].powi(2)).rev().collect();
    inner_integral = integrate::cumulative_simpson(&inner_integral, grid.h_x);
    outer_integral = integrate::cumulative_simpson(&outer_integral, grid.h_x);
    // create the hartree potential and the integrand for the hartree energy
    let mut integrand = vec![0.0; grid.n];
    for i in 0..grid.n {
        // the outer integral was reversed to cumulatively integrate it from the
        // back, it needs to be accessed in reverse to fix the order again
        let v_hartree_i = 4.0 * consts::PI *
            (inner_integral[i] / grid.r[i] + outer_integral[grid.n - i - 1]);
        v_eff[i] += v_hartree_i;
        integrand[i] = v_hartree_i * rho[i] * grid.r[i].powi(3);
    }
    2.0 * consts::PI * integrate::simpson(&integrand, grid.h_x)
}

/// Solves the Kohn-Sham equations for the specified l, grid and v_eff.
/// n_shells determines how many eigenvalues and eigenvectors are calculated.
/// They are returned in ascending order: first the eigenvalues, then all
/// eigenvectors concatenated in a single vector.
fn solve_kohn_sham(grid: &LogGrid, l: usize, v_eff: &Vec<f64>, n_shells: usize)
                   -> (Vec<f64>, Vec<f64>) {
    // second-order finite differences lead to a symmetric tridiagonal matrix
    // the equations were adjusted for the log grid r = e^x
    // initialize stencil matrix and right side diagonal matrix B
    // omit last row for Dirichlet boundary at r_max
    let diag = (0..grid.n - 1).map(|i|
        1.0 / grid.h_x.powi(2) + 0.125 + 0.5 * (l * (l + 1)) as f64 + grid.r[i].powi(2) * v_eff[i]
    ).collect();
    let off = vec![-0.5 / grid.h_x.powi(2); grid.n - 2];
    let mut a = linalg::TridiagonalMatrix { n: grid.n - 1, diag, off };
    let mut b: Vec<f64> = (0..grid.n - 1).map(|i| grid.r[i].powi(2)).collect();
    // Neumann BC at r_min
    a.diag[0] += l as f64 + 0.5 / grid.h_x;
    a.diag[0] /= 2.0;
    b[0] /= 2.0;
    // the LAPACK function does not work with a generalized problem A * x = epsilon * B * x
    // instead solve C * y = epsilon * y, transform C = B^(-1/2) * A * B^(-1/2) in-place
    for i in 0..grid.n - 2 {
        a.diag[i] /= b[i];
        a.off[i] /= (b[i] * b[i + 1]).sqrt();
    }
    a.diag[grid.n - 2] /= b[grid.n - 2];

    // solve eigenvalue problem
    let (mut epsilon, y) = linalg::eigh_tridiagonal(a, n_shells as i32);

    // sort eigenvalues
    let mut sort_index = (0..n_shells).collect::<Vec<_>>();
    sort_index.sort_by(|&a, &b| epsilon[a].partial_cmp(&epsilon[b]).unwrap());
    epsilon = sort_index.iter().map(|&i| epsilon[i]).collect();

    // sort eigenvectors and append 0.0 boundary point
    let mut u = vec![0.0; grid.n * n_shells];
    for dst in 0..n_shells {
        let src = sort_index[dst];
        for i in 0..grid.n - 1 {
            // destination vectors are 1 longer than source vectors
            let src_index = src * (grid.n - 1) + i;
            let dst_index = dst * grid.n + i;
            // transform back from y to u during sort
            u[dst_index] = y[src_index] / b[i].sqrt() * grid.r[i].sqrt();
        }
    }

    // normalize u
    for j in 0..n_shells {
        let integrand = (0..grid.n).map(|i| u[j * grid.n + i].powi(2) * grid.r[i]).collect();
        let norm = integrate::simpson(&integrand, grid.h_x).sqrt();
        for i in 0..grid.n {
            u[j * grid.n + i] /= norm;
        }
    }
    (epsilon, u)
}

/// Constructs and returns the density from the given eigenvectors.
fn construct_rho(grid: &LogGrid, u: &[Vec<f64>; 4], occupations: &Occupations) -> Vec<f64> {
    let mut new_rho = vec![0.0; grid.n];
    for l in 0..occupations.max_l {
        for n in 0..occupations.f[l].len() {
            for i in 0..grid.n {
                new_rho[i] += occupations.f[l][n] as f64 * u[l][n * grid.n + i].powi(2) /
                    (4.0 * consts::PI * grid.r[i].powi(2));
            }
        }
    }
    new_rho
}

/// Performs pulay mixing by consuming the old and new densities and returning
/// the resulting density. The histories of the residuals and densities need to
/// be passed and will be modified to contain the new residual and density.
fn pulay_mixing(grid: &LogGrid, densities: &mut Vec<Vec<f64>>, residuals: &mut Vec<Vec<f64>>,
                old_rho: Vec<f64>, mut new_rho: Vec<f64>, steps: usize, i: usize) -> Vec<f64> {
    // calculate residual and save residual and density in history
    let residual = (0..grid.n).map(|i| new_rho[i] - old_rho[i]).collect();
    residuals[i % steps] = residual;
    densities[i % steps] = new_rho.clone();
    if i == 0 {
        // in first iteration just leave density as is
    } else if i < steps - 1 {
        // use linear mixing until history is filled
        for j in 0..grid.n {
            new_rho[j] = old_rho[j] * 0.8 + new_rho[j] * 0.2;
        }
    } else {
        // pulay mixing using history
        // create linear system Ax=b
        let mut a = vec![-1.0; (steps + 1).pow(2)];
        let mut b = vec![0.0; steps + 1];
        a[(steps + 1).pow(2) - 1] = 0.0;
        b[steps] = -1.0;
        for x in 0..steps {
            // only fill lower triangle of A
            for y in x..steps {
                // vector dot product of residuals, weighted by r^3
                let a_i = x * (steps + 1) + y;
                a[a_i] = (0..grid.n).fold(0.0, |acc, j|
                    acc + grid.r[j].powi(3) * residuals[x][j] * residuals[y][j]
                );
            }
        }
        // solve linear system
        let ls = linalg::LinearSystem { n: steps + 1, a, b };
        let coefficients = linalg::solve_symmetric(ls, false);
        // linear superposition of densities using coefficients
        new_rho.fill(0.0);
        for j in 0..grid.n {
            for k in 0..steps {
                new_rho[j] += coefficients[k] * densities[k][j];
            }
        }
    }
    new_rho
}

/// Returns the kinetic energy for the given density, v_eff and eigenvalues.
fn get_kinetic_energy(grid: &LogGrid, rho: &Vec<f64>, v_eff: &Vec<f64>, eigenvalues: &[Vec<f64>; 4],
                      occupations: &Occupations) -> f64 {
    let mut eigenvalue_sum = 0.0;
    for l in 0..occupations.max_l {
        for n in 0..occupations.f[l].len() {
            eigenvalue_sum += occupations.f[l][n] as f64 * eigenvalues[l][n];
        }
    }
    let integrand = (0..grid.n).map(|i| v_eff[i] * rho[i] * grid.r[i].powi(3) ).collect();
    eigenvalue_sum - 4.0 * consts::PI * integrate::simpson(&integrand, grid.h_x)
}

/// Returns the electron count by integrating over the given density.
fn get_electron_count(grid: &LogGrid, rho: &Vec<f64>) -> f64 {
    let integrand = (0..grid.n).map(|i| rho[i] * grid.r[i].powi(3)).collect();
    4.0 * consts::PI * integrate::simpson(&integrand, grid.h_x)
}

/// Performs DFT for the specified parameters by starting an SCF loop until the
/// convergence criteria are met and returns a struct containing the final
/// density, energies and orbitals.
pub fn single_atom_dft(z: i32, r_min: f64, r_max: f64, n_grid: usize, occupations: Occupations,
                       mixing_steps: usize, max_iter: usize, min_iter: usize, e_tol: f64,
                       rho_tol: f64, xc_functional: XCFunctional, verbose: i32) -> DFTResult {
    // initialize grid
    // uniform grid in x from ln(r_min) to ln(r_max) with spacing h_x
    // log grid in r with r = e^x
    let x: Vec<f64> = (0..n_grid).map(|i|
        r_min.ln() + i as f64 / (n_grid - 1) as f64 * (r_max.ln() - r_min.ln())
    ).collect();
    let r = x.iter().map(|&x| consts::E.powf(x)).collect();
    let h_x = (x[x.len() - 1] - x[0]) / (n_grid - 1) as f64;
    let grid = LogGrid { n: n_grid, x, r, h_x };

    // density and residual history
    let mut density_history = vec![vec![0.0; grid.n]; mixing_steps];
    let mut residual_history = vec![vec![0.0; grid.n]; mixing_steps];
    // initial density and potential
    let mut rho = vec![0.0; n_grid];
    let mut v_eff = vec![0.0; n_grid];

    // to save orbitals and energies each iteration
    let mut orbitals = Orbitals {
        eigenvalues: [vec![], vec![], vec![], vec![]],
        eigenvectors: [vec![], vec![], vec![], vec![]],
    };
    let mut energy = Energy { external: 0.0, hartree: 0.0, xc: 0.0, kinetic: 0.0, total: 0.0 };

    // main SCF loop
    for i in 0..max_iter {
        // solve Kohn-Sham equations and construct new density
        for l in 0..occupations.max_l {
            (orbitals.eigenvalues[l], orbitals.eigenvectors[l]) =
                solve_kohn_sham(&grid, l, &v_eff, occupations.f[l].len());
        }
        let new_rho = construct_rho(&grid, &orbitals.eigenvectors, &occupations);

        // calculate difference to previous density, then apply mixing
        let rho_diff = (0..grid.n).fold(0.0, |diff, j|
            diff + (new_rho[j] - rho[j]).powi(2)
        ).sqrt();
        rho = pulay_mixing(&grid, &mut density_history, &mut residual_history, rho, new_rho,
                           mixing_steps, i);

        // get energies and new effective potential
        v_eff.fill(0.0);
        energy.external = add_v_external(z, &grid, &rho, &mut v_eff);
        energy.hartree = add_v_hartree(&grid, &rho, &mut v_eff);
        energy.xc = xc_functional(&grid, &rho, &mut v_eff);
        energy.kinetic = get_kinetic_energy(&grid, &rho, &v_eff, &orbitals.eigenvalues,
                                            &occupations);

        // calculate energy difference to previous iteration
        let new_e_total = energy.external + energy.kinetic + energy.hartree + energy.xc;
        let e_total_diff = (new_e_total - energy.total).abs();
        energy.total = new_e_total;

        // print
        if verbose >= 2 {
            let electron_count = get_electron_count(&grid, &rho);
            println!("Iteration {i}:\n\
                     \tElectrons: {:12.5}\n\
                     \tε        : {:?}\n\
                     \tE_xc     : {:12.5}\n\
                     \tE_nucleus: {:12.5}\n\
                     \tE_coulomb: {:12.5}\n\
                     \tE_kinetic: {:12.5}\n\
                     \tE_total  : {:12.5}",
                     electron_count, orbitals.eigenvalues, energy.xc, energy.external,
                     energy.hartree, energy.kinetic, energy.total);
        }
        // check convergence criteria
        if i > min_iter && e_total_diff < e_tol && rho_diff < rho_tol {
            if verbose >= 2 {
                println!("Reached desired convergence. (rho = {:.5e}, E = {:.5e})",
                         rho_diff, e_total_diff);
            } else if verbose == 1 {
                println!("Z = {:3}: {:4} {:12.5} {:12.5} {:12.5} {:12.5} {:12.5}    {:?}",
                         z, i, energy.total, energy.kinetic, energy.hartree, energy.external,
                         energy.xc, orbitals.eigenvalues);
            }
            return DFTResult { energy, density: rho, orbitals, grid, occupations, iterations: i,
                success: true };
        }
    }
    if verbose >= 2 { println!("Maximum number of iterations reached.") };
    DFTResult { energy, density: rho, orbitals, grid, occupations, iterations: max_iter,
        success: false }
}
