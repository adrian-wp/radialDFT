use crate::integrate;
use crate::dft;
use std::f64::consts;

#[allow(dead_code)]
pub fn slater_exchange(grid: &dft::LogGrid, rho: &Vec<f64>, v_eff: &mut Vec<f64>) -> f64 {
    for i in 0..grid.n {
        v_eff[i] += -(3.0 / consts::PI).powf(1.0 / 3.0) * rho[i].powf(1.0 / 3.0);
    }
    let integrand = (0..grid.n).map(|i| {
        rho[i].powf(4.0 / 3.0) * grid.r[i].powi(3)
    }).collect();
    let e_x = -3.0 * consts::PI * (3.0 / consts::PI).powf(1.0 / 3.0) * integrate::simpson_unequal(&integrand, &grid.x);
    e_x
}

fn epsilon_vwn(rs: f64) -> (f64, f64) {
    // constants
    let a: f64 = 0.0310907;
    let b: f64 = 3.72744;
    let c: f64 = 12.9352;
    let x0: f64 = -0.10498;
    // helper functions
    let q = (4.0 * c - b.powi(2)).sqrt();
    let f1 = 2.0 * b / q;
    let f2 = b * x0 / (x0.powi(2) + b * x0 + c);
    let f3 = 2.0 * (2.0 * x0 + b) / q;
    let fx = rs + b * rs.sqrt() + c;

    // calculate epsilon_xc
    let epsilon_x = -(3.0 / 4.0) * (3.0 / (2.0 * consts::PI)).powf(2.0 / 3.0) / rs;
    let epsilon_c = a * ((rs / fx).ln() + (f1 - f2 * f3) * (q / (2.0 * rs.sqrt() + b)).atan()
        - f2 * ((rs.sqrt() - x0).powi(2) / fx).ln());
    let epsilon_xc = epsilon_x + epsilon_c;

    // calculate grad_epsilon_xc
    let fx_prime = 1.0 + b / (2.0 * rs.sqrt());
    let epsilon_x_prime = (3.0 / 4.0) * (3.0 / (2.0 * consts::PI)).powf(2.0 / 3.0) / rs.powi(2);
    let epsilon_c_prime = a * ((1.0 / rs) + (f2 - 1.0) * (fx_prime / fx) - f2 / (rs - x0 * rs.sqrt())
        - (f1 - f2 * f3) * q / (rs.sqrt() * (q.powi(2) + (2.0 * rs.sqrt() + b).powi(2))));
    let grad_epsilon_xc = epsilon_x_prime + epsilon_c_prime;

    (epsilon_xc, grad_epsilon_xc)
}

pub fn vwn_xc(grid: &dft::LogGrid, rho: &Vec<f64>, v_eff: &mut Vec<f64>) -> f64 {
    let mut integrand = vec![0.0; grid.n];
    for i in 0..grid.n {
        let rho_i = rho[i].max(1e-12);
        let rs = (3.0 / (4.0 * consts::PI * rho_i)).powf(1.0 / 3.0);
        let epsilon;
        let grad_epsilon;
        (epsilon, grad_epsilon) = epsilon_vwn(rs);
        v_eff[i] += epsilon - (rs / 3.0) * grad_epsilon;
        integrand[i] = rho_i * epsilon * grid.r[i].powi(3);
    }
    4.0 * consts::PI * integrate::simpson(&integrand, grid.h_x)
}
