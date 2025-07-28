//! Contains integration functions based on Simpson's 1/3 rule similar to scipy's
//! implementations. Based on Wikipedia's "Simpson's rule" entry and
//! Cartwright, Kenneth V. (September 2017) http://msme.us/2017-2-1.pdf

/// Uses the composite simpson rule to integrate y(x)
/// samples of y are equally spaced by dx
pub fn simpson(y: &Vec<f64>, dx: f64) -> f64 {
    // at least 3 entries required
    debug_assert!(y.len() >= 3);
    // number of intervals
    let n = y.len() - 1;
    // apply simpson rule on two intervals at a time
    let mut sum = 0.0;
    for i in 0..n / 2 {
        sum += y[2 * i] + 4.0 * y[2 * i + 1] + y[2 * i + 2];
    }
    // if there is an uneven number of intervals (even number of samples) apply correction
    // for last interval according to
    if n % 2 != 0 {
        sum += 1.25 * y[n] + 2.0 * y[n - 1] - 0.25 * y[n - 2];
    }
    dx / 3.0 * sum
}

/// integrates y(x) cumulatively using Simpson's rule adding one sample at a time
/// with an initial value of 0.0
pub fn cumulative_simpson(y: &Vec<f64>, dx: f64) -> Vec<f64> {
    debug_assert!(y.len() >= 3);
    let n = y.len();
    let mut res = vec![0.0; n];
    // equations are based on scipy's _cumulatively_sum_simpson_integrals()
    for i in 1..(n + 1) / 2 {
        // h1 in scipy's source code
        res[2 * i - 1] = res[2 * i - 2]
            + dx / 3.0 * (1.25 * y[2 * i - 2] + 2.0 * y[2 * i - 1] - 0.25 * y[2 * i]);
        // h2 in scipy's source dode
        res[2 * i] = res[2 * i - 1]
            + dx / 3.0 * (1.25 * y[2 * i] + 2.0 * y[2 * i - 1] - 0.25 * y[2 * i - 2]);
    }
    if n % 2 == 0 {
        // last interval uses formula for h2
        res[n - 1] = res[n - 2] + dx / 3.0 * (1.25 * y[n - 1] + 2.0 * y[n - 2] - 0.25 * y[n - 3]);
    }
    res
}

/// same as simpson() but for unevenly sized intervals
pub fn simpson_unequal(y: &Vec<f64>, x: &Vec<f64>) -> f64 {
    debug_assert!(y.len() >= 3);
    // x should be in ascending order
    debug_assert_eq!(y.len(), x.len());
    debug_assert!(x.first().unwrap() < x.last().unwrap());
    let n = y.len() - 1;
    // vector containing the sizes of the intervals
    let h: Vec<f64> = (0..n).map(|i| x[i + 1] - x[i]).collect();
    let mut sum = 0.0;
    for i in 0..n / 2 {
        sum += (h[2 * i] + h[2 * i + 1]) / 6.0
            * ((2.0 - h[2 * i + 1] / h[2 * i]) * y[2 * i]
                + (h[2 * i] + h[2 * i + 1]).powi(2) / (h[2 * i] * h[2 * i + 1]) * y[2 * i + 1]
                + (2.0 - h[2 * i] / h[2 * i + 1]) * y[2 * i + 2]);
    }
    if n % 2 != 0 {
        sum += (2.0 * h[n - 1].powi(2) + 3.0 * h[n - 1] * h[n - 2]) / (6.0 * (h[n - 2] + h[n - 1]))
            * y[n]
            + (h[n - 1].powi(2) + 3.0 * h[n - 1] * h[n - 2]) / (6.0 * h[n - 2]) * y[n - 1]
            - h[n - 1].powi(3) / (6.0 * h[n - 2] * (h[n - 2] + h[n - 1])) * y[n - 2];
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts;

    fn almost_equal_vec_f64(a: &Vec<f64>, b: &Vec<f64>, epsilon: f64) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            if a[i] != b[i] {
                assert!((a[i] - b[i]).abs() < epsilon);
            }
        }
    }

    #[test]
    fn test_simpson_pow() {
        let y: Vec<f64> = (0..10).map(|i: i32| i.pow(3) as f64).collect();
        let res = simpson(&y, 1.0);
        assert_eq!(res, 1640.5);
    }

    #[test]
    fn test_simpson_sin() {
        let dx = 3.0 * consts::PI / 99.0;
        let y = (0..100).map(|i| (i as f64 * dx).sin()).collect();
        let res = simpson(&y, dx);
        assert!((2.0000043205470206 - res).abs() < 1e-8)
    }

    #[test]
    fn test_cumulative_linear() {
        let y: Vec<f64> = (0..10).map(f64::from).collect();
        let result = cumulative_simpson(&y, 1.0);
        let reference = [0.0, 0.5, 2.0, 4.5, 8.0, 12.5, 18.0, 24.5, 32.0, 40.5];
        assert_eq!(result, reference);
    }

    #[test]
    fn test_cumulative_sin() {
        let dx = 2.0 * consts::PI / 29.0;
        let y: Vec<f64> = (0..30).map(|i| (i as f64 * dx).sin()).collect();
        let result = cumulative_simpson(&y, dx);
        let reference = vec![
            0.00000000e+00,
            2.34694019e-02,
            9.24257182e-02,
            2.03982539e-01,
            3.52618056e-01,
            5.31641111e-01,
            7.32480679e-01,
            9.45877707e-01,
            1.16179630e+00,
            1.37012104e+00,
            1.56120628e+00,
            1.72595008e+00,
            1.85688004e+00,
            1.94759014e+00,
            1.99416251e+00,
            1.99407123e+00,
            1.94767715e+00,
            1.85680136e+00,
            1.72601674e+00,
            1.56115476e+00,
            1.37015502e+00,
            1.16178144e+00,
            9.45872736e-01,
            7.32505242e-01,
            5.31598104e-01,
            3.52677497e-01,
            2.03909445e-01,
            9.25090485e-02,
            2.33797321e-02,
            -8.96697799e-05,
        ];
        almost_equal_vec_f64(&result, &reference, 1e-8);
    }
}
