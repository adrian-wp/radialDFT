//! Wrappers of functions of the LAPACK library that we need for solving
//! eigenvalue problems and systems of linear equations

#[allow(unused_imports)]
use lapack_src;

/// Used to pass a linear system of equations to the LAPACK wrapper.
pub struct LinearSystem {
    // size N of the system
    pub n: usize,
    // matrix NxN in column-major format
    pub a: Vec<f64>,
    // right-hand side b of the system, size N
    pub b: Vec<f64>,
}

/// Used to pass a symmetric tridiagonal matrix to the LAPACK wrapper.
pub struct TridiagonalMatrix {
    // size of matrix NxN
    pub n: usize,
    // diagonal of tridiagonal matrix with length N
    pub diag: Vec<f64>,
    // off diagonals of symmetric matrix with length N-1
    pub off: Vec<f64>,
}

/// Solves the linear system Ax=b, where the matrix A is symmetric. Consumes
/// the linear system `ls` and returns the solution vector. `upper` specifies
/// whether the upper or lower triangle of the will matrix should be accessed.
pub fn solve_symmetric(mut ls: LinearSystem, upper: bool) -> Vec<f64> {
    debug_assert_eq!(ls.a.len(), ls.n.pow(2));
    debug_assert_eq!(ls.b.len(), ls.n);
    // whether the upper or lower triangle of A should be accessed
    let uplo = if upper { b'U' } else { b'L' };
    // pivot indices that define permutation matrix P
    let mut ipiv = vec![0; ls.n];
    // workspace size and work array (not optimal but should be sufficient)
    let lwork = ls.n;
    let mut work = vec![0.0; lwork];
    // returns 0 if successful
    let mut info = -1;
    // call LAPACK
    unsafe {
        lapack::dsysv(
            uplo,
            ls.n as i32,
            1,
            ls.a.as_mut_slice(),
            ls.n as i32,
            ipiv.as_mut_slice(),
            ls.b.as_mut_slice(),
            ls.n as i32,
            work.as_mut_slice(),
            lwork as i32,
            &mut info,
        );
    }
    assert_eq!(info, 0, "LAPACKs DSYSV returned non-zero info code.");
    ls.b
}

/// Calculates the first `iu` eigenvalues and eigenvectors of the symmetric
/// tridiagonal matrix A by using LAPACKs DSTEBZ and DSTEIN. Consumes the
/// tridiagonal matrix `a` and returns two vectors, the first contains the
/// eigenvalues, the second all the eigenvectors concatenated in a single
/// vector. Eigenvalues are NOT ORDERED.
pub fn eigh_tridiagonal(mut a: TridiagonalMatrix, iu: i32) -> (Vec<f64>, Vec<f64>) {
    debug_assert_eq!(a.diag.len(), a.n);
    debug_assert_eq!(a.off.len(), a.n - 1);
    debug_assert!(a.n >= iu as usize);

    // output containers W and Z for the eigenvalues and eigenvectors
    let mut w = vec![0.0; a.n];
    let mut z = vec![0.0; a.n * iu as usize];

    // will return the number of eigenvalues found, should be iu in this case
    let mut m = 0;
    // returns the number of diagonal blocks in the matrix
    let mut nsplit = [0];
    // specifies to which block an eigenvalue belongs
    let mut iblock = vec![0; a.n];
    // contains the splitting points of the submatrices
    let mut isplit = vec![0; a.n];

    // workspace arrays of size 4*N and 3*N for DSTEBZ, and 5*N and N for DSTEIN
    let mut work = vec![0.0; 5 * a.n];
    let mut iwork = vec![0; 3 * a.n];
    // will return 0 on success
    let mut info = -1;

    // call LAPACK, request first iu eigenvalues
    unsafe {
        lapack::dstebz(
            b'I',
            b'B',
            a.n as i32,
            0.0,
            0.0,
            1,
            iu,
            0.0,
            a.diag.as_mut_slice(),
            a.off.as_mut_slice(),
            &mut m,
            nsplit.as_mut_slice(),
            w.as_mut_slice(),
            iblock.as_mut_slice(),
            isplit.as_mut_slice(),
            work.as_mut_slice(),
            iwork.as_mut_slice(),
            &mut info,
        );
    }
    assert_eq!(info, 0, "LAPACKs DSTEBZ returned non-zero info code.");

    // contains the indices for which the eigenvectors failed to converge
    let mut ifail = vec![0; iu as usize];

    // get eigenvectors from LAPACK
    unsafe {
        lapack::dstein(
            a.n as i32,
            a.diag.as_slice(),
            a.off.as_slice(),
            iu,
            w.as_slice(),
            iblock.as_slice(),
            isplit.as_slice(),
            z.as_mut_slice(),
            a.n as i32,
            work.as_mut_slice(),
            iwork.as_mut_slice(),
            ifail.as_mut_slice(),
            &mut info,
        );
    }
    assert_eq!(info, 0, "LAPACKs DSTEIN returned non-zero info code.");
    // we only requested iu eigenvalues
    w.truncate(iu as usize);
    (w, z)
}

/// Similar to `eigh_tridiagonal` but it uses LAPACKs DSTEMR which is based on
/// the MRRR algorithm. Eigenvalues are ordered in ascending order.
#[allow(unused)]
pub fn eigh_tridiagonal_mrrr(mut a: TridiagonalMatrix, iu: i32) -> (Vec<f64>, Vec<f64>) {
    // DSTEMR requires the off diagonal to have length N, even though the last value is not used
    if a.off.len() == a.n - 1 {
        a.off.push(0.0);
    }
    debug_assert_eq!(a.diag.len(), a.n);
    debug_assert_eq!(a.off.len(), a.n);
    debug_assert!(a.n >= iu as usize);

    // output containers W and Z for the eigenvalues and eigenvectors
    let mut w = vec![0.0; a.n];
    let mut z = vec![0.0; a.n * iu as usize];

    // will return the number of eigenvalues found, should be iu in this case
    let mut m = 0;
    // columns of Z, should be equal to the number of eigenvectors iu
    let nzc = [iu];
    // the indices indicating the nonzero elements in Z
    let mut isuppz = vec![0; 2 * iu as usize];
    // if true, try to use slower algorithm to get more accurate eigenvalues
    let mut tryrac = 0;
    // workspace sizes when calculating eigenvalues and vectors are LWORK >= 18*n, LIWORK >= 10*n
    let lwork = 18 * a.n;
    let liwork = 10 * a.n;
    let mut work = vec![0.0; lwork];
    let mut iwork = vec![0; liwork];
    // will return 0 on success
    let mut info = -1;

    // call LAPACK, request first iu eigenvalues
    unsafe {
        lapack::dstemr(
            b'V',
            b'I',
            a.n as i32,
            a.diag.as_mut_slice(),
            a.off.as_mut_slice(),
            0.0,
            0.0,
            1,
            iu,
            &mut m,
            w.as_mut_slice(),
            z.as_mut_slice(),
            a.n as i32,
            nzc.as_slice(),
            isuppz.as_mut_slice(),
            &mut tryrac,
            work.as_mut_slice(),
            lwork as i32,
            iwork.as_mut_slice(),
            liwork as i32,
            &mut info,
        );
    }
    // check for success
    assert_eq!(info, 0, "LAPACKs DSTEMR returned non-zero info code.");
    // we only want the first iu eigenvalues
    w.truncate(iu as usize);
    (w, z)
}

/// Some simple tests to check if the wrappers are working correctly. They are
/// not meant to test the LAPACK algorithms.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;

    #[test]
    fn test_solve_symmetric() {
        let ls_upper = LinearSystem {
            n: 5,
            a: vec![5.0, 0.0, 0.0, 0.0, 0.0,
                    11.0, 7.0, 0.0, 0.0, 0.0,
                    9.0, 17.0, 13.0, 0.0, 0.0,
                    4.0, 11.0, 17.0, 12.0, 0.0,
                    20.0, 10.0, 11.0, 14.0, 5.0],
            b: vec![2.0, 9.0, 20.0, 13.0, 11.0],
        };
        let ls_lower = LinearSystem {
            n: 5,
            a: vec![5.0, 11.0, 9.0, 4.0, 20.0,
                    0.0, 7.0, 17.0, 11.0, 10.0,
                    0.0, 0.0, 13.0, 17.0, 11.0,
                    0.0, 0.0, 0.0, 12.0, 14.0,
                    0.0, 0.0, 0.0, 0.0, 5.0],
            b: vec![2.0, 9.0, 20.0, 13.0, 11.0],
        };
        let reference = vec![-0.2740597, 1.61105723, 0.80288824, -0.23590904, -1.03168445];
        let solution_upper = solve_symmetric(ls_upper, true);
        utils::almost_equal_vec_f64(&reference, &solution_upper, 1e-6);
        let solution_lower = solve_symmetric(ls_lower, false);
        utils::almost_equal_vec_f64(&reference, &solution_lower, 1e-6);
    }

    #[test]
    fn test_eigh_tridiagonal() {
        let a = TridiagonalMatrix {
            n: 5,
            diag: vec![2.0, 6.0, 12.0, 12.0, 1.0],
            off: vec![20.0, 20.0, 3.0, 1.0],
        };
        let ref_vals = vec![-22.25694335, 0.90346881, 6.23861442, 12.7355039];
        let ref_vecs = vec![-0.57888644, 0.70210078, -0.41307466, 0.03621985, -0.00155738,
                            -0.02393171, 0.00131209, 0.02359735, -0.09603019, 0.99480999,
                            0.65912035, 0.13968785, -0.65745378, 0.33136256, 0.06325386,
                            -0.24733648, -0.13276409, 0.20262483, 0.93477104, 0.07965325];
        let (e_vals, e_vecs) = eigh_tridiagonal(a, 4);
        // in this case the values are ordered but this shouldn't be assumed
        utils::almost_equal_vec_f64(&e_vals, &ref_vals, 1e-6);
        utils::almost_equal_vec_f64(&e_vecs, &ref_vecs, 1e-6);
    }

    #[test]
    fn test_eigh_tridiagonal_mrrr() {
        let a = TridiagonalMatrix {
            n: 5,
            diag: vec![2.0, 6.0, 12.0, 12.0, 1.0],
            off: vec![20.0, 20.0, 3.0, 1.0],
        };
        let ref_vals = vec![-22.25694335, 0.90346881, 6.23861442, 12.7355039];
        let ref_vecs = vec![-0.57888644, 0.70210078, -0.41307466, 0.03621985, -0.00155738,
                            -0.02393171, 0.00131209, 0.02359735, -0.09603019, 0.99480999,
                            0.65912035, 0.13968785, -0.65745378, 0.33136256, 0.06325386,
                            -0.24733648, -0.13276409, 0.20262483, 0.93477104, 0.07965325];
        let (e_vals, e_vecs) = eigh_tridiagonal_mrrr(a, 4);
        utils::almost_equal_vec_f64(&e_vals, &ref_vals, 1e-6);
        utils::almost_equal_vec_f64(&e_vecs, &ref_vecs, 1e-6);
    }
}
