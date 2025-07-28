//! Wrappers of functions of the LAPACK library that we need for
//! solving eigenvalue problems and systems of linear equations

#[allow(unused_imports)]
use lapack_src;

/// used to pass a linear system of equations to the LAPACK wrapper
pub struct LinearSystem {
    // size N of the system
    pub n: usize,
    // matrix NxN in column-major format
    pub a: Vec<f64>,
    // right-hand side b of the system, size N
    pub b: Vec<f64>,
}

/// used to pass a tridiagonal matrix to the LAPACK wrapper
pub struct TridiagonalMatrix {
    // size of matrix NxN
    pub n: usize,
    // diagonal of tridiagonal matrix with length N
    pub diag: Vec<f64>,
    // off diagonals of symmetric matrix, LAPACK requires length N (the last element is not used)
    pub off: Vec<f64>,
}

/// Solves the linear system Ax=b, the result will be contained in b.
pub fn solve_symmetric(ls: &mut LinearSystem, upper: bool) {
    debug_assert_eq!(ls.a.len(), ls.n.pow(2));
    debug_assert_eq!(ls.b.len(), ls.n);
    // whether the upper or lower triangle of A should be accessed
    let uplo = if upper {b'U'} else {b'L'};
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
            &mut ls.a,
            ls.n as i32,
            &mut ipiv,
            &mut ls.b,
            ls.n as i32,
            &mut work,
            lwork as i32,
            &mut info,
        );
    }
    assert_eq!(info, 0);
    // no return value, the result is saved in ls.b
}

/// Calculates the first `iu` eigenvalues and eigenvectors of the tridiagonal matrix `a`
/// by using LAPACKs DSTEMR. Returns two vectors, first contains the eigenvalues, the
/// second all the eigenvectors in a single vector.
pub fn eigh_tridiagonal(a: &mut TridiagonalMatrix, iu: i32) -> (Vec<f64>, Vec<f64>) {
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
            &mut a.diag,
            &mut a.off,
            0.0,
            0.0,
            1,
            iu,
            &mut m,
            &mut w,
            &mut z,
            a.n as i32,
            &nzc,
            &mut isuppz,
            &mut tryrac,
            &mut work,
            lwork as i32,
            &mut iwork,
            liwork as i32,
            &mut info,
        );
    }
    // check for success
    assert_eq!(info, 0);
    // we only want the first iu eigenvalues
    w.truncate(iu as usize);
    (w, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eigh() {
        let mut a = TridiagonalMatrix {
            n: 5,
            diag: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            off: vec![0.0; 5],
        };
        let e_vals;
        let e_vecs;
        (e_vals, e_vecs) = eigh_tridiagonal(&mut a, 5);
        println!("{:?}", e_vals);
        println!("{:?}", e_vecs);
    }
}
