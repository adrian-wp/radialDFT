//! This module contains helper functions mainly for obtaining the occupations
//! in the format used by the DFT code.

use crate::dft;

pub fn default_occupation(mut z: i32) -> dft::Occupations {
    debug_assert!(0 <= z && z <= 118);
    let f = match z {
        // exceptions
        // d-block
        24 => { [vec![2, 2, 2, 1], vec![6, 6], vec![5], vec![]] }
        29 => { [vec![2, 2, 2, 1], vec![6, 6], vec![10], vec![]] }
        41 => { [vec![2, 2, 2, 2, 1], vec![6, 6, 6], vec![10, 4], vec![]] }
        42 => { [vec![2, 2, 2, 2, 1], vec![6, 6, 6], vec![10, 5], vec![]] }
        44 => { [vec![2, 2, 2, 2, 1], vec![6, 6, 6], vec![10, 7], vec![]] }
        45 => { [vec![2, 2, 2, 2, 1], vec![6, 6, 6], vec![10, 8], vec![]] }
        46 => { [vec![2, 2, 2, 2], vec![6, 6, 6], vec![10, 10], vec![]] }
        47 => { [vec![2, 2, 2, 2, 1], vec![6, 6, 6], vec![10, 10], vec![]] }
        78 => { [vec![2, 2, 2, 2, 2, 1], vec![6, 6, 6, 6], vec![10, 10, 9], vec![14]] }
        79 => { [vec![2, 2, 2, 2, 2, 1], vec![6, 6, 6, 6], vec![10, 10, 10], vec![14]] }
        103 => {
            [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6, 1], vec![10, 10, 10], vec![14, 14]]
        }
        // f-block
        57 => { [vec![2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6], vec![10, 10, 1], vec![]] }
        58 => { [vec![2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6], vec![10, 10, 1], vec![1]] }
        64 => { [vec![2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6], vec![10, 10, 1], vec![7]] }
        89 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 1], vec![14]] }
        90 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 2], vec![14]] }
        91 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 1], vec![14, 2]] }
        92 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 1], vec![14, 3]] }
        93 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 1], vec![14, 4]] }
        96 => { [vec![2, 2, 2, 2, 2, 2, 2], vec![6, 6, 6, 6, 6], vec![10, 10, 10, 1], vec![14, 7]] }
        // use Aufbau principle if no exception
        _ => {
            let mut f = [vec![], vec![], vec![], vec![]];
            // iterate over subshells (n, l) in Aufbau principle order until no electrons left
            let mut n = 1;
            let mut l = 0;
            while z > 0 {
                // fit as many electrons as possible into current subshell
                let shell_capacity = 2 * (2 * l + 1);
                if shell_capacity < z {
                    f[l as usize].push(shell_capacity);
                    z -= shell_capacity;
                } else {
                    f[l as usize].push(z);
                    z = 0;
                }
                // next subshell
                l -= 1;
                n += 1;
                if l < 0 {
                    l = (n - 1) / 2;
                    n = (n + 2) / 2;
                }
            }
            f
        }
    };
    let max_l = f.iter().filter(|vec| { !vec.is_empty() }).count();
    dft::Occupations { f, max_l }
}

pub fn resolve_occupation_string(s: &str) -> dft::Occupations {
    let mut stack: Vec<&str> = s.split_whitespace().rev().collect();
    let mut f = [vec![], vec![], vec![], vec![]];
    while !stack.is_empty() {
        let s = stack.pop().unwrap();
        match s {
            "[He]" => { stack.extend(["1s2"]) }
            "[Ne]" => { stack.extend(["2p6", "2s2", "[He]"]) }
            "[Ar]" => { stack.extend(["3p6", "3s2", "[Ne]"]) }
            "[Kr]" => { stack.extend(["4p6", "3d10", "4s2", "[Ar]"]) }
            "[Xe]" => { stack.extend(["5p6", "4d10", "5s2", "[Kr]"]) }
            "[Rn]" => { stack.extend(["6p6", "5d10", "4f14", "6s2", "[Xe]"]) }
            &_ => {
                let n: usize = s[0..1].parse().unwrap();
                let l = match &s[1..2] {
                    "s" => { 0 }
                    "p" => { 1 }
                    "d" => { 2 }
                    "f" => { 3 }
                    _ => { panic!("invalid subshell for l") }
                };
                let e = s[2..].parse().unwrap();
                while f[l].len() < n - l {
                    f[l].push(0);
                }
                f[l][n - l - 1] = e;
            }
        }
    }
    let max_l = f.iter().filter(|vec| { !vec.is_empty() }).count();
    dft::Occupations { f, max_l }
}

/// used by some tests to compare two vectors of floats
#[allow(unused)]
pub fn almost_equal_vec_f64(a: &Vec<f64>, b: &Vec<f64>, epsilon: f64) {
    assert_eq!(a.len(), b.len());
    for i in 0..a.len() {
        if a[i] != b[i] {
            assert!((a[i] - b[i]).abs() < epsilon,
                    "Vectors are not (almost) equal at position {}", i);
        }
    }
}
