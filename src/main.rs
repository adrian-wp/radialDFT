mod integrate;
mod linalg;
mod dft;
mod functionals;
mod utils;

fn main() {
    let z = 36;
    dft::single_atom_dft(z, 1e-4, 100.0, 3001, utils::default_occupation(z), 5,
                         100, 5, 1e-5, 1e-3, functionals::vwn_xc);
}
