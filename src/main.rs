mod integrate;
mod linalg;
mod dft;
mod functionals;
mod utils;

fn main() -> std::io::Result<()> {
    let z_max = 36;
    let start = std::time::Instant::now();
    let mut results = Vec::with_capacity(z_max);
    for z in 1..z_max + 1 {
        let mut result = dft::single_atom_dft(z as i32, 1e-4, 100.0, 5001,
                                              utils::default_occupation(z as i32), 5,
                                              100, 5, 1e-5, 1e-4,
                                              functionals::vwn_xc, 1);
        // clear unnecessary results before saving in vector
        result.density = vec![];
        result.orbitals.eigenvectors = [vec![], vec![], vec![], vec![]];
        result.grid = dft::LogGrid { n: 0, x: vec![], r: vec![], h_x: 0.0 };
        results.push(result);
    }
    println!("{:?} elapsed.", start.elapsed());
    let z_range = (1..z_max + 1).collect();
    utils::write_energies("energies/energies_rust.csv", &results, &z_range);
    Ok(())
}
