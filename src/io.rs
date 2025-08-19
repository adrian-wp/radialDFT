use serde::Deserialize;
use std::io::Read;
use std::io::Write;
use crate::functionals;
use crate::dft;
use crate::utils;

/// Represents a DFT run configuration as specified in a config.toml.
pub struct ConfigWithDefaults {
    pub use_z_range: bool,
    pub z: i32,
    pub z_range: [i32; 2],
    pub n_grid: usize,
    pub r_min: f64,
    pub r_max: f64,
    pub min_iterations: usize,
    pub max_iterations: usize,
    pub mixing_steps: usize,
    pub e_tol: f64,
    pub rho_tol: f64,
    pub density_path: String,
    pub energies_path: String,
    pub densities_folder: String,
    pub use_default_occupations: bool,
    pub occupations: dft::Occupations,
    pub functional: dft::XCFunctional,
}

// These structs are used to deserialize the toml file and are only used by
// the read_config() function

#[derive(Deserialize)]
struct Config {
    z: Option<i32>,
    z_range: Option<[i32; 2]>,
    occupations: Option<String>,
    functional: Option<String>,
    grid: Option<ConfigGrid>,
    scf: Option<ConfigSCF>,
    output: Option<ConfigOutput>,
}

#[derive(Deserialize)]
struct ConfigGrid {
    n: Option<usize>,
    r_min: Option<f64>,
    r_max: Option<f64>,
}

#[derive(Deserialize)]
struct ConfigSCF {
    min_iterations: Option<usize>,
    max_iterations: Option<usize>,
    mixing_steps: Option<usize>,
    e_tol: Option<f64>,
    rho_tol: Option<f64>,
}

#[derive(Deserialize)]
struct ConfigOutput {
    density: Option<String>,
    energies: Option<String>,
    densities_folder: Option<String>,
}

/// Reads and parses the configuration toml file at the given path. Default
/// values are applied if a values is not present.
pub fn read_config(path: &str) -> ConfigWithDefaults {
    let mut file = std::fs::File::open(path).expect("Could not open config file");
    let mut toml_str = String::new();
    file.read_to_string(&mut toml_str).expect("Could not read config file");
    let config: Config = toml::from_str(&toml_str).expect("Could not parse config file");

    // check if z or z_range is present
    let mut use_z_range = false;
    let mut z = 0;
    let mut z_range = [0, 0];
    if config.z_range.is_some() {
        use_z_range = true;
        z_range = config.z_range.unwrap();
    } else if config.z.is_some() {
        z = config.z.unwrap();
    } else {
        panic!("Either z or z_range must be provided");
    }

    // select the functional, use VWN by default
    let functional = match config
        .functional
        .unwrap_or("VWN".to_string())
        .to_ascii_lowercase()
        .as_str()
    {
        "vwn" => functionals::vwn_xc,
        "slater" => functionals::slater_exchange,
        _ => {
            panic!("Unknown functional");
        }
    };

    // resolve occupations string if specified
    let occupations_str = config.occupations.unwrap_or("default".to_string());
    let use_default_occupations = occupations_str == "default";
    let occupations = if use_default_occupations {
        dft::Occupations {
            f: [vec![], vec![], vec![], vec![]],
            max_l: 0,
        }
    } else {
        utils::resolve_occupation_string(occupations_str.as_str())
    };

    // get values or use default
    let grid_config = config.grid.unwrap_or(ConfigGrid {
        n: None,
        r_min: None,
        r_max: None,
    });
    let scf_config = config.scf.unwrap_or(ConfigSCF {
        min_iterations: None,
        max_iterations: None,
        mixing_steps: None,
        e_tol: None,
        rho_tol: None,
    });
    let output_config = config.output.unwrap_or(ConfigOutput {
        density: None,
        energies: None,
        densities_folder: None,
    });
    let n_grid = grid_config.n.unwrap_or(3001);
    let r_min = grid_config.r_min.unwrap_or(1e-4);
    let r_max = grid_config.r_max.unwrap_or(100.0);
    let max_iterations = scf_config.min_iterations.unwrap_or(5);
    let min_iterations = scf_config.max_iterations.unwrap_or(100);
    let mixing_steps = scf_config.mixing_steps.unwrap_or(5);
    let e_tol = scf_config.e_tol.unwrap_or(1e-5);
    let rho_tol = scf_config.rho_tol.unwrap_or(1e-4);
    let density_path = output_config.density.unwrap_or_default();
    let energies_path = output_config.energies.unwrap_or_default();
    let densities_folder = output_config.densities_folder.unwrap_or_default();
    ConfigWithDefaults {
        use_z_range,
        z,
        z_range,
        n_grid,
        r_min,
        r_max,
        min_iterations,
        max_iterations,
        mixing_steps,
        e_tol,
        rho_tol,
        density_path,
        energies_path,
        densities_folder,
        use_default_occupations,
        occupations,
        functional,
    }
}

/// Writes the energies and atomic numbers of the given DFTResults into a
/// csv-file at the specified path.
pub fn write_energies(path: &str, results: &Vec<dft::DFTResult>, z: &Vec<i32>) {
    debug_assert_eq!(results.len(), z.len());
    let mut file = std::fs::File::create(path).expect("Failed to create output file");
    writeln!(file, "z,iterations,e_total,e_kinetic,e_hartree,e_nucleus,e_xc,eigenvalues")
        .expect("Failed to write output file");
    for i in 0..results.len() {
        let e = &results[i].energy;
        writeln!(file, "{},{},{},{},{},{},{},\"{:?}\"", z[i], results[i].iterations, &e.total,
                 &e.kinetic, &e.hartree, &e.external, &e.xc, results[i].orbitals.eigenvalues)
            .expect("Failed to write output file");
    }
}

/// Writes the grid, density and orbitals of the DFTResult into a csv file.
pub fn write_orbitals(path: &str, result: &dft::DFTResult) {
    // verify length of all vectors
    debug_assert_eq!(result.grid.n, result.grid.r.len());
    debug_assert_eq!(result.grid.n, result.density.len());
    for l in 0..4 {
        debug_assert_eq!(result.orbitals.eigenvectors[l].len() * result.grid.n,
                         result.occupations.f[l].len());
    }
    // move orbital vectors into more convenient format, create header for orbitals
    let mut orbitals = Vec::new();
    let mut orbital_header = String::new();
    for l in 0..4 {
        orbitals.extend(result.orbitals.eigenvectors[l].chunks_exact(result.grid.n));
        for i in 0..result.occupations.f[l].len() {
            orbital_header.push_str(",");
            orbital_header.push_str((l + i + 1).to_string().as_str());
            orbital_header.push_str(match l {
                0 => "s",
                1 => "p",
                2 => "d",
                3 => "f",
                _ => "?"
            });
        }
    }
    // create file
    let mut file = std::fs::File::create(path)
        .expect(format!("Failed to create output file: {}", path).as_str());
    writeln!(file, "r,density{}", orbital_header)
        .expect(format!("Failed to write to output file: {}", path).as_str());
    // write data
    for i in 0..result.grid.n {
        // list orbital values
        let mut orbital_string = String::new();
        for j in 0..orbitals.len() {
            orbital_string.push_str(",");
            orbital_string.push_str(orbitals[j][i].to_string().as_str());
        }
        // write line
        writeln!(file, "{},{}{}", result.grid.r[i], result.density[i], orbital_string)
            .expect(format!("Failed to write to output file: {}", path).as_str());
    }
}
