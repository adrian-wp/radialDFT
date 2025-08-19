mod integrate;
mod linalg;
mod dft;
mod functionals;
mod utils;
mod io;

/// The main function of the executable. It parses the input arguments and then
/// calls dft::single_atom_dft() which performs the DFT calculations.
fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    // The two options for the arguments are specifying a config file, if the
    // first argument is "-c" or specifying some parameters via the arguments
    // in that case the first argument is the atomic number z
    if args.len() <= 1 {
        // no arguments -> print help
        println!(concat!(
            "No arguments provided.\n",
            "Usage: RadialDFT -c <config path>\n",
            "       RadialDFT <Z> [<r_min> <r_max> <n_grid> <output path>]\n",
            "Parameters\n",
            "    config path   path to toml file with parameters\n",
            "    Z             atomic number of atom\n",
            "    r_min         smallest value of the logarithmic grid (default: 1e-4)\n",
            "    r_max         largest value of the logarithmic grid (default: 100)\n",
            "    n_grid        number of grid points (default: 3001)\n",
            "    output path   path for file to save density and orbitals (optional)"
        ));
    } else if args[1] == "-c" {
        // Option 1: config file specified -> parse parameters from toml
        assert!(args.len() > 2, "Expected 2 arguments to specify config file");
        let config = io::read_config(args[2].as_str());
        // check whether z_range or z is specified, then run DFT
        if config.use_z_range {
            // z_range specified
            assert!(config.z_range[0] <= config.z_range[1],
                    "Start of z_range needs to be smaller than end");
            let mut results = Vec::new();
            let start = std::time::Instant::now();
            for z in config.z_range[0]..config.z_range[1] + 1 {
                // get occupations for current z or from string
                let occupations = if config.use_default_occupations {
                    utils::default_occupation(z)
                } else {
                    config.occupations.clone()
                };
                let mut result = dft::single_atom_dft(
                    z,
                    config.r_min,
                    config.r_max,
                    config.n_grid,
                    occupations,
                    config.mixing_steps,
                    config.max_iterations,
                    config.min_iterations,
                    config.e_tol,
                    config.rho_tol,
                    config.functional,
                    1,
                );
                // save density if a directory is specified
                if !config.densities_folder.is_empty() {
                    let current_path = format!("{}/{}.csv", config.densities_folder, z);
                    io::write_orbitals(current_path.as_str(), &result);
                }
                // clear unused vectors and add to results if energies will be saved later
                if !config.energies_path.is_empty() {
                    result.density = vec![];
                    result.orbitals.eigenvectors = [vec![], vec![], vec![], vec![]];
                    result.grid = dft::LogGrid {
                        n: 0,
                        x: vec![],
                        r: vec![],
                        h_x: 0.0,
                    };
                    results.push(result);
                }
            }
            println!("{:?} elapsed.", start.elapsed());
            if !config.energies_path.is_empty() {
                let z_vector = (config.z_range[0]..config.z_range[1] + 1).collect();
                io::write_energies(config.energies_path.as_str(), &results, &z_vector);
            }
        } else {
            // no z_range -> only calculate single z
            // use default occupations or from config
            let occupations = if config.use_default_occupations {
                utils::default_occupation(config.z)
            } else {
                config.occupations
            };
            // start calculation
            let start = std::time::Instant::now();
            let result = dft::single_atom_dft(
                config.z,
                config.r_min,
                config.r_max,
                config.n_grid,
                occupations,
                config.mixing_steps,
                config.max_iterations,
                config.min_iterations,
                config.e_tol,
                config.rho_tol,
                config.functional,
                2,
            );
            println!("{:?} elapsed.", start.elapsed());
            // save density/orbitals and energy if output path was specified
            if !config.density_path.is_empty() {
                io::write_orbitals(config.density_path.as_str(), &result);
            }
            if !config.energies_path.is_empty() {
                let results = vec![result];
                let z_vector = vec![config.z];
                io::write_energies(config.energies_path.as_str(), &results, &z_vector);
            }
        }
    } else {
        // Option 2: z, r_min, r_max, n_grid specified via args
        // parse values or use defaults
        let z = args[1].parse::<i32>().expect("Could not parse argument for Z to integer");
        let r_min: f64 = if args.len() > 2 {
            args[2].parse::<f64>().expect("Could not parse argument for r_min to float")
        } else { 1e-4 };
        let r_max: f64 = if args.len() > 3 {
            args[3].parse::<f64>().expect("Could not parse argument for r_max to float")
        } else { 100.0 };
        let n_grid: usize = if args.len() > 4 {
            args[4].parse::<usize>().expect("Could not parse argument for n_grid to integer")
        } else { 3001 };
        assert!(z >= 1 && z <= 118, "Z needs to be between 1 and 118");
        assert!(r_min.is_finite() && r_max.is_finite(),
                "r_min and r_max must be finite");
        assert!(r_min < r_max, "r_min must be smaller than r_max");
        assert!(n_grid > 2, "n_grid must be > 2");
        // start run
        let start = std::time::Instant::now();
        let result = dft::single_atom_dft(
            z,
            r_min,
            r_max,
            n_grid,
            utils::default_occupation(z),
            5,
            100,
            5,
            1e-5,
            1e-4,
            functionals::vwn_xc,
            2,
        );
        println!("{:?} elapsed.", start.elapsed());
        if args.len() > 5 {
            let path = args[5].clone();
            println!("Writing density and orbitals to {}", path);
            io::write_orbitals(path.as_str(), &result);
        }
    }
    Ok(())
}
