extern crate rand;
extern crate csv;

use rand::Rng;
use rand_distr::{Normal, Distribution};
use core::time;
use std::fs;
use csv::Writer;
use chrono::Local;

const A: f64 = 0.8935e-9;                  // distance between sites (m)
const V_0: f64 = 1.0e12;                  // base hopping attemp rate (eV)

//**!? in papers stated as the "inverse localization radius" where 2ya = 10, this results in y = 5.596*10^9, however they have y listed has 1/5.596*10^-9
// const GAMMA: f64 = 10.0/(2.0 * A);      // satisfy 2*gamma*a = 10
const GAMMA: f64 = 5.596e-9;            // m-1
const K_B: f64 = 8.617333262e-5;        // Boltzmann constant (eV/K)
const MATERIAL_LENGTH: f64 = 100e-9;            // 100nm material

struct Lattice {
    size: usize,                        // size of the Cubic lattice
    site_energy: Vec<Vec<Vec<f64>>>,    // vector containing energy level for each site
    site_state: Vec<Vec<Vec<u8>>>,      // vector containing information of site status
    temp: f64,                          // Temperature (kelvin)
    field: f64,                         // volts/m
}

impl Lattice {
    fn new(size: usize, temp: f64, field: f64, integrity: f64) -> Self {
        let site_energy= generate_energy_levels(size);
        let site_state = generate_site_states(size, integrity);
        Lattice { size, site_energy, site_state, temp, field}
    }

    // check if the value in the site_state vetor is 0 or 1 for the xyz coordinates
    fn is_valid_site(&self, x: usize, y: usize, z: usize) -> bool {
        if x < self.size && y < self.size && z < self.size {
            self.site_state[x][y][z] == 1
        } else {
            false
        }
    }

    // hopping probability (v_ij) for jumping from one site to another
    fn charge_hopping_rate(&self, (x1,y1,z1): (usize,usize,usize), (x2,y2,z2): (usize,usize,usize), target_site_distance: f64, delta_v: f64) -> f64 {

        // println!("site energies target {} origin {}, delta v {}, target dist {}", self.site_energy[x2][y2][z2], self.site_energy[x1][y1][z1], delta_v, target_site_distance);

        // find difference in energy between the two sites
        let delta_e = self.site_energy[x2][y2][z2] - self.site_energy[x1][y1][z1] + delta_v;

        // println!("delta_e {}", delta_e);

        // let delta_e = -1.0;
        let energy_factor = if delta_e > 0.0 {
            f64::exp((-delta_e)/(K_B*self.temp))
        } else {
            1.0
        };

        // println!("energy factor: {}", energy_factor);
        // println!("gamma {}, V_O {}, target dist {}", GAMMA, V_0, target_site_distance);

        V_0 * (-2.0 * GAMMA * target_site_distance).exp() * energy_factor
    }

    // MARK: Simulate Hop
    // calculate hopping probabilites for all possible sites
    fn simulate_hop(&self, current_site: (u8, u8, u8), grid_size: i8) -> ((u8, u8, u8), f64, f64) {
        let (x0, y0, z0) = current_site;
        
        let number_of_neighbours = (2 * grid_size as usize + 1).pow(3) - 1; // exclude the current site
        let mut site_coord: Vec<(u8, u8, u8)> = Vec::with_capacity(number_of_neighbours);
        let mut site_prob: Vec<f64> = Vec::with_capacity(number_of_neighbours);

        for x in -grid_size..=grid_size {
            let new_x: u8 = (x0 as i8 + x).rem_euclid(self.size as i8) as u8;
            for y in -grid_size..=grid_size {
                let new_y: u8 = (y0 as i8 + y).rem_euclid(self.size as i8) as u8;
                for z in -grid_size..=grid_size {
                    let new_z: u8 = (z0 as i8 + z).rem_euclid(self.size as i8) as u8;

                    // Skip the current site location and sites out of lattice boundary
                    if ( new_x == x0 && new_y == y0 && new_z == z0 ) || !self.is_valid_site(new_x as usize,new_y as usize,new_z as usize) {
                        // println!("Skipped site x:{} y:{} z:{} site does not exist, or is current location", new_x,new_y,new_z);
                        continue;
                    }

                    // euclidian distance to target site
                    let target_site_distance: f64 = A * ((x * x + y * y + z * z) as f64).sqrt();
        
                    // voltage difference between sites
                    let delta_v = self.field * x as f64 * A;   // volts

                    let prob = self.charge_hopping_rate((x0.into(),y0.into(),z0.into()), (new_x.into(),new_y.into(),new_z.into()), target_site_distance, delta_v);

                    site_coord.push((new_x,new_y,new_z));
                    site_prob.push(prob);
                }
            }
        }

        //sum of total event probabilities 
        let vtot: f64 = site_prob.iter().sum();

        //get normalised hop probabilities
        let event_probabilities = site_prob.iter()
            .map(|prob| prob / vtot)
            .collect::<Vec<_>>(); // normalise probabilities
        
        //random event selection
        let mut event_selector = rand::thread_rng();                                 // declare random number generator
        let random_site_index: f64 = event_selector.gen_range(0.0..=1.0) + f64::EPSILON;        // randomise an index of event to execute

        //execute hopping event for randomly selected site
        let mut current_probability_total = event_probabilities[0]; // initialise total to first entry in event probs
        let mut i = 0;

        while current_probability_total < random_site_index  && i < event_probabilities.len() {
            i += 1;
            current_probability_total += event_probabilities[i];
        }

        let new_site = site_coord[i];
        
        let random_time_index: f64 = event_selector.gen_range(0.0..=1.0) + f64::EPSILON; 
        let time_elapsed = -random_time_index.ln() / vtot;  // seconds

        // rust cannot declare as global since exp() is a runtime variable and it cannot handle it.
        let t_0: f64 = 1.0/(6.0 * V_0 * (-2.0 * GAMMA * A).exp());

        let time_increment = time_elapsed/t_0;                // normalised to periods of t_0
        // println!("new site: {:?}, new sitex: {}, x0: {}, A: {}", new_site, new_site.0, x0, A);
        let x_difference = new_site.0 as f64 - x0 as f64;
        let x_difference_correction = (x_difference + 30.0) % 60.0 - 30.0;
        let x_displacement = x_difference_correction * A;

        // println!("time_elapsed {}, time_inc {}, x_displacement {}", time_elapsed, time_increment, x_displacement);

        return (new_site, time_increment, x_displacement)            // return new site, elapsed time, and displacement in the x vector
    }

}

// MARK: Main 
fn main() {
    // declare lattice of size, temp and integrity (1.0 = 100% all sites present, 0.0 = 0% no sites present)
    let size: u8 = 60;
    let grid_size = 3;
    let voltage = 1.0;                  // default voltage
    let temperature_var = 298.0;        // default temperature
    let material_integrity = 1.0;       // default integrity (1.0 = 100%)

    // temperature Variance
    // for temperature_var in (285..=315).step_by(5) {
    //     let field = voltage/MATERIAL_LENGTH;        // v/m 
    //     let lattice = Lattice::new(size as usize, temperature_var as f64, field, material_integrity);

    //     // random selection algorithm for starting site
    //     let mut rng = rand::thread_rng();
    //     let starting_site: (u8,u8,u8) = (rng.gen_range(0..size), rng.gen_range(0..size), rng.gen_range(0..size));

    //     let log_vector = simulate_hops(lattice,starting_site, grid_size, 1000, 100000.0);

    //     //Export the logged data
    //     export_data(log_vector, temperature_var, voltage, material_integrity)
    // }

    // electric field variance
    // let voltage_step = 0.5;
    // let mut voltage_var = 0.5;
    // while voltage_var < 3.1 {
    //     let field = voltage_var/MATERIAL_LENGTH;        // v/m 
    //     let lattice = Lattice::new(size as usize, temperature_var as f64, field, material_integrity);

    //     println!("{}", lattice.field);

    //     // random selection algorithm for starting site
    //     let mut rng = rand::thread_rng();
    //     let starting_site: (u8,u8,u8) = (rng.gen_range(0..size), rng.gen_range(0..size), rng.gen_range(0..size));

    //     let log_vector = simulate_hops(lattice,starting_site, grid_size, 1000, 100000.0);

    //     //Export the logged data
    //     export_data(log_vector, temperature_var as i32, voltage_var, material_integrity);

    //     voltage_var += voltage_step;
    // }

    //integtriy variance
    let integrity_decimal = [1.0, 1.0, 1.0];
    for material_integrity in integrity_decimal.iter() {
        let field = voltage/MATERIAL_LENGTH;        // v/m 
        let lattice = Lattice::new(size as usize, temperature_var, field, *material_integrity);

        // random selection algorithm for starting site
        let mut rng = rand::thread_rng();
        let starting_site: (u8,u8,u8) = (rng.gen_range(0..size), rng.gen_range(0..size), rng.gen_range(0..size));

        let log_vector = simulate_hops(lattice,starting_site, grid_size, 1000, 100000.0);

        //Export the logged data
        export_data(log_vector, temperature_var as i32, voltage, *material_integrity)
    }

    println!("Simulation Complete!")
}


fn simulate_hops(lattice: Lattice, random_site: (u8, u8, u8), grid_size: i8, snapshot_freq: usize, runtime: f64) -> Vec<(f64, f64)> {
    // Declare the lattice
    let mut t0_elapsed_total = 0.0;
    let mut x_displacement_total = 0.0;
    let mut log_vector: Vec<(f64, f64)> = Vec::new(); // Vector to store time and displacement
    let mut iteration_counter = 0;

    // For reporting simulation speed
    let mut timer = std::time::Instant::now();
    let mut percent_complete_reported = 5.0;
    let mut iterations_at_last_report = 0;

    let (mut x0, mut y0, mut z0) = random_site;

    while t0_elapsed_total < runtime {
        // println!("MAIN: Simulating Hop: {:?}", (x0,y0,z0));
        let (new_site, t0_elapsed_inc, x_displacement_inc) = lattice.simulate_hop((x0,y0,z0), grid_size);

        (x0, y0, z0) = new_site;
        t0_elapsed_total += t0_elapsed_inc;
        x_displacement_total += x_displacement_inc;

        iteration_counter += 1;

        // Store every Nth snapshot in the vector
        if iteration_counter % snapshot_freq == 0 {
            log_vector.push((t0_elapsed_total, x_displacement_total));

            let percent_complete = (t0_elapsed_total/runtime)*100.0;
            if percent_complete >= percent_complete_reported {
                let elapsed_time = timer.elapsed().as_secs_f64();
                let iterations = iteration_counter - iterations_at_last_report;
                println!("{:.0}% complete, simulation speed = {:.0} iter/sec", percent_complete, iterations as f64/elapsed_time);
                timer = std::time::Instant::now();
                percent_complete_reported += 5.0;
                iterations_at_last_report = iteration_counter;
            }
        }
        // println!("TIME LOG BASE: {}", t0_elapsed_total)
        // println!("MAIN: new site: {:?}", new_site);

    }
    println!("Average Displacement = {}", x_displacement_total/t0_elapsed_total);

    log_vector
}

// MARK: Gen. energy levels
// Generates Guassian/Normal Distribution of starting energy levels for each site
fn generate_energy_levels(size: usize) -> Vec<Vec<Vec<f64>>> {
    let std_dev: f64 = 25.68e-3; // eV
    let mean: f64 = 0.0;    //conduction band level for material
    let normal_dist = Normal::new(mean, std_dev).unwrap();
    let mut rng = rand::thread_rng();

    let mut site_energy = vec![vec![vec![0.0; size]; size]; size];

    for x in 0..size {
        for y in 0..size {
            for z in 0..size {
                site_energy[x][y][z] = normal_dist.sample(&mut rng);
            }
        }
    }

    site_energy
}

// MARK: Gen. site states
fn generate_site_states(size: usize, integrity: f64) -> Vec<Vec<Vec<u8>>> {
    let mut rng = rand::thread_rng();
    let mut site_states: Vec<Vec<Vec<u8>>> = vec![vec![vec![1; size]; size]; size];
    let (mut x, mut y, mut z);


    if integrity < 1.0 {

        let damage_sites: usize = ((size*size*size) as f64 * (1.0 - integrity)) as usize; // number of damaged sites

        // flag sites as damaged
        for _i in 0..damage_sites {
            
            //pick random number between 0 and size inclusive
            loop {
                x = rng.gen_range(0..size);
                y = rng.gen_range(0..size);
                z = rng.gen_range(0..size);

                // check that site is present
                if site_states[x][y][z] == 1 {
                    break
                }
            }

            // flag site as removed
            site_states[x][y][z] = 0;

        }
    }

    site_states
}

// TODO
// fn generate_site_distances(grid_size: i8) -> Vec<((i8,i8,i8), f64)> {
//     let number_of_neighbours = (2 * grid_size as usize + 1).pow(3) - 1; // exclude the current site
//     let mut site_distance_index: Vec<((i8,i8,i8), f64)> = Vec::with_capacity(number_of_neighbours);

//     for x in -grid_size..=grid_size {
//         for y in -grid_size..=grid_size {
//             for z in -grid_size..=grid_size {
//                 let site_distance = A * ((x * x + y * y + z * z) as f64).sqrt();

//                 site_distance_index.push(((x,y,z), site_distance));
//             }
//         }
//     }
    
//     site_distance_index
// }

// Write Lattice data to a CSV file
fn export_data(log_vector: Vec<(f64, f64)>, temp: i32, voltage: f64, integrity: f64) {
    // create data subfolder
    fs::create_dir_all("data").expect("Unable to create directory or directory already exists");

    let date = Local::now();
    let formatted_date = date.format("%d-%m-%Y-%H;%M;%S").to_string();
    let file_name = format!("data/log_data-Temp;{}_Voltage;{:.1}V_Integrity;{}%-{}.csv", temp as i32, voltage, (integrity*100.0) as i8, formatted_date);

    //write data to csv
    let file = fs::File::create(&file_name).expect("Unable to create file");
    let mut wtr = Writer::from_writer(file);
    
    // Write headers
    // wtr.write_record(&["Time", "Distance"]).expect("Unable to write record");

    for (time, distance) in log_vector.iter() {
        wtr.write_record(&[time.to_string(), distance.to_string()]).expect("Unable to write record");
    }

    wtr.flush().expect("Unable to flush CSV writer");
    println!("export successful to: {}", file_name);
}