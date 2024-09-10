extern crate rand;
extern crate csv;

use rand::Rng;
use rand_distr::{Normal, Distribution};
use std::fs::File;
use csv::Writer;

const A: f64 = 0.8935e-9;                  // distance between sites (m)
const V_0: f64 = 1e12;                  // base hopping attemp rate (eV)

//**!? in papers stated as the "inverse localization radius" where 2ya = 10, this results in y = 5.596*10^9, however they have y listed has 1/5.596*10^-9
// const GAMMA: f64 = 10.0/(2.0 * A);      // satisfy 2*gamma*a = 10
const GAMMA: f64 = 1.0/(5.596e-9);
const K_B: f64 = 8.617333262e-5;        // Boltzmann constant (eV/K)
const T: f64 = 295.0;                   // Temperature (kelvin)

struct Lattice {
    size: usize,
    site_energy: Vec<Vec<Vec<f64>>>,
    site_state: Vec<Vec<Vec<u8>>>,
    temp: f64,
}

impl Lattice {
    fn new(size: usize, temp: f64) -> Self {
        let site_energy= generate_energy_levels(size);
        let site_state = vec![vec![vec![1; size]; size]; size];
        Lattice { size, site_energy, site_state, temp }
    }

    // check if the value in the site_state vetor is 0 or 1 for the xyz coordinates
    fn is_valid_site(&self, x: usize, y: usize, z: usize) -> bool {
        if x < self.size && y < self.size && z < self.size {
            self.site_state[x][y][z] == 0 || self.site_state[x][y][z] == 1
        } else {
            false
        }
    }

    // hopping probability (v_ij) for jumping from one site to another
    fn charge_hopping_rate(&self, x1: usize, y1: usize, z1: usize, x2: usize, y2: usize, z2: usize) -> f64 {
        let target_site_distance: f64 = self.euclidean_distance(x1, y1, z1, x2, y2, z2);

        // find difference in energy between the two sites
        let delta_e = self.site_energy[x2][y2][z2] - self.site_energy[x1][y1][z1];
        // let delta_e = -1.0;
        let energy_factor = if delta_e > 0.0 {
            f64::exp((-delta_e)/(K_B*self.temp))
        } else {
            1.0
        };
        
        V_0 * (-2.0 * GAMMA * target_site_distance).exp() * energy_factor
    }
    
    // Euclidean distance between two sites
    fn euclidean_distance(&self, x_cur: usize, y_cur: usize, z_cur: usize, x_tgt: usize, y_tgt: usize, z_tgt: usize) -> f64 {

        let dx: f64 = x_tgt as f64 - x_cur as f64;
        let dy: f64 = y_tgt as f64 - y_cur as f64;
        let dz: f64 = z_tgt as f64 - z_cur as f64;
        
        A * (dx * dx + dy * dy + dz * dz).sqrt()    // euclidean distance factor
    }

    // calculate hopping probabilites for all possible sites
    fn simulate_hop(&self, current_site: (u8, u8, u8), grid_size: u8) -> ((u8, u8, u8), f64, f64) {
        let (x0, y0, z0) = current_site;

        let mut site_coord: Vec<(u8, u8, u8)> = vec![];
        let mut site_prob: Vec<f64> = vec![];
        let mut event_probabilities: Vec<f64> = vec![];

        println!("{:?}", current_site);
        for x in (x0.saturating_sub(grid_size))..=x0+grid_size {
            println!("{}", x);
            // for y in y0-grid_size..=y0+grid_size {
            for y in (y0.saturating_sub(grid_size))..=y0+grid_size {
                println!("{}", y);
                for z in (z0.saturating_sub(grid_size))..=z0+grid_size {
                    println!("{}", z);
                    // Skip the current site location and sites out of lattice boundary
                    if (x == x0 && y == y0 && z == z0) || !self.is_valid_site(x as usize,y as usize,z as usize) {
                        println!("Skipped site x:{} y:{} z:{} site does not exist, or is current location", x,y,z);
                        continue;
                    }

                    let prob = self.charge_hopping_rate(x0.into(), y0.into(), z0.into(), x.into(), y.into(), z.into());

                    site_coord.push((x,y,z));
                    site_prob.push(prob);

                }            
            }
        }

        //sum of total event probabilities 
        let vtot: f64 = site_prob.iter().sum();

        //get normalised hop probabilities
        for &prob in site_prob.iter() {
            let event_prob = prob / vtot;
            event_probabilities.push(event_prob);
        }
        
        //random event selection
        let mut event_selector = rand::thread_rng();                                 // declare random number generator
        let random_site_index: f64 = event_selector.gen_range(0.0..=1.0) + f64::EPSILON;        // randomise an index of event to execute

        // ^^ *? Epsilon the best way to not allow 0? would also mean site_index could go beyond 1
        //       however it does get caught by the "i < event_probabilities.len()" in while loop providing the same outcome as if i > 1 is floored

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

        let t_0 = 1.0/(6.0 * V_0 * (-2.0 * GAMMA * A).exp());

        let time_increment = time_elapsed/t_0;                // normalised to periods of t_0
        let x_displacement = (new_site.0 as f64 - x0 as f64) * A;

        return (new_site, time_increment, x_displacement)            // return new site, elapsed time, and displacement in the x vector
    }

}

fn main() {
    // declare lattice of size and temp
    let lattice = Lattice::new(60, T);

    // **? what determines starting point - random selection
    // let log_vector = simulate_hops(lattice,(20,20,20), 10, 100.0);

    // println!("{:?}",lattice.simulate_hop((x0,y0,z0), 3));
    lattice.simulate_hop((0,0,0), 3);
    // lattice.simulate_hop((30,30,30), 3);

    //Export the logged data
    // export_data(log_vector)
}

fn simulate_hops(lattice: Lattice, random_site: (u8, u8, u8), snapshot_freq: usize, runtime: f64) -> Vec<(f64, f64)> {
    // Declare the lattice
    let mut t0_elapsed_total = 0.0;
    let mut x_displacement_total = 0.0;
    let mut log_vector: Vec<(f64, f64)> = Vec::new(); // Vector to store time and displacement
    let mut iteration_counter = 0;

    let (mut x0, mut y0, mut z0) = random_site;

    while t0_elapsed_total < runtime {
        println!("MAIN: Simulating Hop: {:?}", (x0,y0,z0));
        let (new_site, t0_elapsed_inc, x_displacement_inc) = lattice.simulate_hop((x0,y0,z0), 3);

        (x0, y0, z0) = new_site;
        t0_elapsed_total += t0_elapsed_inc.log10();
        x_displacement_total += x_displacement_inc;

        iteration_counter += 1;

        // Store every Nth snapshot in the vector
        if iteration_counter % snapshot_freq == 0 {
            log_vector.push((t0_elapsed_total, x_displacement_total));
        }


        // println!("TIME LOG BASE: {}", t0_elapsed_total)
        // println!("MAIN: new site: {:?}", new_site);

    }

    return log_vector
}

// Generates Guassian/Normal Distribution of starting energy levels for each site
fn generate_energy_levels(size: usize) -> Vec<Vec<Vec<f64>>> {
    let std_dev: f64 = 162e-3; // eV
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

// Write Lattice data to a CSV file
fn export_data(log_vector: Vec<(f64, f64)>) {
    //write data to csv
    let file = File::create("log_data.csv").expect("Unable to create file");
    let mut wtr = Writer::from_writer(file);
    
    // Write headers (optional)
    // wtr.write_record(&["Time", "Distance"]).expect("Unable to write record");

    for (time, distance) in log_vector.iter() {
        wtr.write_record(&[time.to_string(), distance.to_string()]).expect("Unable to write record");
    }

    wtr.flush().expect("Unable to flush CSV writer");
    println!("export successful to: log_data.csv");
}