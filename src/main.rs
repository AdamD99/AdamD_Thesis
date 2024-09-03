extern crate rand;
extern crate csv;

use rand::Rng;
use rand_distr::{Normal, Distribution};
use std::fs::File;
use csv::Writer;


struct Lattice {
    size: usize,
    site_energy: Vec<Vec<Vec<f64>>>,
    site_state: Vec<Vec<Vec<u8>>>,
}

impl Lattice {
    fn new(size: usize) -> Self {
        let site_energy= generate_energy_levels(size);
        let site_state = vec![vec![vec![1; size]; size]; size];
        Lattice { size, site_energy, site_state }
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
        let a: f64 = 0.6e-9;                // distance between sites       //** since uniform can be made obsolete?
        let v_0: f64 = 1e12;               // base hopping attemp rate
        let gamma: f64 = 10.0/(2.0 * a);    // satisfy 2*gamma*a = 10
        const K_B: f64 = 1.380649e-23; // J/K Boltzmann constant
        const T: f64 = 295.0; //kelvin
        let target_site_distance: f64 = a * self.euclidean_factor(x1, y1, z1, x2, y2, z2);

        // find difference in energy between the two sites
        let delta_e = self.site_energy[x2][y2][z2] - self.site_energy[x1][y1][z1];
        // let delta_e = -1.0;
        let energy_factor = if delta_e > 0.0 {
            f64::exp((-delta_e)/(K_B*T))
        } else {
            1.0
        };
        
        v_0 * (-2.0 * gamma * target_site_distance).exp() * energy_factor
    }
    
    // Euclidean distance between two sites
    fn euclidean_factor(&self, x_cur: usize, y_cur: usize, z_cur: usize, x_tgt: usize, y_tgt: usize, z_tgt: usize) -> f64 {

        let dx: f64 = x_tgt as f64 - x_cur as f64;
        let dy: f64 = y_tgt as f64 - y_cur as f64;
        let dz: f64 = z_tgt as f64 - z_cur as f64;
        
        (dx * dx + dy * dy + dz * dz).sqrt()    // euclidean distance factor
    }

    // calculate hopping probabilites for all possible sites
    fn simulate_hop(&self, current_site: (u8, u8, u8), grid_size: u8) {
        let (x0, y0, z0) = current_site;

        let mut site_coord: Vec<(u8, u8, u8)> = vec![];
        let mut site_prob: Vec<f64> = vec![];
        let mut event_probabilities: Vec<f64> = vec![];

        for x in x0-grid_size..=x0+grid_size {
            for y in y0-grid_size..=y0+grid_size {
                for z in z0-grid_size..=z0+grid_size {
                    // Skip the current site location and sites out of lattice boundary
                    if (x == x0 && y == y0 && z == z0) || !self.is_valid_site(x as usize,y as usize,z as usize) {
                        println!("Skipped site x:{} y:{} z:{} site does not exist, or is current location", x,y,z);
                        continue;
                    }

                    let prob = self.charge_hopping_rate(x0.into(), y0.into(), z0.into(), x.into(), y.into(), z.into());

                    site_coord.push((x,y,z));
                    site_prob.push(prob);

                    // println!("hopping prob: {}", self.charge_hopping_rate(x0.into(), y0.into(), z0.into(), x.into(), y.into(), z.into()));
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
        let random_site_index: f64 = event_selector.gen_range(0.0..=1.0) + f64::EPSILON;    // randomise an index of event to execute

        // ^^ *? Epsilon the best way to not allow 0? would also mean site_index could go beyond 1
        //       however it does get caught by the "i < event_probabilities.len()" in while loop providing the same outcome as if i > 1 is floored

        println!("site index: {}", random_site_index);

        //execute hopping event for randomly selected site
        let mut current_probability_total = event_probabilities[0]; // initialise total to first entry in event probs
        let mut i = 0;

        while current_probability_total < random_site_index  && i < event_probabilities.len() {
            i += 1;
            current_probability_total += event_probabilities[i];
        }

        // ^^ *? All of the site_prob and event_probabilities for sites where delta_e is > 0, are equal to 0

        let new_site = site_coord[i];
        
        let random_time_index: f64 = event_selector.gen_range(0.0..=1.0) + f64::EPSILON; 
        println!("time index: {}", random_time_index);
        let time_increment = -random_time_index.ln() / vtot;

        let a: f64 = 0.6e-9;                // distance between sites       //** since uniform can be made obsolete?
        let v_0: f64 = 1e12;               // base hopping attemp rate
        let gamma: f64 = 10.0/(2.0 * a);    // satisfy 2*gamma*a = 10

        let t_0 = 1.0/(6.0 * v_0 * (-2.0 * gamma * a).exp());

        println!("Target Hopping Site: {:?}", new_site);
        println!("time incremented: {}", time_increment);

        println!("time scaler: {}", t_0);
        println!("scaled time: {}", time_increment/t_0);

        // ^^* time increment is huge
        

        // println!("hopping prob: {} Total", _vtot);
        // println!("normalised hopping prob: {:?}", event_probabilities);
        // println!("{:?}", site_coord);
        // println!("{:?}", site_prob);
    }

}

fn main() {
    // Declare the lattice
    let lattice = Lattice::new(60);

    // **? what determines starting point

    // let mut hop_count = 10;
    // while hop_count > 0 {
    //     lattice.simulate_hop((20,20,20), 3);
    //     hop_count -= 1;
    // }

    lattice.simulate_hop((20,20,20), 3);
    lattice.simulate_hop((30,30,30), 3);

    //Export the site_state data
    //export_data(&lattice)
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
fn export_data(lattice: &Lattice) {
    //write data to csv
    let file = File::create("lattice.csv").expect("Unable to create file");
    let mut wtr = Writer::from_writer(file);
    
    for x in 0..(lattice.size) {
        for y in 0..(lattice.size) {
            for z in 0..(lattice.size) {
                wtr.write_record(&[x.to_string(), y.to_string(), z.to_string(), lattice.site_state[x][y][z].to_string()])
                    .expect("Unable to write record");
            }
        }
    }

    wtr.flush().expect("Unable to flush CSV writer");
    println!("export successful to: lattice.csv");
}