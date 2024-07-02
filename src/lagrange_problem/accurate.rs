use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use itertools::izip;


pub struct Task {
    piston_mass: f64,
    tube_len: f64,
    tube_diameter: f64,
    gas_const: f64,
    gas_k: f64,
    p0: f64,
    temper0: f64,
    x0: f64
}


impl Task {
    pub fn new(piston_mass: f64,
           tube_len: f64,
           tube_diameter: f64,
           gas_const: f64,
           gas_k: f64,
           p0: f64,
           temper0: f64,
           x0: f64) -> Task
    {
        Task {
            piston_mass,
            tube_len,
            tube_diameter,
            gas_const,
            gas_k,
            p0,
            temper0,
            x0
        }
    }
}


pub struct Solution
{
    t: Vec<f64>,
    x: Vec<f64>,
    v: Vec<f64>
}


impl Solution
{
    fn new(t: Vec<f64>,
           x: Vec<f64>,
           v: Vec<f64>) -> Self
    {
        Self {
            t, x, v
        }
    }

    pub fn time(&self) -> &Vec<f64>
    {
        &self.t
    }

    pub fn coordinate(&self) -> &Vec<f64>
    {
        &self.x
    }

    pub fn velocity(&self) -> &Vec<f64>
    {
        &self.v
    }

    pub fn save_csv(&self, path: &Path)
    {
        let display = path.display();
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file
        };
        file.write("t,x,v\n".as_bytes()).expect("writing error");
        for (ti, xi, vi) in izip!(&self.t, &self.x, &self.v) {
            write!(&mut file, "{},{},{}\n", ti, xi, vi)
                .expect("writing error");
        }
    }
}


pub fn solve(task: &Task, n: usize) -> Solution
{
    let k = task.gas_k;

    let rho0 = task.p0 / (task.temper0 * task.gas_const);
    let s = 0.25 * PI * task.tube_diameter * task.tube_diameter;
    let m_gas = rho0 * s * task.x0;
    let c0 = (k * task.p0 / rho0).sqrt();
    let gamma = (k + 1.) / (2.*k);
    let m_ratio = m_gas / task.piston_mass;
    let tau = task.x0 / c0;

    let bracket = move |t: f64| {
        1. + gamma*m_ratio * t/tau
    };
    let fx = move |t: f64| {
        task.x0 + 2.*c0*t/(k - 1.) + 2.*k/((k - 1.)*m_ratio)*task.x0 * (
            1. - bracket(t).powf(2./(k + 1.))
        )
    };
    let fv = move |t: f64| {
        2.*c0/(k - 1.) * (
            1. - bracket(t).powf((1. - k)/(k + 1.))
        )
    };

    let mut t = vec![0.0; n];
    let duration = tau*(2. + gamma*m_ratio);
    let dt = duration / (n as f64);
    for i in 1..n {
        t[i] = (i as f64) * dt;
    }
    let mut x = Vec::from_iter(t.iter().map(|&t| fx(t)));
    let i = x.iter().position(|&xi| xi > task.tube_len);
    if i != None {
        let indx = i.unwrap();
        x = Vec::from(&x[..indx]);
        t = Vec::from(&t[..indx]);
    }
    let v = Vec::from_iter(
        t.iter().map(|&t| fv(t))
    );

    Solution::new(t, x, v)
}
