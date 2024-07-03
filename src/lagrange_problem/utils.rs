use std::fs::File;
use std::io::Write;
use std::path::Path;

use itertools::{izip, Itertools};


pub struct Results
{
    t: Vec<f64>,
    p_piston: Vec<f64>,
    p_bottom: Vec<f64>,
    piston_v: Vec<f64>,
    piston_x: Vec<f64>
}


impl Results
{
    pub fn new(t: Vec<f64>,
               p_piston: Vec<f64>,
               p_bottom: Vec<f64>,
               piston_v: Vec<f64>,
               piston_x: Vec<f64>) -> Self
    {
        Self {
            t, p_piston, p_bottom, piston_v, piston_x
        }
    }

    pub fn time(&self) -> &Vec<f64>
    {
        &self.t
    }

    pub fn pressure_on_piston(&self) -> &Vec<f64>
    {
        &self.p_piston
    }

    pub fn pressure_on_bottom(&self) -> &Vec<f64>
    {
        &self.p_bottom
    }

    pub fn piston_speed(&self) -> &Vec<f64>
    {
        &self.piston_v
    }

    pub fn piston_coordinate(&self) -> &Vec<f64>
    {
        &self.piston_x
    }

    pub fn save_csv(&self, path: &Path)
    {
        let display = path.display();
        let mut file = match File::create(path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file
        };
        file.write_all("t,x,v,p_piston,p_bottom\n".as_bytes()).expect(
            "writing error"
        );
        for (t, x, v, pp, pb) in izip!(
                self.time(),
                self.piston_coordinate(),
                self.piston_speed(),
                self.pressure_on_piston(),
                self.pressure_on_bottom()
        )
        {
            writeln!(&mut file, "{},{},{},{},{}", t, x, v, pp, pb).expect(
                "writing error"
            );
        }
    }
}


pub fn linspace(start: f64, stop: f64, n: usize) -> (Vec<f64>, f64)
{
    let step = (stop - start) / ((n - 1) as f64);
    (
        (0..n).map(|i| start + (i as f64)*step).collect_vec(),
        step
    )
}


pub fn min(v: &[f64]) -> f64
{
    v.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}
