use itertools::Itertools;


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
