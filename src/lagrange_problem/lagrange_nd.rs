use std::f64::consts::PI;

use itertools::izip;
use ndarray::prelude::*;
use ndarray::Array;

use super::lagrange::Task;
use super::utils::{Results, min};


pub fn solve(task: &Task) -> Results
{
    let k = task.gas_k;
    let volumes = task.n_volumes;
    let nodes = task.n_volumes + 1;

    let mut mesh =
        Array::linspace(0., task.x0, nodes);
    let dx = mesh[1] - mesh[0];

    let mut u_borders =
        Array::zeros(nodes.f());
    let mut p =
        Array::from_elem(volumes.f(), task.p0);
    let rho0 = task.p0 / (task.gas_const * task.temper0);
    let mut rho =
        Array::from_elem(volumes.f(), rho0);
    let mut inner_e = (&p / &rho) / (k - 1.);
    
    let s = 0.25 * PI * task.tube_diameter * task.tube_diameter;
    let q =
        Array::from_elem(volumes.f(), rho0 * s * dx);
    let mut q_interface = 
        Array::zeros(nodes.f());
    if let Some(last) = q_interface.last_mut() {
        *last = task.piston_mass;
    }

    let mut t = 0.;
    let mut u =
        Array::zeros(volumes.f());

    let capacity = 1000;
    let mut t_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut p_piston_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut p_bottom_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut piston_x_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut piston_v_store: Vec<f64> = Vec::with_capacity(capacity);

    t_store.push(t);
    p_piston_store.push(task.p0);
    p_bottom_store.push(task.p0);
    piston_x_store.push(task.x0);
    piston_v_store.push(0.);

    while *mesh.last().unwrap() < task.tube_len {
        let mut c = k * (&p / &rho);
        c.mapv_inplace(f64::sqrt);
        let dt = calc_time_step(&mesh, &u, &c, task.cfl);
        let pl = p.slice(s![..p.len()-1]);
        let pr = p.slice(s![1..]);
        let delta_p = &pr - &pl;

        for (i, ub) in u_borders.iter_mut()
            .enumerate()
            .take(volumes)
            .skip(1)
        {
            *ub -= dt * delta_p[i-1] * s / (
                q[i] + q_interface[i]
            );
        }
        if let Some(last) = u_borders.last_mut() {
            *last += p.last().unwrap() * s * dt / (
                0.5*q.last().unwrap() + task.piston_mass
            );
        }

        mesh = mesh + &u_borders * dt;

        let ubl = u_borders.slice(
            s![..u_borders.len()-1]
        );
        let ubr = u_borders.slice(s![1..]);
        let delta_ub = &ubr - &ubl;
        let e_num =
            2.*(&inner_e * &q) - dt*s*(&p * &delta_ub);
        let xl = mesh.slice(
            s![..mesh.len()-1]
        );
        let xr = mesh.slice(s![1..]);
        let delta_x = &xr - &xl;
        for (i, r) in rho.iter_mut().enumerate() {
            *r = q[i] / (s * delta_x[i]);
        }
        let e_denom =
            2. * &q + (k - 1.) * s * dt * (&rho * &delta_ub);
        for (ei, e1, e2)
        in izip!(inner_e.iter_mut(), e_num, e_denom)
        {
            *ei = e1 / e2;
        }
        for (pi, ei, rhoi)
        in izip!(p.iter_mut(), &inner_e, &rho)
        {
            *pi = (k - 1.)*ei*rhoi;
        }

        u = 0.5*(
            &u_borders.slice(s![..u_borders.len()-1])
            + &u_borders.slice(s![1..])
        );
        t += dt;

        t_store.push(t);
        p_piston_store.push(*p.last().unwrap());
        p_bottom_store.push(*p.first().unwrap());
        piston_x_store.push(*mesh.last().unwrap());
        piston_v_store.push(*u_borders.last().unwrap());
    }

    Results::new(
        t_store,
        p_piston_store,
        p_bottom_store,
        piston_v_store,
        piston_x_store
    )
}


fn calc_time_step(mesh: &Array1<f64>,
                  u: &Array1<f64>,
                  c: &Array1<f64>,
                  cfl: f64) -> f64
{
    let xl = mesh.slice(s![..-1]);
    let xr = mesh.slice(s![1..]);
    let abs_u = u.mapv(f64::abs);
    let dt = (&xr - &xl) / (c + &abs_u);
    cfl * min(&dt.to_vec())
}
