use std::f64::consts::PI;

use itertools::{izip, zip_eq, Itertools};

use super::utils::{Results, linspace, min};


pub struct Task {
    piston_mass: f64,
    tube_len: f64,
    tube_diameter: f64,
    gas_const: f64,
    gas_k: f64,
    p0: f64,
    temper0: f64,
    x0: f64,
    n_volumes: usize,
    cfl: f64
}


impl Task {
    pub fn new(piston_mass: f64,
           tube_len: f64,
           tube_diameter: f64,
           gas_const: f64,
           gas_k: f64,
           p0: f64,
           temper0: f64,
           x0: f64,
           n_volumes: usize,
           cfl: f64) -> Task
    {
        Task {
            piston_mass,
            tube_len,
            tube_diameter,
            gas_const,
            gas_k,
            p0,
            temper0,
            x0,
            n_volumes,
            cfl
        }
    }
}


pub fn solve(task: &Task) -> Results
{
    let k = task.gas_k;
    let volumes = task.n_volumes;
    let nodes = task.n_volumes + 1;

    let (mut mesh, dx): (Vec<f64>, f64) = linspace(
        0., task.x0, nodes
    );
    let mut u_borders = vec![0.0; nodes];
    let mut p = vec![task.p0; volumes];
    let rho0 = task.p0 / (task.gas_const * task.temper0);
    let mut rho = vec![rho0; volumes];
    let mut inner_e =
        zip_eq(&p, &rho)
        .map(|(pi, rhoi)| pi / ((k - 1.) * rhoi))
        .collect_vec();
    
    let s = 0.25 * PI * task.tube_diameter * task.tube_diameter;
    let q = vec![rho0 * s * dx; volumes];
    let mut q_interface = vec![0.0; nodes];
    if let Some(last) = q_interface.last_mut() {
        *last = task.piston_mass;
    }

    let mut t = 0.;
    let mut u = vec![0.; volumes];

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
        let c = calc_sonic(&p, &rho, k);
        let dt = calc_time_step(&mesh, &u, &c, task.cfl);
        let delta_p =
            zip_eq(&p[..p.len() - 1], &p[1..])
            .map(|(pl, pr)| pr - pl)
            .collect_vec();
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

        for (mx, ub) in mesh.iter_mut().zip(&u_borders) {
            *mx += dt * ub;
        }
        
        let delta_u =
            zip_eq(&u_borders[..volumes], &u_borders[1..])
            .map(|(ul, ur)| ur - ul)
            .collect_vec();
        let e_num =
            izip!(&inner_e, &q, &p, &delta_u)
            .map(
                |(ei, qi, pi, dui)|
                2.*ei*qi - dt*s*pi*dui
            )
            .collect_vec();
        let delta_x =
            zip_eq(&mesh[..volumes], &mesh[1..])
            .map(|(xl, xr)| xr - xl)
            .collect_vec();
        for (i, r) in rho.iter_mut().enumerate() {
            *r = q[i] / (s * delta_x[i]);
        }
        let e_denom =
            izip!(&q, &rho, &delta_u)
            .map(
                |(qi, ri, dui)|
                2.*qi + (k - 1.)*ri*s*dui*dt 
            )
            .collect_vec();
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

        u = zip_eq(&u_borders[..volumes], &u_borders[1..])
            .map(|(ul, ur)| 0.5*(ul + ur))
            .collect_vec();
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


fn calc_sonic(p: &[f64], rho: &[f64], k: f64) -> Vec<f64>
{
    zip_eq(p, rho).map(
        |(pi, rhoi)| (k * pi / rhoi).sqrt()
    ).collect_vec()
}


fn calc_time_step(mesh: &[f64], u: &[f64], c: &[f64], cfl: f64) -> f64
{
    cfl * min(
        &izip!(
            &mesh[..mesh.len() - 1], &mesh[1..], u, c
        )
        .map(
            |(xl, xr, ui, ci)|
            (xr - xl) / (ci + ui.abs())
        )
        .collect_vec()
    )
}
