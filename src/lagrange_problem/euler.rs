use std::f64::consts::PI;
use std::iter::zip;

use itertools::{izip, Itertools};

use super::utils::{Results, min, linspace};


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
    cfl: f64,
    alpha: f64,
    beta: f64
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
           cfl: f64,
           alpha: f64,
           beta: f64) -> Task
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
            cfl,
            alpha,
            beta
        }
    }
}


pub fn solve(task: &Task) -> Results
{
    let k = task.gas_k;

    let (mut mesh, mut dx) = linspace(
        0., task.x0, task.n_volumes + 1
    );
    let (mut q, p) = init(
        &mesh, task.p0, task.temper0, task.gas_const, k
    );
    let (mut piston_x, mut piston_v) = (*mesh.last().unwrap(), 0.);
    let s = 0.25 * PI * task.tube_diameter * task.tube_diameter;
    let mut t = 0.;

    let capacity = 1000;
    let mut t_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut p_piston_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut p_bottom_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut piston_v_store: Vec<f64> = Vec::with_capacity(capacity);
    let mut piston_x_store: Vec<f64> = Vec::with_capacity(capacity);

    t_store.push(t);
    p_piston_store.push(*p.last().unwrap());
    p_bottom_store.push(*p.first().unwrap());
    piston_v_store.push(piston_v);
    piston_x_store.push(piston_x);

    let nodes = mesh.len();
    while *mesh.last().unwrap() < task.tube_len {
        let p = equation_of_state(&q, task.gas_k);
        let c = calc_sonic(&p, &q[0], k);
        let dt = calc_time_step(&q, &c, dx, task.cfl);

        piston_v += dt * p[p.len() - 2] * s / task.piston_mass;
        piston_x += dt * piston_v;

        let (new_mesh, new_dx) = linspace(
            0., piston_x, nodes
        );
        let (u_borders, _) = linspace(
            0., piston_v, nodes
        );

        update_boundaries(&mut q, piston_v);

        let flux = ausm_plus(
            &q, &u_borders, &c, &p, task.alpha, task.beta
        );
        let nf = flux[0].len();
        let f_left = [
            &flux[0][..nf], &flux[1][..nf], &flux[2][..nf]
        ];
        let f_right = [
            &flux[0][1..], &flux[1][1..], &flux[2][1..]
        ];
        q = step_up_euler(
            new_dx, dx, &q, &f_left, &f_right, dt
        );

        t += dt;
        (mesh, dx) = (new_mesh, new_dx);

        t_store.push(t);
        p_piston_store.push(p[p.len() - 2]);
        p_bottom_store.push(p[1]);
        piston_v_store.push(piston_v);
        piston_x_store.push(piston_x);
    }

    Results::new(
        t_store, p_piston_store, p_bottom_store, piston_v_store, piston_x_store
    )
}


fn init(mesh: &[f64],
    p0: f64,
    temp0: f64,
    gas_const: f64,
    k: f64) -> (Vec<Vec<f64>>, Vec<f64>)
{
let n = mesh.len() + 1;
let p = vec![p0; n];
let rho = p
    .iter()
    .map(|pi| pi / (temp0 * gas_const))
    .collect_vec();
let e = zip(&p, &rho)
    .map(|(pi, rhoi)| pi / ((k - 1.) * rhoi));

let mut q = vec![vec![0.0; n]; 3];
for (q1, rhoi) in zip(q[0].iter_mut(), &rho) {
    *q1 = *rhoi;
}
for (q3, rhoi, ei)
    in izip!(q[2].iter_mut(), &rho, e) {
    *q3 = rhoi * ei;
}

(q, p)
}


fn update_boundaries(q: &mut [Vec<f64>], piston_v: f64)
{
    q[0][0] = q[0][1];
    q[1][0] = -q[1][1];
    q[2][0] = q[2][1];

    let n = q[0].len();
    q[0][n - 1] = q[0][n - 2];
    q[1][n - 1] = -q[1][n - 2] + 2.*q[0][n - 2]*piston_v;
    q[2][n - 1] = q[2][n - 2];
}


fn equation_of_state(q: &[Vec<f64>], k: f64) -> Vec<f64>
{
    let u = zip(&q[0], &q[1])
        .map(|(rho, rho_u)| rho_u / rho);
    let e = izip!(&q[0], &q[2], u)
        .map(
            |(rho, rho_e, ui)|
            rho_e / rho - 0.5 * ui * ui
        );
    zip(&q[0], e)
        .map(|(rho, ei)| (k - 1.) * ei * rho)
        .collect_vec()
}


fn calc_sonic(p: &Vec<f64>, rho: &Vec<f64>, k: f64) -> Vec<f64>
{
    izip!(p, rho)
        .map(|(pi, rhoi)| (k * pi / rhoi).sqrt())
        .collect_vec()
}


fn ausm_plus(q: &[Vec<f64>],
             u_borders: &[f64],
             c: &[f64],
             p: &[f64],
             alpha: f64,
             beta: f64) -> Vec<Vec<f64>>
{
    let n = q[0].len() - 1;

    let interface_c =
        zip(&c[..n], &c[1..])
        .map(|(cl, cr)| 0.5*(cl + cr))
        .collect_vec();
    let u = zip(&q[0], &q[1]).map(
        |(rho, rho_u)| rho_u / rho
    ).collect_vec();

    let mach_left =
        izip!(&u[..n], u_borders, &interface_c)
        .map(
            |(ui, u_bi, c_ii)|
            (ui - u_bi) / c_ii
        )
        .collect_vec();
    let mach_right =
        izip!(&u[1..], u_borders, &interface_c)
        .map(
            |(ui, u_bi, c_ii)|
            (ui - u_bi) / c_ii
        )
        .collect_vec();

    let q1_left = &q[0][..n];
    let q2_left = &q[1][..n];
    let q3_left = &q[2][..n];
    let p_left = &p[..n];
    let flux_left = calc_ausm_flux(
        &[q1_left, q2_left, q3_left], p_left
    );

    let q1_right = &q[0][1..];
    let q2_right = &q[1][1..];
    let q3_right = &q[2][1..];
    let p_right = &p[1..];
    let flux_right = calc_ausm_flux(
        &[q1_right, q2_right, q3_right], p_right
    );

    let interface_mach = calc_interface_mach(
        &mach_left, &mach_right, beta
    );
    let interface_p = calc_interface_pressure(
        &mach_left, &mach_right, p, alpha
    );

    let m = interface_c.len();
    let mut flux = vec![vec![0.0; m]; 3];
    for i in 0..m {
        flux[0][i] = 0.5*interface_c[i]*(
            interface_mach[i]*(flux_right[0][i] + flux_left[0][i])
            - interface_mach[i].abs()*(flux_right[0][i] - flux_left[0][i])
        );
        flux[1][i] = interface_p[i] + 0.5*interface_c[i]*(
            interface_mach[i]*(flux_right[1][i] + flux_left[1][i])
            - interface_mach[i].abs()*(flux_right[1][i] - flux_left[1][i])
        );
        flux[2][i] = interface_p[i]*u_borders[i] + 0.5*interface_c[i]*(
            interface_mach[i]*(flux_right[2][i] + flux_left[2][i])
            - interface_mach[i].abs()*(flux_right[2][i] - flux_left[2][i])
        );
    }
    
    flux
}


fn calc_ausm_flux(q: &[&[f64]], p: &[f64]) -> Vec<Vec<f64>>
{
    let mut flux = vec![vec![0.0; p.len()]; 3];
    
    for (f1, rho) in flux[0].iter_mut().zip(q[0]) {
        *f1 = *rho;
    }
    for (f2, rho_u) in flux[1].iter_mut().zip(q[1]) {
        *f2 = *rho_u;
    }
    for (f3, rho_e, pi)
        in izip!(flux[2].iter_mut(), q[2], p) {
        *f3 = rho_e + pi;
    }

    flux
}


fn calc_interface_mach(mach_left: &[f64],
                       mach_right: &[f64],
                       beta: f64) -> Vec<f64>
{
    zip(
        f_beta(mach_left, '+', beta),
        f_beta(mach_right, '-', beta)
    )
    .map(|(fl, fr)| fl + fr)
    .collect_vec()
}


fn f_beta(mach: &[f64], sign: char, beta: f64) -> Vec<f64>
{
    if sign == '+' {
        return mach
        .iter()
        .map(
            |m| {
                let abs_m = m.abs();
                if abs_m >= 1. {
                    0.5*(m + abs_m)
                } else {
                    0.25*(m + 1.).powi(2) * (1. + 4.*beta*(m - 1.).powi(2))
                }
            }
        )
        .collect_vec();
    }

    mach
    .iter()
    .map(
        |m| {
            let abs_m = m.abs();
            if abs_m >= 1. {
                0.5*(m - abs_m)
            } else {
                -0.25*((m - 1.).powi(2) * (1. + 4.*beta*(m + 1.).powi(2)))
            }
        }
    )
    .collect_vec()
}


fn calc_interface_pressure(mach_left: &[f64],
                           mach_right: &[f64],
                           p: &[f64],
                           alpha: f64) -> Vec<f64>
{
    let gl = g_alpha(mach_left, '+', alpha);
    let gr = g_alpha(mach_right, '-', alpha);
    let pl = &p[..p.len() - 1];
    let pr = &p[1..];
    
    Vec::from_iter(
        izip!(gl, pl, gr, pr)
        .map(
            |(gl_, pl_, gr_, pr_)|
            gl_*pl_ + gr_*pr_
        )
    )
}


fn g_alpha(mach: &[f64], sign: char, alpha: f64) -> Vec<f64>
{
    if sign == '+' {
        return mach
        .iter()
        .map(
            |m| {
                let abs_m = m.abs();
                if abs_m >= 1. {
                    0.5*(m + abs_m) / m
                } else {
                    (m + 1.).powi(2) * (
                        0.25*(2. - m) + alpha*m*(m - 1.).powi(2)
                    )
                }
            }
        )
        .collect_vec();
    }
    
    mach
    .iter()
    .map(
        |m| {
            let abs_m = m.abs();
            if abs_m >= 1. {
                0.5*(m - abs_m) / m
            } else {
                (m - 1.).powi(2) * (
                    0.25*(2. + m) - alpha*m*(m + 1.).powi(2)
                )
            }
        }
    )
    .collect_vec()
}


fn calc_time_step(q: &[Vec<f64>],
                  c: &[f64],
                  dx: f64,
                  cfl: f64) -> f64
{
    let mut u = vec![0.0; q[0].len()];
    for (ui, ci, mi, pi)
        in izip!(u.iter_mut(), c, &q[0], &q[1])
    {
        *ui = dx / (ci + pi/mi);
    }
    cfl * min(&u)
}


fn step_up_euler(new_dx: f64,
                 dx: f64,
                 q: &[Vec<f64>],
                 f_left: &[&[f64]],
                 f_right: &[&[f64]],
                 dt: f64) -> Vec<Vec<f64>>
{
    let n = q[0].len();
    let mut new_q = q.to_owned();
    
    for i in 1..n - 1 {
        new_q[0][i] = dx/new_dx * (
            q[0][i] - dt/dx * (f_right[0][i-1] - f_left[0][i-1])
        );
        new_q[1][i] = dx/new_dx * (
            q[1][i] - dt/dx * (f_right[1][i-1] - f_left[1][i-1])
        );
        new_q[2][i] = dx/new_dx * (
            q[2][i] - dt/dx * (f_right[2][i-1] - f_left[2][i-1])
        );
    }
    
    new_q
}
