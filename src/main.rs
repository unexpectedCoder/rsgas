use std::fs::create_dir_all;
use std::time::Instant;
use std::path::Path;

use crate::lagrange_problem::{
    accurate::{
        solve as solve_accurately,
        Task as AccTask
    },
    euler::{
        solve as solve_euler,
        Task as EulerTask
    },
    lagrange::{
        solve as solve_lagrange,
        Task as LagrangeTask
    },
    lagrange_nd::solve as solve_lagrange_nd
};

pub mod lagrange_problem;


fn main()
{
    let piston_mass = 0.1;
    let tube_len = 2.;
    let tube_diameter = 0.03;
    let gas_const = 287.;
    let gas_k = 1.4;
    let p0 = 5e6;
    let temper0 = 300.;
    let x0 = 0.5;
    let n_volumes = 300;
    let cfl = 0.8;
    let alpha = 0.25;
    let beta = 0.1875;

    let acc_task = AccTask::new(
        piston_mass,
        tube_len,
        tube_diameter,
        gas_const,
        gas_k,
        p0,
        temper0,
        x0
    );
    let euler_task = EulerTask::new(
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
    );
    let lagrange_task = LagrangeTask::new(
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
    );

    let mut timer = Instant::now();
    let acc_sol = solve_accurately(&acc_task, 200);
    let mut dt = timer.elapsed();
    println!("Accurate time is {:.2?}", dt);
    timer = Instant::now();
    let euler_sol = solve_euler(&euler_task);
    dt = timer.elapsed();
    println!("Euler time is {:.2?}", dt);
    timer = Instant::now();
    let lagrange_sol = solve_lagrange(&lagrange_task);
    dt = timer.elapsed();
    println!("Lagrange time is {:.2?}", dt);
    timer = Instant::now();
    let lagrange_ndsol = solve_lagrange_nd(&lagrange_task);
    dt = timer.elapsed();
    println!("LagrangeND time is {:.2?}", dt);

    let res_dir = Path::new("results");
    if !res_dir.try_exists().unwrap() {
        if let Err(why) = create_dir_all(res_dir) {
            panic!(
                "error while creating dirs ({})", why
            )
        };
    }

    acc_sol.save_csv(res_dir.join("acc_results.csv").as_path());
    euler_sol.save_csv(res_dir.join("euler_results.csv").as_path());
    lagrange_sol.save_csv(res_dir.join("lagrange_results.csv").as_path());
    lagrange_ndsol.save_csv(
        res_dir.join("lagrange_results_nd.csv").as_path()
    );
}
