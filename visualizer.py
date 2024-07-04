#!.venv/bin/python

import csv
import matplotlib.pyplot as plt
import os


RES_DIR = "results"
PICS_DIR = os.path.join("results", "pics")
if not os.path.isdir(PICS_DIR):
    os.makedirs(PICS_DIR)

with open(os.path.join(RES_DIR, "acc_results.csv"), "r") as f:
    data = csv.DictReader(f)
    acc_sol = {k: [] for k in data.fieldnames}
    for line in data:
        for k, v in line.items():
            acc_sol[k].append(float(v))

with open(os.path.join(RES_DIR, "euler_results.csv"), "r") as f:
    data = csv.DictReader(f)
    euler_sol = {k: [] for k in data.fieldnames}
    for line in data:
        for k, v in line.items():
            euler_sol[k].append(float(v))

with open(os.path.join(RES_DIR, "lagrange_results.csv"), "r") as f:
    data = csv.DictReader(f)
    lagrange_sol = {k: [] for k in data.fieldnames}
    for line in data:
        for k, v in line.items():
            lagrange_sol[k].append(float(v))


with open(os.path.join(RES_DIR, "lagrange_results_nd.csv"), "r") as f:
    data = csv.DictReader(f)
    lagrange_sol = {k: [] for k in data.fieldnames}
    for line in data:
        for k, v in line.items():
            lagrange_sol[k].append(float(v))

plt.style.use("bmh")

# Euler
fig, ax = plt.subplots(num="Euler x")
ax.plot(euler_sol["t"], euler_sol["x"], label="Численное решение")
ax.plot(acc_sol["t"], acc_sol["x"], ls="--", label="Точное решение")
ax.set(xlabel="$t$, с", ylabel="$x$, м")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Euler x.png"), dpi=300)

fig, ax = plt.subplots(num="Euler v")
ax.plot(euler_sol["t"], euler_sol["v"], label="Численное решение")
ax.plot(acc_sol["t"], acc_sol["v"], ls="--", label="Точное решение")
ax.set(xlabel="$t$, с", ylabel="$v$, м/с")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Euler v.png"), dpi=300)

fig, ax = plt.subplots(num="Euler p")
ax.plot(euler_sol["t"], euler_sol["p_piston"], label="На поршень")
ax.plot(euler_sol["t"], euler_sol["p_bottom"], label="На дно трубы")
ax.set(xlabel="$t$, с", ylabel="$p$, Па")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Euler p.png"), dpi=300)

# Lagrange
fig, ax = plt.subplots(num="Lagrange x")
ax.plot(lagrange_sol["t"], lagrange_sol["x"], label="Численное решение")
ax.plot(acc_sol["t"], acc_sol["x"], ls="--", label="Точное решение")
ax.set(xlabel="$t$, с", ylabel="$x$, м")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Lagrange x.png"), dpi=300)

fig, ax = plt.subplots(num="Lagrange v")
ax.plot(lagrange_sol["t"], lagrange_sol["v"], label="Численное решение")
ax.plot(acc_sol["t"], acc_sol["v"], ls="--", label="Точное решение")
ax.set(xlabel="$t$, с", ylabel="$v$, м/с")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Lagrange v.png"), dpi=300)

fig, ax = plt.subplots(num="Lagrange p")
ax.plot(lagrange_sol["t"], lagrange_sol["p_piston"], label="На поршень")
ax.plot(lagrange_sol["t"], lagrange_sol["p_bottom"], label="На дно трубы")
ax.set(xlabel="$t$, с", ylabel="$p$, Па")
ax.legend()
fig.savefig(os.path.join(PICS_DIR, "Lagrange p.png"), dpi=300)

plt.show()
