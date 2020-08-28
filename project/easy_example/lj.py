import numpy as np
import argparse, sys
import matplotlib.pyplot as plt

from integrator1 import SimpleIntegrator, VerletIntegrator
parser = argparse.ArgumentParser(description='Dynamika molekularna jednowymiarowa dwóch atomów helu')

parser.add_argument('-n', '--n', type=int, required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float, required=True, help='krok czasowy')
parser.add_argument('-A', '--A', type=float, required=True, help='parametr A')
parser.add_argument('-B', '--B', type=float, required=True, help='parametr B')
parser.add_argument('-m', '--mass', type=float,  required=True, help='masa atomu helu')
parser.add_argument('-x', '--x', type=float, required=True, help='połozenie poczatkowe atomow')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

A = args.A
B = args.B
x = args.x
m = args.mass

integrator = VerletIntegrator(dt, m, x)

trajectory = []
velocity = []
energy = []

for _ in range(n):
    r = x
    r7 = r**7
    r13 = r**13
    F = 12*A/r13-6*B/r7
    integrator.next(F)
    x = integrator.x
    trajectory.append(x)
    energy.append(B/(r**12) - A/(r**6))
    velocity.append(integrator.v)

plt.plot(range(n), trajectory, range(n), energy)
plt.legend(['Połozenie', 'Energia'])
plt.show()
