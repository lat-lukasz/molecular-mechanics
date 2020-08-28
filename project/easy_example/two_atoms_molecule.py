import numpy as np
import argparse, sys
import matplotlib.pyplot as plt

from integrator1 import SimpleIntegrator, VerletIntegrator
parser = argparse.ArgumentParser(description='Dynamika molekularna jednowymiarowa cząsteczki dwuatomowej')

parser.add_argument('-n', '--n', type=int, required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float, required=True, help='krok czasowy')
parser.add_argument('-k', '--k', type=float, required=True, help='stała siłowa')
parser.add_argument('-r0', '--r0', type=float, required=True, help='położenie równowagowe')
parser.add_argument('-m', '--mass', type=float,  required=True, help='masa atomu')
parser.add_argument('-x', '--x', type=float, required=True, help='połozenie poczatkowe atomu')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

k = args.k
r0 = args.r0
x = args.x
m = args.mass

integrator = VerletIntegrator(dt, m, x)

trajectory = []
velocity = []
energy = []

for _ in range(n):
    r = x
    F = -k*(r-r0)
    integrator.next(F)
    x = integrator.x
    trajectory.append(x)
    energy.append((k/2)*(r-r0)**2)
    velocity.append(integrator.v)

plt.plot(range(n), trajectory, range(n), energy)
plt.legend(['Połozenie', 'Energia'])
plt.show()
