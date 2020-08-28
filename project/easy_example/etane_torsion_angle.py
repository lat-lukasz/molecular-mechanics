import numpy as np
import argparse, sys
import matplotlib.pyplot as plt

from integrator1 import SimpleIntegrator, VerletIntegrator
parser = argparse.ArgumentParser(description='Dynamika molekularna jednowymiarowa rotacji cząsteczki')

parser.add_argument('-n', '--n', type=int, required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float, required=True, help='krok czasowy')
parser.add_argument('-k0', '--k0', type=float, required=True, help='równowagowy kąt torsyjny')
parser.add_argument('-A', '--A',type=float, required=True, help='wysokość bariery rotacyjnej')
parser.add_argument('-m', '--m',type=float, required=True, help='krotność kąta równowagowego')
parser.add_argument('-x', '--x',type=float, required=True, help='kąt początkowy')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

k0 = np.deg2rad(args.k0)
A = args.A
m = args.m
x = np.deg2rad(args.x)

integrator = VerletIntegrator(dt, m, x)

trajectory = []
energy = []

for _ in range(n):
    F = -m*A*np.sin(k0-m*x)
    integrator.next(F)
    x = integrator.x
    trajectory.append(np.rad2deg(x))
    energy.append(A*(1+np.cos(m*x-k0)))

plt.plot(range(n), trajectory, range(n), energy)
plt.legend(['Połozenie', 'Energia'])
plt.show()
