import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la

from integrator_3d import VerletIntegrator

parser = argparse.ArgumentParser(description='Czasteczka dwuatomowa')

parser.add_argument('-n', '--n', type=int, required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float, required=True, help='krok czasowy')
parser.add_argument('-k', '--k', type=float, required=True, help='stała siłowa')
parser.add_argument('--r0', type=float, required=True, help='odległosc równowagowa')
parser.add_argument('-m', '--mass', type=float, nargs=2, required=True, help='masy atomów')
parser.add_argument('-x', '--x', type=float, nargs=6, required=True, help='połozenia poczatkowe atomów')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

k = args.k
r0 = args.r0

m = np.array(args.mass)
x = np.empty([3, 2])
x[:, 0] = np.array(args.x[0:3])
x[:, 1] = np.array(args.x[3:6])

integrator = VerletIntegrator(dt, m, x)

trajectory = []

for _ in range(n):
    u = x[:, 1] - x[:, 0]
    r = la.norm(u)
    F = np.empty([3, 2])
    F[:, 0] = k * (r - r0) * u / r
    F[:, 1] = -k * (r - r0)* u / r
    integrator.next(F)
    x = integrator.x
    trajectory.append(r)

plt.plot(range(n), trajectory)
plt.legend(['Odległosc'])
plt.show()
