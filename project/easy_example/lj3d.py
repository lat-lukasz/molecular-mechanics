import argparse
import math
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
from integrator_3d import VerletIntegrator

parser = argparse.ArgumentParser(
description='Dynamika molekularna w potencjale LJ')

parser.add_argument('-n', '--n', type=int,required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float, required=True, help='krok czasowy')
parser.add_argument('-A', '--A', type=float,required=True, help='parametr A')
parser.add_argument('-B', '--B', type=float, required=True, help='parametr B')
parser.add_argument('-m', '--mass', type=float, nargs=2,required=True, help='masy atomów')
parser.add_argument('-x', '--x', type=float, nargs=6,required=True, help='połozenia poczatkowe atomów')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

A = args.A
B = args.B

m = np.array(args.mass)
x = np.empty([3, 2])
x[:, 0] = np.array(args.x[0:3])
x[:, 1] = np.array(args.x[3:6])

integrator = VerletIntegrator(dt, m, x)

trajectory = []

for _ in range(n):
    u = x[:, 1] - x[:, 0]
    r2 = np.dot(u, u)
    r6 = r2**3
    r8 = r6 * r2
    r14 = r8 * r6
    F = np.empty([3, 2])
    F[:, 0] = -(12 * A / r14 - 6 * B / r8) * u
    F[:, 1] = (12 * A / r14 - 6 * B / r8) * u
    integrator.next(F)
    x = integrator.x
    
    trajectory.append(math.sqrt(r2))

plt.plot(range(n), trajectory)
plt.legend(['Odległosc'])
plt.show()
