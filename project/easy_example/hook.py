import argparse, sys
import matplotlib.pyplot as plt

from integrator1 import SimpleIntegrator, VerletIntegrator
parser = argparse.ArgumentParser(description='Ciezarek na sprezynie')

parser.add_argument('-n', '--n', type=int,
required=True, help='liczba iteracji')
parser.add_argument('--dt', type=float,
required=True, help='krok czasowy')
parser.add_argument('-k', '--k', type=float,
required=True, help='stała siłowa sprezyny')
parser.add_argument('--x0', type=float,
required=True, help='długosc równowagowa sprezyny')
parser.add_argument('-m', '--mass', type=float,
required=True, help='masa ciezarka')
parser.add_argument('-x', '--x', type=float,
required=True, help='połozenie poczatkowe ciezarka')

args = parser.parse_args(sys.argv[1:])

n = args.n
dt = args.dt

k = args.k
x0 = args.x0
m = args.mass
x = args.x

integrator = SimpleIntegrator(dt, m, x)

trajectory = []
velocity = []
energy = []

for _ in range(n):
    F = -k * (x - x0)
    integrator.next(F)
    x = integrator.x
    trajectory.append(x)
    velocity.append(integrator.v)
    energy.append(m * integrator.v**2 + k / 2 * (x - x0)**2)

plt.plot(range(n), trajectory, range(n), velocity, range(n), energy)
plt.legend(['Połozenie', 'Predkosc', 'Energia'])
plt.show()
