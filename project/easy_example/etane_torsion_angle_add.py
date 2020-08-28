
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Dynamika molekularna rotacji cząsteczki etanu')

parser.add_argument('-k0', '--k0', type=float, required=True, help='równowagowy kąt torsyjny')
parser.add_argument('-A', '--A',type=float, required=True, help='wysokość bariery rotacyjnej')
parser.add_argument('-n', '--n',type=float, required=True, help='krotność kąta równowagowego')

args = parser.parse_args(sys.argv[1:])

k0 = args.k0
A = args.A
n = args.n

x = np.linspace(0,360,360)
y = A*(1+np.cos(n*np.deg2rad(x)-np.deg2rad(k0)))

plt.plot(x,y)
plt.show()
