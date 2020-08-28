import numpy as np

class Integrator:
    def __init__(self, dt, m, x):
        self.m = m
        self.x = x
        self.dt = dt

class SimpleIntegrator(Integrator):
    def __init__(self, dt, m, x):
        super().__init__(dt, m, x)
        self.v = np.zeros_like(x)

    def next(self, F):
        x = self.x
        a = np.empty_like(F)
        for i in range(F.shape[1]):
            a[:, i] = F[:, i] / self.m[i]
        self.x = self.x + self.v * self.dt + 0.5 * a * self.dt**2
        self.v = self.v + 0.5 * a * self.dt

class VerletIntegrator(Integrator):
    def __init__(self, dt, m, x):
        super().__init__(dt, m, x)
        self._prev_x = x

    def next(self, F):
        x = self.x
        a = np.empty_like(F)
        for i in range(F.shape[1]):
            a[:, i] = F[:, i] / self.m[i]
        self.x = 2 * self.x - self._prev_x + 0.5 * a * self.dt**2
        self.v = (self.x - self._prev_x) / (2 * self.dt)
        self._prev_x = x
