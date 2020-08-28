class Integrator:
    def __init__(self, dt, m, x):
        self.m = m
        self.x = x
        self.dt = dt
        
class SimpleIntegrator(Integrator):
    def __init__(self, dt, m, x):
        super().__init__(dt, m, x)
        self.v = 0

    def next(self, F):
        x = self.x
        self.x = self.x + self.v * self.dt + F / (2 * self.m) * self.dt**2     
        self.v = self.v + F / (2 * self.m) * self.dt                           
        
class VerletIntegrator(Integrator):
    def __init__(self, dt, m, x):
        super().__init__(dt, m, x)
        self._prev_x = x

    def next(self, F):
        x = self.x
        self.x = 2 * self.x - self._prev_x + F / (2 * self.m) * self.dt**2
        self.v = (self.x - self._prev_x) / (2 * self.dt)
        self._prev_x = x
