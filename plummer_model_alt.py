import numpy as np

G = 6.67e-11

class Body:
    def __init__(self, idno, pos = np.array([0.0,0.0,0.0]), vel=np.array([0.0,0.0,0.0]), mass=0.0):
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.id = idno
    
    def E_k(self):
        return 0.5 * self.mass * np.linalg.norm(self.vel)**2
    
    def E_p(self, obj_array):
        p = 0
        for x in obj_array.body:
            if x.id != self.id:
                dist = np.linalg.norm(x.pos - self.pos)
                p += (-1 * G * self.mass * x.mass)/dist
        return p
    
class NBody:
    def __init__(self, n):
        self.body = []
        for i in range(n):
            x = Body(idno = i)
            self.body.append(x)
            
    def E_k(self):
        e = 0
        for x in self.body:
            e += x.E_k()
        return e
    
    def E_p(self):
        p = 0
        for x in self.body:
            p += x.E_p(self)
        return p/2
    
    def energy_vals(self):
        e_k = self.E_k()
        e_p = self.E_p()
        e_t = e_k + e_p
        return e_k, e_p, e_t
        
        
def ret_sph(r):
    v = np.array([0.0,0.0,0.0])
    theta = np.arccos(np.random.uniform(-1, 1))
    phi = np.random.uniform(0, 2*np.pi)
    v[0] = r * np.sin(theta) * np.cos(phi)
    v[1] = r * np.sin(theta) * np.sin(phi)
    v[2] = r * np.cos(theta)
    return v

"""
n = number of particles
M = total mass of system
a = structural length scale
"""
     
def make_plummer(n, M, a):
    nb = NBody(n)
    for b in nb.body:
        b.mass = M/n
        radius = a/np.sqrt((M/np.random.uniform(0.1, M))**(2.0/3.0) - 1.0)
        b.pos = ret_sph(radius)
        x = 0.0
        y = 0.1
        while y>(x**2)*((1 - x**2)**3.5):
            x = np.random.uniform(0, 1)
            y = np.random.uniform(0, 0.1)
        velocity = x * np.sqrt(2.0*G*M) * (a**2 + radius**2)**(-0.25)
        b.vel = ret_sph(velocity)
    return nb

