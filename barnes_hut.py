import numpy as np
class cel_obj:
    def __init__(obj, x, y, z, v_x, v_y, v_z, m):
        obj.pos = np.array([x, y, z])
        obj.vel = np.array([v_x, v_y, v_z])
        obj.m = m

class Octtree:
    def __init__(self, center, size, points, masses, ids, leaves = []):
        self.center = center
        self.size = size
        self.children = []
        
        n = len(points)
        
        if n==1:
            leaves.append(self)
            self.COM = points[0]
            self.mass = masses[0]
            self.id = ids[0]
            self.g = np.zeros(3)
        else:
            self.generate_children(points, masses, ids, leaves)
            
            com_total = np.zeros(3)
            m_total = 0.
            for c in self.children:
                m, com = c.mass, c.COM
                m_total += m
                com_total += com * m
            self.mass = m_total
            self.COM = com_total/self.mass
        
    def generate_children(self, points, masses, ids, leaves):
        octant_index = (points > self.center)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    in_octant = np.all(octant_index == np.bool_([i,j,k]), axis=1)
                    if not np.any(in_octant): continue
                    dx = 0.5*self.size*(np.array([i,j,k])-0.5)
                    self.children.append(Octtree(self.center+dx, self.size/2, points[in_octant], masses[in_octant], ids[in_octant], leaves))
        
def TreeWalk(node, node0, thetamax=0.7, G=6.67e-11, soft=0.005):
    dx = node.COM - node0.COM
    r = np.linalg.norm(dx)
    a = node.size/(r + soft)
    if r>0:
        if len(node.children)==0 or a < thetamax:
            node0.g += G * node.mass * dx/r**3
        else:
            for c in node.children: TreeWalk(c, node0, thetamax, G)       

def GravAccel(points, masses, thetamax=0.7, G=6.67e-11):
    center = (np.max(points,axis=0)+np.min(points,axis=0))/2
    topsize = np.max(np.max(points,axis=0)-np.min(points,axis=0))
    leaves = []
    topnode = Octtree(center, topsize, points, masses, np.arange(len(points)), leaves)
    accel = np.zeros_like(points)
    for i,leaf in enumerate(leaves):
        TreeWalk(topnode, leaf, thetamax, G)
        accel[leaf.id] = leaf.g
    return accel