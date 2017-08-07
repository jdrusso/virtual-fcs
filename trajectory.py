#!/usr/bin/env/python
import numpy as np

class Atom:
    
    def __init__(self, index, coords):
        self.xyz = coords
        self.index = index
    
    def setIndex(self, i):
        self.index = i
        
class Topology:
    
    def __init__(self, residues):
        self.n_residues = residues
        
class FakeTrajectory:
    
    def initialize(self, coords, walkers):
        self.topology = Topology(walkers)

        self.topology.atoms = list([Atom(i, coords[i]) for i in range(walkers)])

        self.topology.n_residues = walkers

        self.xyz = coords
        
    def reduce(self, sampled_ratio):
        
        new_tot = int(self.topology.n_residues * sampled_ratio)
        
        self.topology.n_residues = new_tot
        
        kept_atoms = np.random.choice(self.topology.atoms, new_tot, replace=False)
        print("XYZ: %d, N_ATOMS: %d" % (len(self.xyz), len(self.topology.atoms)))
        self.topology.atoms = kept_atoms
        print("XYZ: %d, N_ATOMS: %d" % (len(self.xyz), len(self.topology.atoms)))
        
        self.xyz = self.xyz[:,[atom.index for atom in kept_atoms]]
        
        # Reindex
        for i in range(new_tot):
            self.topology.atoms[i].setIndex(i)
            
    def reindex(self):
        # Reindex
        for i in range(self.xyz.shape[1]):
            self.topology.atoms[i].setIndex(i)
        
    def __len__(self):
        
        return len(self.xyz)