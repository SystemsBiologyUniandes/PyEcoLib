import numpy as np

class Bacteria:
    def __init__(self, idx: int, size: float, mom_idx: int, added: float, mu: float):
        self.idx = idx
        self.size = size
        self.base_size = size
        self.mom_idx = mom_idx
        self.added = added
        self.mu = mu
        
    def evolve(self, time: float) -> None:
        self.added += self.size * (np.exp(self.mu * time) - 1)
        self.size = self.size * np.exp(self.mu * time)
        
    def divide(self) -> None:
        self.size /= 2
        self.added = 0
        self.base_size = self.size