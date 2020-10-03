from random import random
import numpy as np

class Cell:
    def __init__(self, idx, V0, num_steps, gr, divpar, k):
        
        
        self.dp = divpar
        self.gr = gr
        self.idx = idx
        self.V = V0
        self.Vb = V0
        self.Vd = V0
        self.total_steps = num_steps
        self.num_steps = 0
        self.rv = np.random.rand()
        self.ndiv = 0
        self.nextt = 100/gr
        self.k = k

    def change(self, Vn):
        self.V = Vn
        self.rv = np.random.rand()


    def division(self, Vn, dp, gr, k):
        self.gr = gr
        self.k = k
        self.Vb = self.Vd*self.dp
        self.Vd = Vn
        self.dp = dp
        self.V = Vn*self.dp
        self.rv = np.random.rand()
        self.ndiv += 1
        self.num_steps = 0


    def add_growth(self, Vn):
        self.V = Vn

    def get_size(self):
        return self.V
    
    def get_divpar(self):
        return self.dp

    def __str__(self):
        return "Cell: {\n   idx: "+str(self.idx)+"\n    V0: "+str(self.V)+"\n}"