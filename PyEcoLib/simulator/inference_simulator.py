import json
import numpy as np
from tqdm import tqdm

from ..models.bacteria import Bacteria

# mu0: np.log(2)
# cv2: 0.1
# delta_time: 0.1
# division_time: 0.001
# noise: 5.597947996804339
# base_size: 1

class InferenceSimulator:
    def __init__(self, max_time: float, mu0: float, cv2: float, delta_time: float, division_time: float, noise: float, base_size: float):
        self.mu0 = mu0
        self.added = 0
        self.cv2 = cv2
        
        self.max_time = max_time 
        self.delta_time = delta_time
        self.dt = division_time
        self.data_pop = []
        self.cv2s = self.cv2 / 3
        self.noise = noise
        self.base_size = base_size
        self.simulation_data = []
        
    def simulate(self) -> json:
        for cell in tqdm(range(100)):
            time = 0
            size = np.random.gamma(shape=1/self.cv2s, scale=1*self.cv2s)
            bacterium = Bacteria(cell, size, -1, 0, self.mu0)
            bacterium.added = np.random.exponential(0.5)
            self.data_pop.append([cell, time, bacterium.size, bacterium.added])
            pop = [bacterium]
            counter = 0
            is_division = False
            while time < self.max_time:
                temp_pop = []
                for i in range(len(pop)):
                    pop[i].evolve(self.dt)
                    divisions = pop[i].added
                    size = pop[i].size
                    mu = self.mu0
                    numerator = np.exp(self.noise * size - np.exp(self.mu0 * self.dt) * self.noise * size) * (1 + np.exp((-divisions + self.base_size) * self.noise))
                    denominator = 1 + np.exp(self.noise * (-divisions + self.base_size + size - np.exp(self.base_size * self.dt) * size))
                    rate = numerator / denominator
                    if np.random.rand() < 1 - rate:
                        self.simulation_data.append({cell, time, pop[i].size, pop[i].added})
                        pop[i].divide()
                        
                counter += 1
                time += self.dt
                if counter >= 100:
                    for cell in pop:
                        self.data_pop.append({cell.idx, time, cell.size, cell.added})
                    counter = 0
        return self.simulation_data