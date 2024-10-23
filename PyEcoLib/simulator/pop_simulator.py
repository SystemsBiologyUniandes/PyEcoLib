import json
import numpy as np

from src.models.cell import Cell
from src.utils.utils import Utils

class PopSimulator:
    def __init__(self, n_cells: int, growth_rate: float, size_base: float, division_steps: int, cv2_div: float = 0, cv2_gr: float = 0, lamb: float = 1, initial_array: list[Cell] = None, nu: int = 1):
        """
        PopSimulator simulate the entire population or the Simulator changing
        one parameter (nu)
        
        Args:
            n_cells (int):  Number of cells in simulation
            growth_rate (float): Rate of growthing of cells in simulation
            size_base (float): Size base of the cells
            division_steps (int): Number of steps before division
            cv2_div (float): Variability size division to division
            cv2_gr (float): Variability growth rate
            lamb (float): Lambda value
            nu (int): Offspring tracking
            initial_array (list): List of initial size of bacteriums
        
        Returns:
            None
        """
        self.__check_errors(n_cells, growth_rate, size_base, cv2_div, cv2_gr, lamb, nu)
        
        self.n_cells = n_cells
        self.growth_rate = growth_rate
        self.size_base = size_base
        self.division_steps = division_steps
        self.lamb = lamb
        self.nu = nu
        self.cv2_div = cv2_div
        self.cv2_gr = cv2_gr
        
        self.utils = Utils(self.growth_rate, self.total_steps, self.size_base, self.lamb)
        
        self.sample_time = 0
        self.cells: list[Cell] = []
        self.initial_array = initial_array
        
        if self.lamb == 1:
            self.k = self.division_steps * self.growth_rate / self.size_base
        else:
            self.k = self.division_steps * self.utils.get_k()
            
        self.initialize_cells()
        
    def size_dynamics(self, max_time: float, sample_time: float) -> json:
        """
        Method to do size dynamics
        
        Args:
            max_time (float): Max time of simulation
            sample_time (float): Sample time of simulation
            
        Returns:
            json: Returns json with simulation data
        """
        self.initialize_cells()
        
        if self.nu == 2:
            if max_time > 9 * np.log(2) / self.growth_rate:
                raise NameError('Please select max time < 9 * doubling time')
            elif self.n_cells >= 2000:
                raise NameError('Please select a smaller number of cells')
        
        self.time = 0
        
        next_array = np.empty(len(self.cells))
        for i in range(len(self.cells)):
            next_array[i] = self.cells[i]
            self.cells[i].idx = i
        
        ref_time = sample_time
        output = []
        while self.time < max_time:
            next_time = np.min(next_array)
            if next_time > ref_time - self.time:
                target_time = ref_time - self.time
                self.time += target_time
                for i in range(len(self.cells)):
                    cell = self.cells[i]
                    cell.next_time -= target_time
                    next_array[i] -= target_time
                    growth_rate = cell.growth_rate * self.growth_rate
                    cell.size *= np.exp(growth_rate * target_time)
                    
                    for cell in self.cells:
                        output.append(Utils.json_pop_output(cell.idx, self.time, Utils.truncate(cell.size, 4), cell.step))
                        
                ref_time += self.sample_time
            else:
                m = np.argmin(next_array)
                cell = self.cells[m]
                cell.step += 1
                for i in range(len(self.cells)):
                    self.cells[i].next_time -= next_time
                    next_array[i] -= next_time
                    growth_rate = self.cells[i].growth_rate * self.growth_rate
                    self.cells[i].size *= np.exp(growth_rate * next_time)
                
                if cell.step >= self.division_steps:
                    mom_size = cell.size
                    size = cell.size * cell.div_par
                    size2 = cell.size * (1 - cell.div_par)
                    cell.size = size
                    cell.step = 0
                    cell.div_par = Utils.new_division_par(self.cv2_div)
                    cell.growth_rate = Utils.new_growth_rate(self.cv2_gr)
                    growth_rate = cell.growth_rate
                    div_par = cell.div_par
                    cell.next_time = self.utils.get_next_time(size, np.random.rand(), cell.growth_rate, cell.k)
                    next_array[m] = cell.next_time
                    if self.nu == 2:
                        growth_rate2 = Utils.new_growth_rate(self.cv2_gr)
                        div_par2 = Utils.new_division_par(self.cv2_div)
                        idx = len(self.cells)
                        cell2 = Cell(idx, initial_size=size2, num_steps=self.division_steps, growth_rate=growth_rate2, div_par=div_par2, k = growth_rate2)
                        cell2.idx = cell.idx
                        cell2.next_time = self.utils.get_next_time(size2, np.random.rand(), cell2.growth_rate, cell2.k)
                        next_array = np.concatenate((next_array, [cell2.next_time]), axis=None)
                        self.cells.append(cell2)
                    
                    output.append({
                        "idx": cell.idx,
                        "mother_idx": cell.idx,
                        "mother_size": mom_size,
                        "time": next_time + self.time,
                        "size": np.round(size, 8),
                        "growth_rate": np.round(growth_rate * self.growth_rate),
                        "div_par": div_par
                    })
                    output.append({
                        "idx": cell2.idx,
                        "mother_idx": cell.idx,
                        "mother_size": mom_size,
                        "time": next_time + self.time,
                        "size": np.round(size2, 8),
                        "growth_rate": np.round(growth_rate2 * self.growth_rate, 8),
                        "div_par": div_par2
                        
                    })
                
                else:
                    cell.random_value = np.random.rand()
                    cell.next_time = self.utils.get_next_time(cell.size, cell.random_value, cell.growth_rate, cell.k)
                    next_array[m] = cell.next_time
                
                self.time += next_time
        return output
    
    def division_strategy(self, max_time: float) -> json:
        """
        Method for simulate division strategy
        
        Args:
            max_time (float): Max time of simulation
            
        Returns:
            json: Returns json with the data of the simulation
        """
        
        if self.nu == 2: raise NameError('This function was designed for nu = 1')
        
        self.initialize_cells()
        self.time = 0
        sample_time: float = 0.01 * max_time
        
        division_array = np.array([])
        output = []
        
        for cell in self.cells:
            division_array = np.concatenate((division_array, [0]), axis=0)
            
        while self.time < max_time:
            growth_rate_array = np.array([])
            for cell in self.cells:
                growth_rate_array = np.concatenate((growth_rate_array, [cell.growth_rate]), axis=0)
            self.simulate(sample_time)
            self.time += sample_time
            counter = 0
            
            for cell in self.cells:
                if cell.num_divisions > division_array[counter]:
                    mu: float = self.growth_rate * growth_rate_array[counter]
                    td: float = (1 / mu) * np.log(cell.division_size / cell.base_size)
                    if self.time > 0.3 * max_time:
                        result = {
                            'time': Utils.truncate(self.time, 4),
                            'size_at_birth': Utils.truncate(cell.base_size, 4),
                            'division_size': Utils.truncate(cell.division_size, 4),
                            'mu': Utils.truncate(mu, 4),
                            'td': Utils.truncate(td, 4)
                        }
                        output.append(result)
                    division_array[counter] = cell.num_divisions
                counter += 1
        return output
                    
            
        
    def simulate(self, max_time: float) -> None:
        """
        Function that simulate cells in a environment during a
        period of time
        
        Args:
            max_time (float): Max time of simulation
        
        Returns:
            None 
        """
        
        if self.nu == 2:
            raise NameError('This function was designed for Nu = 1')
        
        for cell in self.cells:
            time: float = 0
            while time < max_time:
                tt: float = cell.next_time
                if (time + tt) <= max_time:
                    cell.step += 1
                    size: float = cell.size * np.exp(self.growth_rate * cell.growth_rate * tt)
                    if cell.step >= cell.num_steps:
                        div_par = Utils.new_division_par(self.cv2_div)
                        growth_rate = Utils.new_growth_rate(self.cv2_gr)
                        cell.division(size, div_par, growth_rate, growth_rate)
                    else:
                        cell.change(size)
                    cell.next_time = self.utils.get_next_time(cell.size, cell.random_value, cell.growth_rate, cell.k)
                else:
                    size: float = cell.size * np.exp(self.growth_rate * cell.growth_rate * (max_time - time))
                    cell.change(size)
                    cell.next_time = cell.next_time - (max_time - time)
                
                time += tt
        
    def initialize_cells(self) -> str:
        """
        Function that initialize cells for simulation
        
        Args:
            None
        Returns:
            str: Returns a json with the data of the cells
        """
        
        self.cells = []
        if len(self.initial_array) > 0:
            idx: int = 0
            for size in self.initial_array:
                growth_rate: float = Utils.new_growth_rate(self.cv2_gr)
                division_par: float = Utils.new_division_par(self.cv2_div)
                cell: Cell = Cell(idx, size, num_steps=self.total_steps, growth_rate=growth_rate, div_par=division_par, k=growth_rate)
                cell.next_time = self.utils.get_next_time(size, cell.random_value, cell.growth_rate, cell.k)
                self.cells.append(cell)
                idx += 1
            
            json_cells = [cell.json() for cell in self.cells]
            return json.dumps(json_cells)
        
        for i in range(self.n_cells):
            growth_rate: float = Utils.new_growth_rate(self.cv2_gr)
            div_par: float = Utils.new_division_par(self.cv2_div)
            cell: Cell = Cell(i, self.size_base, num_steps=self.total_steps, growth_rate=growth_rate, div_par=div_par, k=growth_rate)
            cell.next_time = self.utils.get_next_time(self.size_base, cell.random_value, cell.growth_rate, cell.k)
            self.cells.append(cell)
        
        json_cells = [cell.json() for cell in self.cells]
        return json.dumps(json_cells)
    
    def __check_errors(self, n_cells: int, growth_rate: float, size_base: float, division_steps: int, cv2_div: float, cv2_gr: float, lamb: float, nu: float) -> None:
        """
        Validate input parameters
        
        Args:
            n_cells (int): Number of cells in the simulation
            growth_rate (float): Rate of growthing of cells in simulation
            size_base (float): Size base of the cells
            division_steps (int): Number of steps before division
            cv2_div (float): Target division size must be positive or zero
            cv2_gr (float): Target growth rate must be positive or zero
            lamb (float): Lambda value
            nu (float): Offspring tracking
        
        Returns:
            None
        """
        
        if n_cells <= 0: raise NameError('The number of cells must be positive')
        if growth_rate < 0: raise NameError('The growth rate must be positive')
        if size_base < 0: raise NameError('The size base of the cells must be positive or zero')
        if division_steps < 0: raise NameError('The number of division steps must be positive or zero')
        if cv2_div < 0: raise NameError('Target division size must be positive or zero')
        if cv2_gr < 0: raise NameError('Target growth rate must be positive or zero')
        if lamb < 0.5 or lamb > 2: raise NameError('Lambda value must be higher than 0.5 and less than 2.0')
        if not(nu in [1, 2]): raise NameError('Mu must be 1 or 2')