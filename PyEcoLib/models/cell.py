import numpy as np

class Cell:
    def __init__(self, idx: int, initial_size: float, num_steps: int, growth_rate: float, div_par: float, k: float) -> None:
        
        """
        Model of cell in simulation
        
        Args:
            idx (int): Identifier of cell
            initial_size (float): Initial size of the cell
            num_steps (int): Number of steps for division
            growth_rate (float): Rate of groething of the cell
            div_par (float): Target division size of the cell
            k (float): k value
        Returns:
            self
        """
        
        self.idx = idx
        self.size = initial_size
        self.base_size = initial_size
        self.division_size = initial_size
        self.num_steps = num_steps
        self.growth_rate = growth_rate
        self.div_par = div_par
        self.k = k
        
        self.step = 0
        self.random_value = np.random.rand()
        self.num_divisions = 0
        self.next_time = 100 / self.growth_rate
        
    def change(self, size: float) -> None:
        
        """
        Change size of the cell
        
        Args:
            size (float): New size of the cell
        Returns:
            None
        """
        
        self.size = size
        self.random_value = np.random.rand()
        
    def division(self, size: float, div_par: float, growth_rate: float, k: float) -> None:
        
        """
        Divide the cell
        
        Args:
            size (float): New size of the cell
            div_par (float): New target division size of the cell
            growth_rate (float): New growth rate of the cell
            k (float): New k value of the cell
            
        Returns:
            None
        """
        
        self.base_size = self.division_size * self.div_par
        self.size = size * self.div_par
        self.division_size = size
        self.div_par = div_par
        self.growth_rate = growth_rate
        self.k = k
        
        self.random_value = np.random.rand()
        self.num_divisions += 1
        self.step = 0
        
    def json(self):
        return {
            'idx': self.idx,
            'size': self.size,
            'size_at_birth': self.base_size,
            'division_size': self.division_size,
            'num_steps': self.num_steps,
            'growth_rate': self.growth_rate,
            'div_par': self.div_par,
            'k': self.k,
            'step': self.step,
            'rv': self.random_value,
            'num_divisions': self.num_divisions,
            'next_t': self.next_time
        }
        
