import json
import math
import numpy as np
from scipy import integrate
from scipy import optimize as opt
from scipy.stats import gamma

class Utils:
    
    def __init__(self, growth_rate: float, total_steps: int, base_size: float, lamb: float) -> None:
        """
        Utils class to calculate operations for simulation
        
        Args:
            growth_rate (float): Rate of growthing of the cells
            total_steps (int): Number of steps before division
            k (float): k value
            lamb (float): Lambda value
            
        Returns:
            None
        """
        
        self.growth_rate = growth_rate
        self.total_steps = total_steps
        self.base_size = base_size
        self.lamb = lamb
        
        if self.lamb == 1:
            self.k = self.total_steps * self.growth_rate / self.base_size
        else:
            self.k = self.total_steps * self.get_k()
        
        
    def get_next_time(self, initial_size: float, random_value: float, growth_rate: float, k: float) -> float:
        """
        Method that calculate the cell next time division
        
        Args:
            initial_size (float): Initial size of the cell
            random_value (float): Random value
            growth_rate (float): Growth rate of the cell
            k (float): k value of the cell
            
        Returns:
            float: Next time of division
        
        """
        
        mu: float = self.growth_rate * growth_rate
        k: float = self.k * k
        return (1/(self.lamb * mu)) * np.log(1-((self.lamb * mu) / (k * initial_size ** self.lamb)) * np.log(random_value))
    
    def get_size_base(self, k: float) -> float:
        """
        Returns size base of the cells based on a K value
        
        Args:
            k (float): k value
        Returns:
            float: Size base of the cells
        """
        
        def root(tt):
            return self.multimean(tt, k) - 2 * tt
        
        def mean_sb():
            return opt.bisect(root, 0.00001, 100000)
        return mean_sb() 
    
    def get_k(self, size_base: float) -> float:
        """
        Return k when it cannot be calculate with equation growth_rate / size_base
        
        Args:
            size_base (float): Size base of the cells
        Returns:
            float: Returns k value
    
        """
        return opt.bisect(self.opti, 0.001, 1.5)
    
    def multimean(self, size: float, k: float) -> float :
        """
        Method calculate the probability density in different moments
        
        Args:
            size (float): Size of cells
            k (float): k value
        
        Returns:
            float: Returns the calculation of the probability
        """
        
        def moment(division_size: float):
            return self.rhomulti(size, division_size, k) * division_size
        return integrate.quad(moment, size, np.inf)[0]
        
    
    def rhomulti(self, size_base: float, division_size: float, k: float) -> float:
        """
        Method that calculates the probability density function of 
        the Gamma distribution based of base size and division size
        
        Args:
            size_base (float): Size base of the cells
            division_size (float): Target size to cell division
            k (float): k value
        
        Returns:
            float: Returns probability density of the Gamma distribution
        
        """  
        
        c = self.total_steps * k / self.growth_rate
        x = c * ((division_size ** self.lamb - size_base ** self.lamb) / self.lamb)
        return gamma.pdf(x, self.total_steps) * c * division_size ** (self.lamb - 1) 
    
    @staticmethod
    def new_growth_rate(cv2: float) -> float :
        """
        Calculate the new growth rate based on the target growth rate
        
        Args:
            cv2 (float): Target growth rate
        
        Returns:
            float: New growth rate
        """
        
        if cv2 == 0: return 1
        return np.random.gamma(shape=1/cv2, scale=cv2)
    
    @staticmethod
    def new_division_par(cv2: float) -> float:
        """
        Calculate a new division size
        
        Args:
            cv2 (float): Target division size
        Returns:
            float: New division size
        """
        
        if cv2 >= 1: raise NameError('The param target division size must be less than 1')
        if cv2 < 0: raise NameError('The param target division size has to be greater than 0')
        
        if cv2 == 0: return 0.5
        
        beta = 0.5 * ((1/cv2) - 1)
        return np.random.beta(a=beta, b=beta)
    
    @staticmethod
    def json_output(idx: int, time: float, size: float) -> json:
        """
        Method to map the json output
        
        Args:
            idx (int): Index of cell
            time (float): Time of simulation
            size (float): Size of cell
            
        Returns:
            json: Get json output of size cells
        """
        
        return {
            'idx': idx,
            'time': time,
            'size': size
        }
        
    @staticmethod
    def json_pop_output(idx: int, time: float, size: float, num_steps: int) -> json:
        """
        Method to map the json output
        
        Args:
            idx (int): Index of cell
            time (float): Time of simulation
            size (float): Size of cell
            num_steps (int): Number of steps before division
            
        Returns:
            json: Get json output of size cells
        """
        
        return {
            'idx': idx,
            'time': time,
            'size': size,
            'steps': num_steps
        }
    
    @staticmethod
    def truncate(num: float, ciphers: int) -> float:
        """
        Function that truncate sizes
        
        Args:
            num (float): Number to truncate
            ciphers (int): Number of digits in the number
        
        Returns:
            float: Number truncated
        """
        
        pos: float = pow(10.0, ciphers)
        return math.trunc(pos * num) / pos