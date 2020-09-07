# PyEcoLib

![logo](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/logo.png)

PyEcoLib (Python Library for E. coli size dynamics estimation) is library to estimate bacterial cell size stochastic dynamics including time-continuous growth process and division events.

Some of the basic uses PyEcoLyb:
* Estimating stochastic division times of E. coli bacteria with arbitrary precision from measurable parameters like growth rate, mean cell size and number of division steps. These times can be coupled to any stochastic simulation algorithm for gene expression.
* Estimating numerically the dynamics of the mean and the coefficient of variance of the size-distribution for continuously growing and dividing bacteria.
* Computing the trends in added size versus the size at division for different division strategies (Adder, timer-like, sizer-like) depending on the parameter ![equation](https://latex.codecogs.com/gif.latex?%5Clambda). 
* Including variability in the septal position during cell splitting. This can be used to model the molecule segregation.
* Set an arbitrary distribution of initial sizes.
* Set variability in cell-to-cell growth rate.


## Language
python version 3.8.1

## Libraries
PyEcoLib uses as libraries:
* [Numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [math](https://docs.python.org/3/library/math.html) 
* [Platform](https://docs.python.org/2/library/platform.html)


## Get Started
* ### Object "Cell"
The object Cell is defined in cell.py. To initialize a Cell it must be defined the following parameters:

>*  idx: The cell index.
>*  total_steps:  The number of steps to trigger the division (A typical value is 15).
>*  V0: Size at the beginning of the simulation (A typical value is 1 femptoliter).
>*  gr: Growth rate in this cycle (In optimal growth conditions, a typical value is ln(2)/18 min-1).
>*  divpar: division parameter. This corresponds to the ratio between the size at the end of the cycle and the newborn cell (it is 0.5 by default.)
>*  k: Rate of division steps occurrence (We consider it to have the same value than the growth rate).   

* ### Object "Simulator"

Simulator is the main object in the library. This object has most of the function that can be used to study cell size dynamics.
```
from PyEcoLib.simulator import Simulator
simulator = Simulator(ncells, gr, sb, steps, CV2div = 0, CV2gr = 0, lamb=1, V0array=None)

```

> ### Simulator parameters
To start a new simulator, the user must define the following parameters:

Required parameters:

>* ncells: Number of cells to simulate.
>* gr: Growth Rate (log(2)/doubling time).
>* sb: Mean newborn cell Size.
>* steps: Division Steps that trigger the cell splitting.

Optional parameters:

>* V0array: Array with the initial sizes of the cells (By default, all cells start at size sb). This array has to have the same dimension than the number of cells.
>* CV2div: Squared coefficient of variation of septal position. By default, this parameter is zero (A typical value is 0.001).
>* CV2gr: CV2gr. Squared coefficient of variation of cell-to-cell growth rate (zero by default).  (A typical value is 0.02)
>* lamb: Parameter of division strategy. ![equation](https://latex.codecogs.com/gif.latex?%5Clambda%3C0%0D%0A). ![equation](https://latex.codecogs.com/gif.latex?%5Clambda%3D1%0D%0A%0D%0A)defines the adder (by default).  ![equation](https://latex.codecogs.com/gif.latex?0%3C%5Clambda%3C1%0D%0A) is timer-like and .  ![equation](https://latex.codecogs.com/gif.latex?1%3C%5Clambda%3C%5Cinfty%0D%0A) is sizer-like. In optimal growth conditions lamb=1 and in slow growing E. coli cells, lamb is close to 1.5. Some other rod-shaped cells show different lamb between 0.5 and 2.

### Implemented functions:

* ### Obtaining the stochastic size dynamics for all the cells in the simulation.

```
Simulation.szdyn(tmax, sample_time, nameCRM = "./dataCRM.csv")
```
Defining a maximum time tmax and a sampling time sample_time with units of inverse growth rate, the function returns a file with default name "./dataCRM.csv". 

The first row is the time from 0 to tmax, sampled periodically with period sample_time. Subsequent columns correspond to the size of each cell at those times.

| time      | Cell1          | Cell2  |Cell3  |
|----------|----------|----------|----------|
|0|	3|	3|	3|
|1.8	|3.2153	|3.2153	|3.2153|
|3.6	|3.446	|3.446	|3.446|
|5.4	|3.6934	|3.6934	|3.6934|


* ### Estimating numerically the trends of cell size dynamics

```
Simulation.szdynFSP(tmax, CV2sz = 0, nameFSP = "./dataFSP.csv")
```
Estimate numerically the dynamics of the mean and variance of the size distribution with default name "./dataFSP.csv". The variability in the starting cell size can be set by the parameter CV2sz corresponding to the square coefficient of variation of the size.


|time	|Meansize|	VarSize|
|-----|-----|-----|
|0.9	|3.105794772|	1.60E-11|
|1.08	|3.127397282|	9.57E-11|
|1.26	|3.149150051|	4.34E-10|
|1.44	|3.171054121|	1.61E-09|

* ### Simulating the division strategy 

```
Simulation.divstrat(tmax, sample_time, nameDSM = "./dataDSM.csv")
```
This function runs a simulation similar to szdyn producing a file with default name "./dataDSM.csv" the first row is the size at birth, the second row is the size at division and the third column corresponds to the time when division occurs


|S_b	|S_d	|time|
|-----|-----|-----|
|2.8486	|6.0399	|37.8|
|1.5	|7.0832	|37.8|
|3.5025|6.2172	|37.8|

* ### Estimating numerically the division strategy
```
Added,cv2=Simulation.SdStat(sb)
```
Returns an array consisting on the mean added size at division ("Added") and the squared coefficient of variation of this added size ("cv2"). 

## Examples included in the library

* ### How run an example?
#### Required libraries to run the examples
* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [pandas](https://pandas.pydata.org/)
* [scipy](https://www.scipy.org/)

### Size Statistics

```
from examples.SizeStatistics import SizeStatistics
```

Found in "SizeStatistics" folder, this example shows how to use the library to plot the dynamics of the cell size. 

Using
```
Simulation.szdyn(tmax, sample_time, nameCRM = "./dataCRM.csv")
```
where you include the time to simulate the cells "tmax" in units of inverse growth rate, the samling time "sample_time" and the name of the output file ("./dataCRM.csv" by default), you can obtain a file that, after some statistical analysis, can be used to estimate the mean size and its coefficient of variation. 

Using

```
Simulation.szdynFSP(tmax, CV2sz = 0, nameFSP = "./dataFSP.csv")
```
Where you have to set the maximum time to simulate and, optionally, the variability in the initial size distribution "CV2sz" (by default is zero and a typical value is 0.015)

You obtain a file with exactly the mean and variance of the sizes. You can also obtain the size dynamics of single cell.


Gathering both results, you should obtain something like this plot:


![SizeStatistics](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/SizeStatistics/figures/size_statistics.png)


On the other hand, if you are interested not in the size dynamics but in the relationship between the added size and size at birth, you can use

```
Simulation.divstrat(tmax, sample_time, nameDSM = "./dataDSM.csv")
```

To obtain the data of size at division vs size at birth for some cycles. This data, after some statistical analysis, can be used to estimate these trends.

But if you prefer to make the estimation numerically, using  

```
Added,cv2=Simulation.SdStat(sb)
```
you can obtain directly these trends with arbitrary precision. The result of merging the simulated data (dots) and the numeric estimation (lines) should look something like this:


![DivStrategy](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/SizeStatistics/figures/div_strategy.png)



### GeneExpressionDirectSSA

```
from examples.GeneExpressionSSA import GeneExpressionDirectSSA
```

Found in "GeneExpressionDirectSSA" folder, this example shows how to implement a basic simulation of gene transcription and RNA translation considering the cell division.

Together with the reaction times that you simulation calculates, with PyEcoLib, you can estimate the time to the next division step to happen using:

```
Simulation.get_next_t(n)
```

Before division, you can estimate the division parameter, this is, how small will be the descendant cell regarding the size at division of the current cell. 

```
Simulation.get_dp(n)
```

This division parameter is centered at 0.5 and has a stochastic variability quantified by its coefficient of variation "CV2div" (by default is zero) which is set once the user defines the simulator.
Hence, during division, molecules can segregate following a binomial distribution with parameter equal to this division parameter.

This example shows you how to calculate the molecule number and the molecule concentrations for many cells.  

The following plot shows you how the dynamics of the main statistics of both, RNA and protein look like: 


![MergeStatistics](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/GeneExpressionDirectSSA/figures/merge_statistics2.png)




Other functions

The following functions are the basic function implemented in the examples mentioned above but they are not already optimized for using by an user with no experience. We recomend to use them once you are familiar with the library:

Basic functions:

* ### Initialize cells
```
Simulation.initialize_cells(V0array)
```
This function initialize cells. If V0array is not given, all cells start with size sb. the length of V0array must be equal to the number of cells in the simulator.

* ### Simulate cells
```
Simulation.simulate(tmax)
```
Before run this function, make sure you have already initialized the cells. Simulation is performed during a time interval with duration tmax>0.

* ### Obtain the time to the next division step
```
Simulation.get_next_t(n)
```
Returns the time to the next division step for the cell with index n.

* ### Obtain the growth rate of any cell
```
Simulation.get_gr(n)
```
Returns the growth rate of the cell with index n.

* ### Estimate how any cell will get split
```
Simulation.get_dp(n)
```
Returns the division parameter (the inverse between the size at division and the size of the new cell) for the cell with index n (0.5 by default).




## Advanced Examples

### NoisyDiv
```
from PythonExamples.NoisyDiv import NoisyDiv
```

Located in "NoisyDiv" folder, using this example you can learn how to introduce noise in either the septum position and the cell-to-cell growth rate:

```
CV2sz = 0.015
v0 = mean_size*np.random.gamma(shape=1/CV2sz,scale=CV2sz,size=ncells)
simulator = Simulator(ncells, gr, sb, steps, CV2div = 0.001, CV2gr = 0.01, lamb=1, V0array=v0, sample_time = 0)

```
where we have included a noise in septum position (with squared coefficient of variation of 0.001) and noise in growth rate (with magnitude of 0.01). Also, an initial size distribution (Beta distributed with CV2=0.015).

Unfortunately, you cannot obtain the numerical approximation to this trend since our theory cannot consider these additional noise sources.

As a result, after a statistical analysis, you should obtain the following figure: 

![SizeStatisticNoisy](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/NoisyDiv/figures/size_statisticsnoisy.png)


### Sizerlike


```
from PythonExamples.SizerLike import Sizerlike
```

Located in "SizerLike" folder, using this example, you will learn how to study the division strategy changing the division parameter lamb.


The main change in this script is to define the Simulator with lamb = 2 (Sizer-like division strategy in this case). This parameter, lamb, shows stability if 0.5<lamb<2.

```
sim=Simulator(ncells=ncells, gr = gr, sb=mean_size, steps = div_steps, lamb = 2)

```
The simulator, now called "sim", will last a few additional seconds trying to calculate the adequate rate of division steps to obtain the desired sb.


The division strategy can be simulated using the function:

```
sim.divstrat(tmax = tmax, sample_time = 0.1*doubling_time, nameDSM = "./data/dataDSM.csv")

```
where the user has to select the sampling time and the path of the file to export the data. 

You can also find the added size and the stochastic variability around this value, as function of the size at birth sb, using the function:

```
Adder,cv2=sim.SdStat(sb)
```
As a result, after some statistical analysis, you should obtain a plot like this:

![DivStrategySizerlike](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/SizerLike/figures/div_strategy_sizerlike.png)

You can also simulate the size dynamics of a given number of cells using the function
```
sim.szdyn(tmax = tmax, sample_time= 0.1*doubling_time, nameCRM = "./data/dataCRM.csv")
```
that gives you a file called "./data/dataCRM.csv" with the time and the size of all the cells at each instant. 

Numerically, you can also estimate the numerical values of the mean size and its variance using the function:

```
sim.szdynFSP(tmax = tmax, nameFSP = "./data/dataFSP.csv")
```
That makes a file "./data/dataFSP.csv" showing directly the mean and the variance of the cell size along the time.

When you merge these both results, you should obtain a plot like this:


![SizeStatisticsSizer](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/SizerLike/figures/size_statisticssizer.png)


### AdvancedSizeStatistics

```
from examples.AdvancedSizeStatistics import AdvancedSizeStatistics
```
Let's combine most of the properties of cell-size that can be studied using PyEcoLib!

This example includes the calculations made in the examples:

* Size Statistics
* SizeLike
* NoisyDiv

Exploring different division strategies, you can obtain different patterns in added size vs size at birth and its noise some characteristic values (lab=0.5 timer-like, lamb=1 adder and lamb=2 sizer-like) should seem like this:

![FullDivStrategy](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/AdvancedSizeStatistics/figures/full_div_strategy.png)

The comparison on size dynamics including different sources of noise should look like this:

![SizeStatisticsComp](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/AdvancedSizeStatistics/figures/size_statistics_comp1.png)

Combining all together, you can plot the following graph:

![FullSizeStatisticsComparison](https://raw.githubusercontent.com/SystemsBiologyUniandes/PyEcoLib/master/examples/AdvancedSizeStatistics/figures/full_size_statistics_comparison.png)

