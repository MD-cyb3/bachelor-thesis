# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 08:48:32 2017

@author: Michi

preliminarities
"""

'''
import all important packages
'''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
# packages of pybrn
import brn
from popsim import population
from popsim import models


# settings for font of plots and axes properties
font = {'family': 'Palatino Linotype', 'style': 'normal', 
        'weight': 'normal', 'size': 24}
font2 = {'family': 'Arial', 'style': 'normal', 'weight': 'normal', 'size': 24}
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['xtick.labelsize'] = 15 #set the value globally
mpl.rcParams['ytick.labelsize'] = 15 #set the value globally


''' 
inital setup
'''
# load cell cycle model
net = models.CellCycle_MD()
    
# pdf for growth of cells
growth_pdf = None

# initial sizes (Volume)
size_pdf = None 

# set inital states of species
initial_states = None

# simulation method --> no division, no stop condition
method = brn.simulation.VodeSimulator

# since no stimulus is needed here, its an empty dictionary
stimulus = {}

# no hetparams necessary
het_params = {}

# no division, no inheritance
max_div_fac = 1

# approx. maximal number of cells that are simulated
maximal_number = 1

''' 
set simulation parameters #1
first simulation: no division, no death --> determination of Mb treshold
'''
number = 1 # number of mother cells at t=0
t_0 = 0
t_end = 400 # in hours
n = 400 # simulation steps

# initialization of population
pop = population.Population(net, number, hetparameters = het_params, 
                            growth = growth_pdf, size_pdf = size_pdf, 
                            max_divergence_factor = max_div_fac)

# set initial states of cells
# all equal to zero --> converge to limit cycle
# just for the simulation to create all variables           
for i in pop.mother_cells:
    i.state[0:6] *= 0 
        
# timer to get simulation time
t1 = time.clock()

''' 
simulate cells #1
'''
# start simulating the cell population
pop.Population_Simulator(t_0 = t_0, t_end = t_end, n = n, 
                         initial_states = initial_states,
                         maximal_number = maximal_number, 
                         method = method, stimulus = stimulus)
                         
# get list of all cells
cells = pop.get_cells()
                         
'''
set relevant output
'''
print "First simulation to create needed variables!\n"
# second timer
t2 = time.clock()
print "time for population-simulation:  " + str(t2-t1) + "\n"
cells = pop.get_cells()
print "Number of cells simulated with:  " + str(len(cells)) + "\n"

## plot Volume
#plt.figure(figsize=(8,7))
#for i in cells:
#    plt.plot(i.time_points, [x[6] for x in i.trajectory], linewidth=2)
#    plt.plot(range(26), [1.95] *26, linewidth=4, label='Schwelle')
#plt.xlim(0, 25)
#plt.xlabel('Zeit in Stunden', fontdict=font2)
#plt.ylabel('Hilfsvolumen', fontdict=font2)
#plt.legend(loc=4)
#



'''
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
'''

# get species indices
V_index = net.species.index('V')
Ma_index = net.species.index('Ma')
Mb_index = net.species.index('Mb')
Md_index = net.species.index('Md')
Me_index = net.species.index('Me')
E2F_index = net.species.index('E2F')
Cdc20_index = net.species.index('Cdc20')


# store values of the state Mb
Mb_values = []

for i in cells:
    Mb_values.append([x[Mb_index] for x in i.trajectory])
    
Ma_values = []
Md_values = []
Me_values = []
E2F_values = []
Cdc20_values = []

for i in cells:
    Ma_values.append([x[Ma_index] for x in i.trajectory])
    Md_values.append([x[Md_index] for x in i.trajectory])
    Me_values.append([x[Me_index] for x in i.trajectory])
    E2F_values.append([x[E2F_index] for x in i.trajectory])
    Cdc20_values.append([x[Cdc20_index] for x in i.trajectory])
      
'''                
compute minimal (y-)values of Mb for each oscillaton (local minima)
'''
Mb_minima = [] # array of all local minima

# cycle_time and minimum_range ensure, that each interval contains only one
# local minimum
cycle_time = 23 # slightly less then 25, compensated by interval from
                # [starting_point - starting_point + minimum_range]
minimum_range = 23 # range in which local minimum appears

starting_points = [65] # starting point for first local minimum
                       # 65 is chosen, as the cell has to converge to the limit
                       # cycle first

# set starting points for all intervals
for i in range(1, (t_end / cycle_time)):
    if (59 + i * cycle_time) > t_end:
        # otherwise last point of simulation is taken,
        # which is not necessary a minimum   
        break 
    else:
        starting_points.append(65 + i * cycle_time)
    
# store minima of y-values of Mb of the cell in array 'Mb_minima'
for i in range(len(starting_points)-1): # -1 because of last starting_point
    Mb_minima.append(np.min(Mb_values[0][starting_points[i]:
                                         starting_points[i] + minimum_range]))
    
Mb_threshold = np.max(Mb_minima) # value for stopcondition

# print result
print "Maximum value of all local minima of Mb:  " + str(Mb_threshold) + "\n\n"
    
''' 
set simulation parameters #2
second simulation, no division, no death --> creation of limit cycle trajectory
'''
number = 1 # number of mother cells at t=0
t_0 = 0
t_end = 200 # in hours
n = 50000 # simulation step (40 * t_end to reach desired resolution)

# initialization of population
pop = population.Population(net, number, hetparameters = het_params, 
                            growth = growth_pdf, size_pdf = size_pdf, 
                            max_divergence_factor = max_div_fac)
                            
# set initial states of cells
# all equal to zero --> converge to limit cycle
# just for the simulation to create all variables           
for i in pop.mother_cells:
    i.state[0:5] *= 0 

# timer to get simulation time
t1 = time.clock()

''' 
simulate cells #2
'''
# start simulating the cell population
pop.Population_Simulator(t_0 = t_0, t_end = t_end, n = n, 
                         initial_states = initial_states,
                         maximal_number = maximal_number, 
                         method = method, stimulus = stimulus)
                         
'''
set relevant output
'''
print "Second simulation to create limit cycle trajectory!\n"
# second timer
t2 = time.clock()
print "time for population-simulation:  " + str(t2-t1) + "\n"
cells = pop.get_cells()
print "Number of cells simulated with:  " + str(len(cells)) + "\n"

'''
store values of state Mb
'''
Mb_values = []

for i in cells:
    Mb_values.append([x[Mb_index] for x in i.trajectory])
    
# t_0 is chosen at an arbitrary minimum of Ma,
# t_end is the subsequent minimum of Ma --> one cell cycle
t_0 = 39960
t_end = 45852

# time for one limit cycle
limit_cycle_period = t_end-t_0

# store limit cycle trajectories of all variables
for i in cells:
    Ma_trajectory = [x[Ma_index] for x in i.trajectory[t_0:t_end]]
    Mb_trajectory = [x[Mb_index] for x in i.trajectory[t_0:t_end]]
    Md_trajectory = [x[Md_index] for x in i.trajectory[t_0:t_end]]
    Me_trajectory = [x[Me_index] for x in i.trajectory[t_0:t_end]]
    E2F_trajectory = [x[E2F_index] for x in i.trajectory[t_0:t_end]]
    Cdc20_trajectory = [x[Cdc20_index] for x in i.trajectory[t_0:t_end]]
    V_trajectory = [x[V_index] for x in i.trajectory[0:t_end-t_0]]
        
    
# save values of one limit cycle of all five states (for reduced phase model)
limit_cycle_trajectory = zip(Ma_trajectory, Mb_trajectory, Me_trajectory,
                             E2F_trajectory, Cdc20_trajectory)
                             
# save values of one limit cycle of all six variables
# for initialization of cell population
limit_trajectories = zip(Ma_trajectory, Mb_trajectory, Md_trajectory, 
                         Me_trajectory, E2F_trajectory, Cdc20_trajectory)