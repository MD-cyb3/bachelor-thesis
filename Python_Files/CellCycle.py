# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 08:48:32 2017

@author: Michi

cell cycle simulation
"""
import numpy as np
#import math
import brn
from brn import brneval


from popsim import cell
from popsim import population
from popsim import propfun
from popsim import pdf
from popsim import models

import matplotlib.pyplot
import copy
import time
import shelve


''' 
inital setup
'''

net = models.CellCycle_MD()
s = len(net.species)
#evl = brneval.PythonBrnEvaluator(net, x0=np.ones(s))
#reactions = evl.ratefun(np.ones(s), evl.p)

# growth: after about 25.1h (Factor 2/3) (in h)
#my_g_pdf = pdf.normal(0.02656042, 0.00005)

# growth: after about 25.1h (Factor 2/3) (in min)
my_g_pdf = pdf.normal(0.00044267, 0.00005)

# initial sizes  V, trapezoidal distribution
size_pdf_trapez = pdf.trapez(0., 2., 1., 0.5)  


"""
# heterogeneous parameters
pdf0=pdf.lognormal(8000, .1, numpoints=400)  # C8 production rate
pdf2=pdf.lognormal(1000, .1, numpoints=400)  # C3
pdf4=pdf.lognormal(4000, .1, numpoints=400)   # IAP production rate,
#pdf6=pdf.lognormal(3000, .1, numpoints=400)   # CARP production rate,
hets = {'kp0':pdf0,'kp2':pdf2,'kp4':pdf4}              #,'kp6':pdf6}  
"""

# set inital states of species
# for Ma and Mb values on limit cycle, for other variables 0.01 (paper)
initialstates = []
for i in xrange(s):
    if i == 0:
        initialstates.append(2.5)
    elif i == 1:
        initialstates.append(0.852)
    else:
        initialstates.append(0.01)
        
net.set(dict(zip(net.species, initialstates)))


''' 
set simulation parameters
'''

s = 0.
number = 50 # number of mother cells at t=0
t_0 = 0.
t_end = 100
n = 100

#stimulus = [np.asarray([0, 100]), [s, 0]]   # timedependent stimulus, in h!!
stimulus = [np.asarray([0, 4220]), [s, 0]] # timedependent stimulus, in min!!

maximal_number = 1000000  # approx. maximal number of cells that are simulated

net.set({'Stimulus':s})   # external Stimulus, constant over time
pop = population.Population(net, number, hetparameters = {}, growth = my_g_pdf,
                            size_pdf = size_pdf_trapez, max_divergence_factor=1)

t1=time.clock()
V_index = net.species.index('V')
Mb_index = net.species.index('Mb')


''' 
simulate cells
set stopcondition with species V and Mb coupled
stopcondition = lambda x: x[V_index] >= 1.95 and x[Mb_index] <= 0.0008 
'''

# set Simulation parameters, time-dependent stimulus
pop.Population_Simulator(t_0 = 0, t_end = t_end, n = n, initial_states = None,
                         maximal_number = maximal_number, 
                         method = brn.simulation.VodeSimulatorStop, 
                         stimulus = {'Stimulus':stimulus}, 
                         stopcondition = lambda x: x[V_index] >= 1.95 and x[Mb_index] <= 0.0008)


'''
set relevant output
'''

t2=time.clock()
print 'time for population-simulation:  ' + str(t2-t1) + ' at s = ' + str(s)
cells = pop.get_cells()
print "Terminal total number of cells:" + str(len(cells))

alive, dead = pop.get_numbers()
 
print "Terminal number of living cells:" + str(alive[-1])
print "Population curve:" + str(alive)


'''
compute minimal (y-)values of Mb for each oscillaton (local minima)
'''
'''
Mb_minima = [] # array of all local minima
Mb = [] # array of y-values of trajectory of Mb
cycle_time = 22 
# slightly less then 25, compensated by interval from
# starting_point until starting_point + minimum_range

minimum_range = 20 # range in which local minimum appears
starting_points = [12] # starting point for first local minimum

 # set starting points for all local minima
for i in range(1, (t_end / cycle_time)):
    starting_points.append(12 + i * cycle_time)
    
# store y-values of trajecotry of Mb in array 'Mb'
for i in cells:
    Mb.append([x[Mb_index] for x in i.trajectory])
    
# store minima of y-values of Mb of all cells in array 'Mb_minima'
for ind in range(len(Mb)):
    for i in range(len(starting_points)-1): # -1 because last minimum is no real minimum (end of simulation)
        Mb_minima.append(np.min(Mb[ind][starting_points[i]:starting_points[i] + minimum_range]))
    
# compute maximum of all minima of all cells for stopcondition
Mb_maxima = []
for i in range(len(Mb_minima)):
    Mb_maxima.append(np.max(Mb_minima[i]))
Mb_maxima = np.max(Mb_maxima) # nearly always < 0.0008 --> value for stopcondition
'''

'''
Create plots
'''

font = {'family':'serif', 'weight':'normal', 'size':35 }
'''
# plot Ma
for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[0] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin A/Cdk2', fontdict=font)
'''    
# plot Mb
matplotlib.pyplot.figure()
    
for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[1] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin B/Cdk1', fontdict=font)
'''
# plot Md
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[2] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin D/Cdk4-6', fontdict=font)
    
# plot Me
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[3] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin E/Cdk2', fontdict=font)
    
# plot Ma, Mb, Md, Me
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[0] for x in i.trajectory],
                           i.time_points, [x[1] for x in i.trajectory],
                           i.time_points, [x[2] for x in i.trajectory],
                           i.time_points, [x[3] for x in i.trajectory])
    matplotlib.pyplot.title('Ma, Mb, Md, Me', fontdict=font)
'''    
# plot Volume
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[6] for x in i.trajectory])
    matplotlib.pyplot.title('V', fontdict=font)

# plot limit cycle
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot([x[0] for x in i.trajectory], [x[1] for x in i.trajectory])
    matplotlib.pyplot.title('limit cycle, Mb over Ma', fontdict=font)
    
    
matplotlib.pyplot.show()
