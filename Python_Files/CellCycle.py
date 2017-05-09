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

net = models.CellCycle_MD()
s = len(net.species)
#evl = brneval.PythonBrnEvaluator(net, x0=np.ones(s))
#reactions = evl.ratefun(np.ones(s), evl.p)

# growth: after about 25.1h (Factor 2/3) (in h)
#my_g_pdf = pdf.normal(0.02656042, 0.00005)

# growth: after about 25.1h (Factor 2/3) (in min)
my_g_pdf = pdf.normal(0.00044267, 0.00005)

# initial sizes  V, trapezoidal distribution
size_pdf_trapez = pdf.trapez(2./3, 4./3, 1.,0.5)  

# initial sizes V, gamma distribution
# size_pdf_gamma = pdf.gamma(0.5, 9.)

"""
# heterogeneous parameters
pdf0=pdf.lognormal(8000, .1, numpoints=400)  # C8 production rate
pdf2=pdf.lognormal(1000, .1, numpoints=400)  # C3
pdf4=pdf.lognormal(4000, .1, numpoints=400)   # IAP production rate,
#pdf6=pdf.lognormal(3000, .1, numpoints=400)   # CARP production rate,
hets = {'kp0':pdf0,'kp2':pdf2,'kp4':pdf4}              #,'kp6':pdf6}  
"""

# set inital states of species
initialstates = []
for i in xrange(s):
    initialstates.append(0)

net.set(dict(zip(net.species, initialstates)))

s = 0.
number = 10 # number of mother cells at t=0
t_0 = 0.
t_end = 200
n = 80

# stimulus = [np.asarray([0, 100]), [s, 0]]   # timedependent stimulus, in h!!
stimulus = [np.asarray([0, 4220]), [s, 0]] # timedependent stimulus, in min!!

maximal_number = 100000  # approx. maximal number of cells that are simulated

net.set({'Stimulus':s})   # external Stimulus, constant over time
pop = population.Population(net, number, hetparameters = {}, growth = my_g_pdf,
                            size_pdf = size_pdf_trapez, max_divergence_factor=1)

"""
# heterogeneous parameters initialized
for i in pop.mother_cells: 
    i.state[0]=i.hetvals['kp0']*i.size/6.15e-2    #C8
    i.state[2]=i.hetvals['kp2']*i.size/4.76e-2   #C3  
    i.state[4]=i.hetvals['kp4']*i.size/1e-1  #IAP
    i.state[6]=i.hetvals['kp6']*i.size/7.5e-2  #CARP
"""

t1=time.clock()
Mb_index = net.species.index('Mb')

# set Simulation parameters, time-dependent stimulus
pop.Population_Simulator(t_0 = 0, t_end = t_end, n = n, initial_states = None,
                         maximal_number = maximal_number, 
                         method = brn.simulation.VodeSimulatorStop, 
                         stimulus = {'Stimulus':stimulus}, 
                         stopcondition = lambda x: x[Mb_index] >= 1000) 

# set relevant output
t2=time.clock()
print 'time for population-simulation:  ' + str(t2-t1) + ' at s = ' + str(s)
cells = pop.get_cells()
print "Terminal total number of cells:" + str(len(cells))

alive, dead = pop.get_numbers()
 
print "Population curve:" + str(alive)

"""
VolumeSorted=pop.get_species_distribution(8,'sim_to_end')

X0valuesend = pop.get_species_distribution(0,'sim_to_end')
X1valuesend = pop.get_species_distribution(1,'sim_to_end')
X7valuesend = pop.get_species_distribution(7,'sim_to_end')
C8All=[sum(x) for x in zip(X0valuesend, X1valuesend, X7valuesend)]

X2valuesend = pop.get_species_distribution(2,'sim_to_end')
X3valuesend = pop.get_species_distribution(3,'sim_to_end')
X5valuesend = pop.get_species_distribution(5,'sim_to_end')
C3All=[sum(x) for x in zip(X2valuesend, X3valuesend, X5valuesend)]

X4valuesend = pop.get_species_distribution(4,'sim_to_end')
IAPAll=[sum(x) for x in zip(X4valuesend, X5valuesend)]

X6valuesend = pop.get_species_distribution(6,'sim_to_end')
CARPAll=[sum(x) for x in zip(X6valuesend, X7valuesend)]
"""

'''
Create plots
'''

font = {'family':'serif', 'weight':'normal', 'size':35 }

# plot Ma
for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[0] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin A/Cdk2', fontdict=font)
    
# plot Mb
matplotlib.pyplot.figure()
    
for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[1] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin B/Cdk1', fontdict=font)

# plot Md
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[2] for x in i.trajectory])
    matplotlib.pyplot.title('cyclin D/Cdk4-6', fontdict=font)
    
# plot Me
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[3] for x in i.trajectory])
    matplotlib.pyplot.title('cycln E/Cdk2', fontdict=font)
    
# plot Ma, Mb, Md, Me
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[0] for x in i.trajectory],
                           i.time_points, [x[1] for x in i.trajectory],
                           i.time_points, [x[2] for x in i.trajectory],
                           i.time_points, [x[3] for x in i.trajectory])
    matplotlib.pyplot.title('Ma, Mb, Md, Me', fontdict=font)
    
# plot Volume
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points, [x[6] for x in i.trajectory])
    matplotlib.pyplot.title('V', fontdict=font)

# plot limit cycle
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot([x[0] for x in i.trajectory], [x[1] for x in i.trajectory])
    matplotlib.pyplot.title('limit cycle, Ma over Mb', fontdict=font)
    
matplotlib.pyplot.show()