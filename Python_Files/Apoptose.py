# -*- coding: utf-8 -*-
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


net=models.apoptose_DI()
s=len(net.species)
evl = brneval.PythonBrnEvaluator(net,x0=np.ones(s))
reactions = evl.ratefun(np.ones(s),evl.p)
f=len(reactions)-2*s
d=reactions[f:f+s]
p=reactions[f+s:f+2*s]





net.set({'C8a':0.})
net.set({'C3a':0.})
net.set({'C8a_CARP':0.})
net.set({'C3a_IAP':0.})





#my_g_pdf=pdf.normal(0.00063492,0.00005)    # growth            after about 17.5h   --> without TRAIL and Pulse
my_g_pdf=pdf.normal(0.00065359,0.00005)    # growth            after about 17h

#my_g_pdf=pdf.normal(0.0006734,0.00005)    # growth            after about 16.5h  --> 500ng/ml  



#my_g_pdf=pdf.normal(0.0005698,0.00005)    # growth            after about 19.5h 
#my_g_pdf=pdf.normal(0.0005848,0.00005)    # growth            after about 19h  
#my_g_pdf=pdf.normal(0.00070015,0.00005)    # growth            after about 16.5h  
#my_g_pdf=pdf.normal(0.0006418,0.00005)    # growth            after about 18h  
#my_g_pdf=pdf.normal(0.00067956,0.00005)    # growth            after about 17h  
#my_g_pdf=pdf.normal(0.00074074,0.00005)    # growth            after about 15h  


size_pdf=pdf.trapez(2./3, 4./3, 1.,0.5)    #pdf.equal(2./3,4./3)   # initial sizes  V
pdf0=pdf.lognormal(8000, .1, numpoints=400)  # C8 production rate
pdf2=pdf.lognormal(1000, .1, numpoints=400)  # C3
pdf4=pdf.lognormal(4000, .1, numpoints=400)   # IAP production rate,
#pdf6=pdf.lognormal(3000, .1, numpoints=400)   # CARP production rate,
hets = {'kp0':pdf0,'kp2':pdf2,'kp4':pdf4}              #,'kp6':pdf6}  

#para_set 1
#net.set({'kp0':17821.516587844522})
#net.set({'kp2':647.99800106228679})
#net.set({'kp4':814.1502890599437})
#net.set({'kp6':1565.427082000798})



initialstates=[]
for i in xrange(s):
    if d[i] != 0:
        initialstates.append(p[i]/d[i])
    else:
        initialstates.append(0)

net.set(dict(zip(net.species,initialstates)))

s=0.
number=10 # number of mother cells at t=0
t_0=0.
t_end=100
#n=10000 #run 3
n=80


stimulus=[np.asarray([0,2880]),[s,0]]   # zeitabhängiger Stimulus!!

maximal_number=100000000  # approx. maximal number of cells that are simulated

#born = []
alive = []
dead = []

#dead = []
#simulated = []       # pop ist Vektor von Populationen  pop[0]=Population für erstes Element in s, nur wenn s Vektor..

cells=[]

net.set({'TRAIL':s})   # external Stimulus, constant over time
pop = population.Population(net,number,hetparameters=hets,growth=my_g_pdf,size_pdf=size_pdf,max_divergence_factor=1)


for i in pop.mother_cells:                              # heterogeneous parameters initialized
    i.state[0]=i.hetvals['kp0']*i.size/6.15e-2    #C8
    i.state[2]=i.hetvals['kp2']*i.size/4.76e-2   #C3  
    i.state[4]=i.hetvals['kp4']*i.size/1e-1  #IAP
    #i.state[6]=i.hetvals['kp6']*i.size/7.5e-2  #CARP


t1=time.clock()
C3aind = net.species.index('C3a')
pop.Population_Simulator(t_0 = 0,t_end=t_end,n=n, initial_states=None,maximal_number=maximal_number,method=brn.simulation.VodeSimulatorStop, stimulus={'TRAIL':stimulus}, stopcondition=lambda x: x[C3aind] >= 1000)  # time-dependent stimulus

t2=time.clock()
print 'time for population-simulation:  ' + str(t2-t1) + ' at s = ' + str(s)
cells = pop.get_cells()
print "Terminal total number of cells:" + str(len(cells))

alive, dead = pop.get_numbers()
 
print "Terminal number of living cells:" + str(alive[-1])
print "Population curve:" + str(alive)


#kp11valuesend = pop.get_parameter_distribution('kp11','sim_to_end')
#kp7valuesend = pop.get_parameter_distribution('kp7','sim_to_end')

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

font = {'family' : 'serif','weight' : 'normal','size'   : 45 }

#X3unsorted = pop.get_unsorted_species_distribution(3,'all')   #C3a! für alle Zellen?!
#kp0value = pop.get_unsorted_parameter_distribution('kp0','sim_to_end')  
#kp2value = pop.get_unsorted_parameter_distribution('kp2','sim_to_end')  
#kp4value = pop.get_unsorted_parameter_distribution('kp4','sim_to_end')  
#kp6value = pop.get_unsorted_parameter_distribution('kp6','sim_to_end')  


#wrong=[]
#Ind=[]
#erg=[]
#for i in X3unsorted:
#    if i > 1 or i < -1:
#        wrong.append(i)

#for i in wrong:
#    z=X3unsorted.index(i)
#    Ind.append(z)
    
#kombi=zip(X3unsorted,kp0value,kp2value,kp4value,kp6value)

#for i in Ind:
#    j=kombi[i]
#    erg.append(j)
    
#print erg


#X3unsorted = pop.get_unsorted_species_distribution(3,'sim_to_end')
#X5unsorted = pop.get_unsorted_species_distribution(5,'sim_to_end')
#C3Allunsorted=[sum(x) for x in zip(X2unsorted, X3unsorted, X5unsorted)]

#VolumeSorted = pop.get_death_distribution(9,'dead')


#for i in cells:
    #print [x[33] for x in i.trajectory][0]
    
"""   
IAPAll0 = (np.ones(alive[-1]))
for j in range(0,int(alive[-1])):
    IAPAll0[j]= pdf4.sample()[0]
IAPAll0=IAPAll0/1.16e-2

matplotlib.pyplot.hist(np.log(IAPAll)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution after TRAIL treatment',color='red')
matplotlib.pyplot.hist(np.log(IAPAll0)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution in untreated population',color='blue')
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 45 }
matplotlib.pyplot.title('Distribution of IAP molecules after 48h \n TRAIL and before treatment',fontdict=font)
matplotlib.rc('xtick',labelsize=45)
matplotlib.rc('ytick',labelsize=45)
matplotlib.pyplot.xlabel('log(Number of IAP molecules)', fontdict=font)
matplotlib.pyplot.ylabel('Cells', fontdict=font)
matplotlib.pyplot.legend(loc='upper left', prop={'size':35})
matplotlib.pyplot.show()

matplotlib.pyplot.figure()


C3All0 = (np.ones(alive[-1]))
for j in range(0,int(alive[-1])):
    C3All0[j]= pdf2.sample()[0]
C3All0=C3All0/3.9e-3

matplotlib.pyplot.hist(np.log(C3All)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution after TRAIL treatment',color='red')
matplotlib.pyplot.hist(np.log(C3All0)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution in untreated population',color='blue')
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 45 }
matplotlib.pyplot.title('Distribution of C3 molecules after 48h \n TRAIL and before treatment',fontdict=font)
matplotlib.rc('xtick',labelsize=45)
matplotlib.rc('ytick',labelsize=45)
matplotlib.pyplot.xlabel('log(Number of C3 molecules)', fontdict=font)
matplotlib.pyplot.ylabel('Cells', fontdict=font)
matplotlib.pyplot.legend(prop={'size':35})
matplotlib.pyplot.show()

matplotlib.pyplot.figure()


C8All0 = (np.ones(alive[-1]))
for j in range(0,int(alive[-1])):
    C8All0[j]= pdf0.sample()[0]
C8All0=C8All0/3.9e-3

matplotlib.pyplot.hist(np.log(C8All)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution after TRAIL treatment',color='red')
matplotlib.pyplot.hist(np.log(C8All0)/np.log(10),bins=130,alpha=0.6,histtype='bar',log=False,label='Distribution in untreated population',color='blue')
matplotlib.rc('xtick',labelsize=45)
matplotlib.rc('ytick',labelsize=45)
matplotlib.pyplot.title('Distribution of C8 molecules after 48h \n TRAIL and before treatment',fontdict=font)
matplotlib.pyplot.xlabel('log(Number of C8 molecules)', fontdict=font)
matplotlib.pyplot.ylabel('Cells', fontdict=font)
matplotlib.pyplot.legend(prop={'size':35})
matplotlib.pyplot.show()




"""


"""
for i in cells:
        matplotlib.pyplot.plot(np.array(i.time_points)/60,[x[9] + x[10] + x[12] for x in i.trajectory]) # alles C3!
        font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 50 }
        matplotlib.rc('xtick',labelsize=25)
        matplotlib.rc('ytick',labelsize=25)
        matplotlib.pyplot.title('Number of C3 molecules, no TRAIL',fontdict=font)
        matplotlib.pyplot.xlabel('Time (h)', fontdict=font)
        matplotlib.pyplot.ylabel('C3 Molecules', fontdict=font)
        
matplotlib.pyplot.figure()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               




"""

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[2] for x in i.trajectory]) # alles inaktives C3
    matplotlib.pyplot.title('Inactive Caspase 3',fontdict=font)
    
matplotlib.pyplot.figure()
    
for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[4] + x[5] for x in i.trajectory]) # alles IAP
    matplotlib.pyplot.title('IAP',fontdict=font)

matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[0] + x[1] + x[7] for x in i.trajectory]) # alles C8
    matplotlib.pyplot.title('Caspase 8',fontdict=font)
    
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[6] + x[7] for x in i.trajectory]) # alles CARP
    matplotlib.pyplot.title('CARP',fontdict=font)
    
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[8] for x in i.trajectory]) # V
    matplotlib.pyplot.title('V',fontdict=font)

matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[1] for x in i.trajectory]) # alles aktiviertes C8
    #matplotlib.pyplot.title('Activated Caspase 8, critical parameter set 1, \n kp0:17821.516587844522, kp2: 647.99800106228679, \n kp4: 814.1502890599437, kp6: 1565.427082000798',fontdict=font)
    #matplotlib.pyplot.title('Activated Caspase 8, critical parameter set 2, \n kp0:25901.436606402542, kp2:518.94908963078876, \n kp4:1187.3979700202974, kp6:951.21690718804678',fontdict=font)
    #matplotlib.pyplot.title('Activated Caspase 8, critical parameter set 3, \n kp0:20386.408300914311, kp2:2053.117121717446, \n kp4:4410.1506651759473, kp6:943.6979674645504',fontdict=font)
    #matplotlib.pyplot.title('Activated Caspase 8, critical parameter set 4, \n kp0:19107.228217953434, kp2:3930.1018871105916, \n kp4:1751.2697709913921, kp6:4034.5569553290293',fontdict=font)
    #matplotlib.pyplot.title('Activated Caspase 8, critical parameter set 5, \n kp0:15514.95093722544, kp2:2429.5010668520067, \n kp4:3409.272781858308, kp6:856.73834482337941',fontdict=font)
    matplotlib.pyplot.title('Activated Caspase 8',fontdict=font)
    
    
    
    matplotlib.rc('xtick',labelsize=45)
    matplotlib.rc('ytick',labelsize=45)


    
matplotlib.pyplot.figure()

for i in cells:
    matplotlib.pyplot.plot(i.time_points,[x[3] + x[2] + x[5] for x in i.trajectory]) # alles C3
    matplotlib.pyplot.title('All Caspase 3',fontdict=font)


matplotlib.pyplot.show()

"""

matplotlib.pyplot.bar(np.array(pop.time_points)/60,np.array(alive)/400,label='Results of simulation', linewidth=4, color='blue')     #500 cells
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 40 }
matplotlib.pyplot.title('Cell number during 48h treatment',fontdict=font)
#matplotlib.pyplot.xlim((0,9800/60))
matplotlib.rc('xtick',labelsize=35)
matplotlib.rc('ytick',labelsize=35)
matplotlib.pyplot.xlabel('Time (h)', fontdict=font)
matplotlib.pyplot.ylabel('Cells (Normalized)', fontdict=font)
matplotlib.pyplot.legend(loc='upper left',prop={'size':35})

matplotlib.pyplot.show()








pdf1=pdf.lognormal(8000, .1, numpoints=400)
stst=pdf1.varax/6.15e-2 
matplotlib.pyplot.plot(np.log(stst)/np.log(10),pdf1.pdf*pdf1.varax) #C8
matplotlib.pyplot.show()

pdf1=pdf.lognormal(1000, .28, numpoints=400)
stst=pdf1.varax/4.76e-2 
matplotlib.pyplot.plot(np.log(stst)/np.log(10),pdf1.pdf*pdf1.varax) #C3
matplotlib.pyplot.show()

pdf1=pdf.lognormal(4000, .3, numpoints=400)
stst=pdf1.varax/1e-1 
matplotlib.pyplot.plot(np.log(stst)/np.log(10),pdf1.pdf*pdf1.varax) #IAP
matplotlib.pyplot.show()

pdf1=pdf.lognormal(3000, .3, numpoints=400)
stst=pdf1.varax/7.5e-2 
matplotlib.pyplot.plot(np.log(stst)/np.log(10),pdf1.pdf*pdf1.varax) #CARP
matplotlib.pyplot.show()





import math

def funct(x):
    myfun = 1.5*np.log(2)*math.exp(2*np.log(2)*(1-(0.75*x)))   #evt Rundungsfehler??
    return myfun


q=[]
v=[]

for V in range(200/3,400/3,1):
   k=float(V)/100
   ps=funct(k)
   q.append(ps)
   v.append(k)
   print q

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 40 }
matplotlib.pyplot.plot(v,70*np.array(q),label='$\overline{p}_s (V)$ calculated', linewidth=4, color='blue') 
matplotlib.pyplot.hist(VolumeSorted,bins=500,alpha=0.6,histtype='bar',log=False,label='Simulated distribution',color='red')
matplotlib.pyplot.xlim((0.65,1.334))
matplotlib.pyplot.title('Time invariant volume distribution',fontdict=font)
matplotlib.pyplot.xlabel('Volume ($V$)', fontdict=font)
matplotlib.pyplot.ylabel('Number of cells', fontdict=font)
matplotlib.pyplot.legend(loc='upper right',prop={'size':35})
matplotlib.pyplot.show()


Vol = []
for j in range(0,500)):
    #Vol[j]= size_pdf.sample()
    
    



"""

#d=shelve.open('apopt_5000min')
#d['new']=pop
#d.close()
