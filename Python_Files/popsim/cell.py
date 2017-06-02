# -*- coding: utf-8 -*-

import brn
import numpy as np
import math
import pdf
import copy
import random
import matplotlib.pyplot
from brn import brneval
from brn import simulation

try:
    import scipy as sp
    import scipy.integrate as spi
    have_scipy=True
except ImportError:
    have_scipy=False

class Cell(object):

    def __init__(self, net, hetparameters={},growth=None, mother=None, d_size=1.95,t0=0, partitionpdf=None,x0=None,size=1.,max_divergence_factor=np.inf):

        """
        net - biochemical reaction network
        hetparameters - key:parameter value:pdf
        growth - pdf for growth rate
        d_size - size at which cell will divide
        partitionpdf - pdf which determines how mother and daughter share the molecules
        x0 - initial state
        size - initial size
        max_divergence_factor - the maximum factor, by which the heterogeneous parameter values of each daughter cell may deviate from the mother cell's value.
        Use np.inf for no corellation between daughter and mother values and 1. to enforce exactly the same values. See also divide(). A factor between 0.999 and 1.001 is treated as 1.

        """

        self.d_size=d_size
        self.net = net
        self.trajectory = []
        self.time_points = []
        self.mother = mother
        self.daughters = []
        self.alive = True
        self.hetparameters = hetparameters
        self.partitionpdf = partitionpdf
        self.max_divergence_factor = max_divergence_factor
        self.state = np.asarray([i for i in x0]) if x0 != None else np.asarray(self.net.eval(self.net.species))
        #self.state = np.asarray([int(i) for i in x0]) if x0 != None else np.asarray(self.net.eval(self.net.species))
        self.division_times = [float(t0)]         # division_times[0] = time of birth
        self.simulated_time = t0


        self.growth = growth
        self.g = self.growth.sample()[0] if growth != None else 0.
        self.size=size
        #else:
        #    self.size=None
        if mother != None:  # daughter inherits hetvals from mother, deviation allowed up to max_divergence_factor
            self.hetvals = dict.fromkeys(hetparameters.keys())
            """
            determine values of heterogeneous parameters
            """
            for j in hetparameters.keys():
                if abs(self.max_divergence_factor - 1) < 0.001:
                    self.hetvals[j] = mother.hetvals[j]
                else:
		    oldval = mother.hetvals[j]
                    self.hetvals[j] = self.hetparameters[j].sample()[0]
                    while (self.hetvals[j] < oldval/self.max_divergence_factor) or (self.hetvals[j] > oldval*self.max_divergence_factor):
                        self.hetvals[j]= self.hetparameters[j].sample()[0]
        else:
            self.hetvals = dict.fromkeys(hetparameters.keys())
            for j in hetparameters.keys():
                self.hetvals[j]= self.hetparameters[j].sample()[0]


        if 'V' in self.net.species:
            self.net.set({'V':size})
        self.removed=False



    def divide(self,t):
        """
        division at time t
        inheritance: daugther=(1-factor)*mother, mother *= factor
        """
        self.division_times.append(t)
        if self.partitionpdf == None:
	    factor = 0.5
	else:
            factor = self.partitionpdf.sample()
            factor=max(factor,1-factor)
        daughter_x= np.asarray([(1-factor)*i for i in self.state])
        #daughter_x= np.asarray([int((1-factor)*i) for i in self.state])
        new_daughter = Cell(self.net, self.hetparameters,self.growth, self, self.d_size,t0=t,partitionpdf=self.partitionpdf,x0=daughter_x,size=self.d_size*(1-factor),max_divergence_factor=self.max_divergence_factor)
        #new_daughter.hetvals = self.hetvals  # daughter has same hetvals as mother

        #new_daughter.initialize(x0=daughter_x, t0=t)
        #new_daughter.size = self.d_size*(1-factor)
        self.daughters.append(new_daughter)
        self.size = self.d_size*factor
        for j in self.hetparameters.keys():    # new hetvals are sampled for mother cell after division
            oldval = self.hetvals[j]
            if abs(self.max_divergence_factor - 1) < 0.001:
                self.hetvals[j] = oldval
            else:
                self.hetvals[j] = self.hetparameters[j].sample()[0]
                while (self.hetvals[j] < oldval/self.max_divergence_factor) or (self.hetvals[j] > oldval*self.max_divergence_factor):
                    self.hetvals[j]= self.hetparameters[j].sample()[0]
        for i in xrange(len(self.state)):
            self.state[i] -= daughter_x[i]




    def kill(self):
        self.alive = False


    def initialize(self, x0=None, t0=None): # not used

        if t0 != None:
            self.simulated_time = t0
        if x0 != None:
            self.state = x0
        else:
            self.state = self.net.eval(self.net.species)


    def simulate_stepwise(self, t0, t_end, time_points, x0, method):
        # simulates the cell's brn over interval [t0, t_end] by splitting the interval to avoid too long simulation periods
        step = 0.1
        nstep = int((t_end-t0)/step)
        times = [t0+step*i for i in xrange(nstep+1)]
        times.extend(time_points)
        times.extend([t_end])
        times.sort()
        i=0
        while i<(len(times)-1):    # remove multiple entries of times
	    if times[i+1]==times[i]:
	        times.pop(i+1)
	    else:
	        i=i+1
        t,res=self.net.simulate(times,method=method,x0=x0)
        return t,res

    def simulate_stepwisestop(self, t0, t_end, time_points, x0, method, stopcondition):
        # simulates the cell's brn over interval [t0, t_end] by splitting the interval to avoid too long simulation periods
        step = 0.1   # maximal length of one simulation interval
        nstep = int((t_end-t0)/step)
        times = [t0+step*i for i in xrange(nstep+1)]
        times.extend(time_points)
        times.extend([t_end])
        times.sort()
        i=0
        while i<(len(times)-1):    # remove multiple entries of times
	    if times[i+1]==times[i]:
	        times.pop(i+1)
	    else:
	        i=i+1
        t,res=self.net.simulate(times,method=method,x0=x0, stopcondition=stopcondition)
        return t,res

    def save_trajectory(self,time_points, t, res):
        # appends those elements of t and res which are specified in time_points to self.time_points and self.trajectory
        indices=[]
        for i in xrange(len(t)):   # find indices of t which are element of time_points
	    for s in time_points:
	        if t[i]==s:
		    indices.extend([i])

	tt = [t[i] for i in indices]    # save values only if desired
	rres = [res[i] for i in indices]
	self.time_points.extend(tt)
	self.trajectory.extend(rres)


    def simulate_interval(self, t0, t_end, time_points, x0, simulator):
        # simulates the brn over a given interval [t0, t_end] with initial-condition x0
        #if len(time_points) != 0:
        if time_points != []:
            j1 = (time_points>=t0).argmax()
            j2 = (time_points>t_end).argmax() if time_points[-1]>t_end else len(time_points)# last index in interval
            t_points = time_points[j1:j2]
        else:
	    t_points = []

        step = 0.1
        nstep = int((t_end-t0)/step)
        times = [t0+step*i for i in xrange(nstep+1)]
        times.extend(t_points)
        times.extend([t_end])
        times.sort()
        i=0
        while i<(len(times)-1):    # remove multiple entries of times
	    if times[i+1]==times[i]:
	        times.pop(i+1)
	    else:
	        i=i+1
        simulator.initialize(x0=x0, p=self.net.eval(self.net.parameters))
        t,res=simulator.run(times)
        self.save_trajectory(time_points, t, res)
        return t[-1], res[-1]


    def Cell_Simulator(self,t_end,time_points=[], t0=None ,x0=None,stimulus={},simulator=None):
        """
        event: 0: simulated until t_end; 1: cell divided; 2: cell died
        simulate cell with simulation method specified in input argument "method" until its next division
        growth rate: s(t)=s(t0)*exp(g*(t-t0))
        stimulus: dictionary, assigns a time-dependent value to one variable of the brn
          key: species or variable or parameter
          value: array, a[0][:]=timepoints (numpy array) at which value of the key changes, a[0][0] must be less then or equal to t0/self.simulated_time
                        a[1][:]=values of species or variable or parameter specified in key

        """
        if self.g != 0:
            self.net.set({'g':self.g})
            next_division = t0 + (self.d_size-self.size)/self.g   # linear cell growth      #  exponential cell growth:   t0 + np.log(self.d_size/self.size)/self.g
            t_end = min(next_division, t_end)
        self.net.set(self.hetvals)

        if 'V' in self.net.species:
            #self.net.set({'V':self.size})
            self.state[self.net.species.index('V')] = self.size

        if simulator is None:
            simulator = simulation.VodeSimulator(self.net)


        x0 = x0 if x0 is not None else self.state
        t0 = t0 if t0 is not None else self.simulated_time
        #j1=(time_points>t0).argmax()
        #j2=min((time_points>=next_division).argmax(),len(time_points)-1) # -1?
        if stimulus != {}:
	    s = stimulus.keys()[0]
	    stim = stimulus[s]
	    maxind=len(stim[0])   # last index of stimulus
	    ind = (stim[0]>t0).argmax()-1 if stim[0][-1] > t0 else maxind-1 # last index <= t0
	    self.net.set({s:stim[1][ind]})
            if s in self.net.species:
                x0[self.net.species.index(s)] = stim[1][ind]
	    ind +=1    # next index of stimulus change
	    while 1: #(ind<maxind) & (stim[0][ind]<t_end) :    # terminates if  last index of stimulus is reached or time_point of stimulus change is >= t_end
	        ###     Watch out, not tested!!!
	        if ind == maxind:
		    break
		elif (stim[0][ind]>=t_end):
		    break
	        t0, x0 = self.simulate_interval(t0, stim[0][ind], time_points, x0, simulator)   # simulated until next stimulus change
	        self.net.set({s:stim[1][ind]})
                if s in self.net.species:
                    x0[self.net.species.index(s)] = stim[1][ind]
	        if t0 < stim[0][ind]:  # simulation stopped due to stopcondition
		    self.simulated_time = t0
		    self.state = x0
		    #self.size = self.state[self.net.species.index('V')]/100 if 'V' in self.net.species else self.size   # reference-size: 100
		    event = 2
		    return [float(self.simulated_time),event]
		ind +=1

        self.simulated_time, self.state = self.simulate_interval(t0, t_end, time_points, x0, simulator)    # if there is no stimulus an simulate the rest after while-loop terminated
        if self.simulated_time < t_end:  # simulation stopped due to stopcondition
	    #self.simulated_time = t0
	    #self.state = x0
	    #self.size = self.state[self.net.species.index('V')]/100 if 'V' in self.net.species else self.size   # reference-size: 100
            self.divide(self.simulated_time)
            event=1
            self.size = 0.5*self.d_size
	    #event = 2

	elif (simulator.teststop == True):# & (self.removed == False):#(next_division == t_end)&(self.removed==False):
            self.divide(next_division)
            event=1
            self.size = 0.5*self.d_size

        else:
            event = 0
            self.size = self.state[self.net.species.index('V')] if 'V' in self.net.species else self.size
        
        return [float(self.simulated_time),event]

        """
	self.size*=np.exp(self.g*(t[-1]-self.simulated_time))
        if eval(self.net.death_condition):
	    event=2
        elif (next_division<t_end)&(self.removed==False):
            self.divide(next_division)
            event=1
        else:
	    event=0
	    self.size*=math.exp(self.g*(self.simulated_time-t0))
	"""





    """
    def Cell_Simulator(self,t0,t_end,time_points,stimulus={},method=simulation.VodeSimulator):
        #
        #event: 0: simulated until t_end; 1: cell divided; 2: cell died
        #simulate cell with SSA_Simulator until its next division
        #growth rate: s(t)=s(t0)*exp(g*(t-t0))
        #stimulus: dictionary, assigns a time-dependent value to one variable of the brn
        #  key: species or variable or parameter
        #  value: array, a[0][:]=timepoints (numpy array) at which value of the key changes, a[1][:]=values of species or variable or parameter specified in key

        #
        #next_division = t0 + np.log(self.d_size/self.size)/self.g
        #self.SSA_Simulator(t0=self.simulated_time,t_end=min(next_division,t_end),x0=self.state,time_points=time_points,stimulus=stimulus)
        #self.ODE_Simulator(t0=self.simulated_time,t_end=min(next_division,t_end),time_points=time_points,stimulus=stimulus)
        #

        self.net.set(self.hetvals)
        self.net.set(dict(zip(self.net.species,self.state)))
        p = [self.net.values[i] for i in self.net.parameters]
        j=max((time_points>=next_division).argmax(),len(time_points)) # -1?
        jj=(time_points>self.simulated_time).argmax()

        times=[self.simulated_time]
        t=[self.simulated_time]
        #
        if stimulus != {}:
	    stim=stimulus.keys()[0]
	    self.net.set({stim:stimulus[stim][1][0]})
	    # ind: next index where stimulus changes
	    ind=(stimulus[stim][0][:]>self.simulated_time).argmax() if max(stimulus[stim][0][:])>self.simulated_time else len(stimulus[stim][0][:]-1)    # determine index that of minimal timepoint >= t
	    while (ind<len(stimulus[stim][0][:])):    # loop over all indices
	        if (stimulus[stim][0][ind]>=time_points[j-1]):   # simulate until time_points[j], then break the while loop

		    times.append(time_points[jj:j])
	            t,res=self.net.simulate(times,method=method,x0=self.state)
	            times=[stimulus[stim][0][ind]]
                    self.time_points.extend(t)
                    self.trajectory.extend(res)
                    self.state=res[-1]
                    break
                k=(time_points>=stimulus[stim][0][ind]).argmax()-1   # last time_point < next change of stimulus
                times.extend(time_points[jj:k])
	        #times.append(stimulus[stim][0][ind])
	        t,res=self.net.simulate(times,method=method,x0=self.state)
	        jj=k+1
	        times=[t[-1]]#[stimulus[stim][0][ind]]
                self.time_points.extend(t)
                self.trajectory.extend(res)
                self.state=res[-1]
                self.net.set({stim:stimulus[stim][1][ind]})
                ind+=1
            if t[-1]<time_points[j-1]:        # times from stimulus ended before time_points[j] -> simulate until end
	        times.extend(time_points[jj-1:j])
	        t,res=self.net.simulate(times,method=method,x0=self.state)
                self.time_points.extend(t)
                self.trajectory.extend(res)
                self.state=res[-1]
        #
        else:
            times.extend(time_points[jj:j])
            t,res=self.net.simulate(times,method=method)
            self.time_points.extend(t)
            self.trajectory.extend(res)
            self.state=res[-1]

	self.size*=np.exp(self.g*(t[-1]-self.simulated_time))
        self.simulated_time=t[-1]
        if eval(self.net.death_condition):
	    event=2
        elif (next_division<t_end)&(self.removed==False):
            self.divide(next_division)
            event=1
        else:
	    event=0
	    self.size*=math.exp(self.g*(self.simulated_time-t0))
        return [float(self.simulated_time),event]
    """


    """
    def ODE_Simulator(self,t0,t_end,time_points,x0=None,stimulus={}):
        if not have_scipy:
            raise ImportError("Could not import scipy, scipy.integrate. These are needed by ODE_Simulator")
	self.net.set(self.hetvals)
	self.evl = brneval.PythonBrnEvaluator(self.net,p=self.net.parameters,x0=self.state)
	self.with_jacobian = self.evl.havejac
        self.ode = spi.ode(lambda t,z: self.evl.ode_rhs(z),lambda t,z: self.evl.ode_jac(z))
        self.ode.set_integrator('vode',method='bdf',with_jacobian=self.with_jacobian)
        self.ode.set_initial_value(self.evl.z0,t0)
        tend = time_points[0]
        t = t0
        i = 0
        while self.ode.successful() and self.ode.t < t_end:
            if self.ode.t >= tend:
                i += 1
                tend = time_points[i]
            self.ode.integrate(tend)
            self.time_points.append(self.ode.t)
            x=self.ode.y
            self.trajectory.append(x)
            if eval(self.net.death_condition):
	        self.kill()

    """




    """
    def SSA_Simulator(self,t0,t_end, x0=None, time_points=[],stimulus={}):

        #time_points = list of time points at which x is returned
        #stimulus: dictionary, assigns a time-dependent value to a variable of the brn
        #  key: species or variable or parameter
        #  value: array, a[0][:]=timepoints, a[1][:]=values of species or variable or parameter specified in key


        t=float(t0)
        if x0 is not None:
            x = x0
        else:
            x = self.state
        if len(time_points) != 0:
            J = True
            i = 0
        else:
            J = False
        tt=t0
        self.hetvals['V']=self.size
        while (t <= t_end)&(self.alive):
            #if t>tt:
	    # print 'time: ' + str(tt)
	    # tt+=10
	    for stim in stimulus.keys():         # evaluation of stimulus at time t
	        ind=(stimulus[stim][0][:]>=t).argmax() if max(stimulus[stim][0][:])>=t else len(stimulus[stim][0][:])    # determine index that of minimal timepoint >= t
	        self.net.set({stim:stimulus[stim][1][ind]})     # evaluate stim at this index

            a=self.net.eval(self.net.reactions,x=x,values=self.hetvals)
            a_sum=sum(a)
            a_=np.cumsum(a)
            j=(a_/a_sum > np.random.random(1)).argmax()   # first index where a_/a_sum > random is true
            tau=np.random.exponential(1/a_sum,1)

            if J  :
                if (t < time_points.max()):# k: number of time_points in (t,t+tau)
                  k =  (time_points > t).argmax() - i
                  for l in xrange(k):
                    self.trajectory.append(x.copy())
                  self.time_points.extend(time_points[i:i+k])
                  i += k
            t += tau

            if J : # append missing last entries
              if (t >= time_points.max()):
                  k=len(time_points[i:])
                  for l in xrange(k):
                    self.trajectory.append(copy.copy(x))
                  self.time_points.extend(time_points[i:])
                  J = 0

            x += self.net.stoich[:,j].transpose()
            #self.hetvals['V']*=math.exp(self.g*tau)
            self.hetvals['V']=self.size*math.exp(self.g*(t-t0))
            if eval(self.net.death_condition):
            #if  x[3]>200:    # & self.hetvals['V']<1    # boolean expression to determine whether cell dies at its current state
                self.kill()




        self.simulated_time = t
        self.state = x

        """



    def evaluate(self,t):
        """
        error if time_points=[]!!
        """
        if t<= max(self.time_points):
	    i=[j>=t for j in self.time_points].index(1)
	return self.trajectory[i]


    def get_cells(self, cell_list=[]):
        cell_list.append(self)
        for i in self.daughters:
	    i.get_cells(cell_list)






    """

    def SSA_Simulator(self,t0,t_end, x0=None, time_points=[]):

        #time_points = list of time points at which x is returned


        t=float(t0)
        if x0 is not None:
            x = x0
        else:
            x = self.state
        if len(time_points) != 0:
            J = True
            i = 0
        else:
            J = False
        while (t <= t_end):
            a=self.net.eval(self.net.reactions,x=x,values=self.hetvals)
            a_sum=sum(a)
            a_=np.cumsum(a)
            j=(a_/a_sum > np.random.random(1)).argmax()   # first index where a_/a_sum > random is true
            tau=np.random.exponential(1/a_sum,1)

            if J  :
                if (t < time_points.max()):# k: number of time_points in (t,t+tau)
                  k =  (time_points > t).argmax() - i
                  for l in xrange(k):
                    self.trajectory.append(x.copy())
                  self.time_points.extend(time_points[i:i+k])
                  i += k
            t += tau

            if J : # append missing last entries
              if (t >= time_points.max()):
                  k=len(time_points[i:])
                  for l in xrange(k):
                    self.trajectory.append(copy.copy(x))
                  self.time_points.extend(time_points[i:])
                  J = 0

            x += self.net.stoich[:,j].transpose()




        self.simulated_time = t
        self.state = x

    """










