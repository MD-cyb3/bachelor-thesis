# -*- coding: utf-8 -*-
import pdf
import brn
from brn import simulation
import cell
import numpy as np
import copy



class Population(object):




    def __init__(self, net, number=1, hetparameters={},partitionpdf=None,growth=None,division_size=1.95,size_pdf=None, max_divergence_factor=np.inf):
        """
        dictionary hetparameters with key: parameter  value: pdf
        population represents list of cells with reaction network net and heterogeneous parameters from hetparameters
        growth: pdf
        division_size: cell divides if size > division_size
        """

        self.current_n = number
        self.net=net
        self.hetparameters = hetparameters
        self.time_points=[]
        self.partitionpdf = partitionpdf
        self.cells=[]
        self.mother_cells=[]
        self.max_divergence_factor = max_divergence_factor
        sizes=[float(size_pdf.sample()) for i in xrange(number)] if size_pdf !=None else np.ones(number)
        for i in xrange(number):
            self.cells.append(cell.Cell(net,hetparameters,growth=growth, mother=None, d_size=division_size,partitionpdf=self.partitionpdf,size=sizes[i],max_divergence_factor=max_divergence_factor))
            self.cells[-1].state *= sizes[i]     #letzte entstehende Zelle jeweils neue Initialcondition
            self.mother_cells.append(self.cells[i])



    def division(self,i):
        """
        division of cell i
        """
        self.cells.append(self.cells[i].daughters[-1])
        self.s_times.append(self.cells[-1].simulated_time)
        self.current_n += 1
        r=-1
        if self.current_n > self.maximal_number:
            r = int(np.random.random(1)*(self.current_n))
            #self.delete(r)
        return r



    def delete(self,i):
         # cell = self.cells[i]
         # cell.mother.daughters.pop(cell.mother.daughters.index(cell)) 	# remove from mother cell
         self.cells[i].removed=True
         self.cells.pop(i)
         self.current_n -= 1
	 self.s_times.pop(i)


    def kill(self, i):
        """
        kills cell i
        """
        self.cells[i].kill



    def Population_Simulator(self,t_0 = 0,t_end=0,division_size=4./3,n=11, initial_states=None, maximal_number=100, stimulus={},method=simulation.VodeSimulator, **simargs):

        """
        while not all cells are simulated until t_end:
        - determine minimal simulated_time of all cells and this cell's index
        - simulate this cell, simulation is terminated if this cell divides
        """

        #
        """
        time_scale=1
        if method==simulation.VodeSimulator:
	    time_scale=10
	self.net.scale=[time_scale*i for i in self.net.scale]
        """
        #
        self.s_times=[0]*self.current_n  #np.zeros(self.number)
        self.maximal_number = maximal_number
        for i in xrange(self.current_n):
            self.s_times[i]=float(self.cells[i].simulated_time)
        st=min(self.s_times)
        k=self.s_times.index(st)
        if initial_states != None:
            for i in xrange(self.current_n):
                self.cells[i].initialize(x_0=initial_states,t_0=t_0)
        time = np.linspace(t_0,t_end,n)
        self.time_points.extend(time)
        self.n_divisions = [0]*(n-1)
        self.n_deaths = [0]*(n-1)
        self.n_removed = [0]*(n-1)
        self.pop_size=[self.current_n]*n
        st = float(t_0)
        simulator = method(self.net, **simargs)
        while (st < t_end):
            i_1 = (time>=st).argmax() if time.max() >= st else n
            [self.s_times[k],event]=self.cells[k].Cell_Simulator(t_end=t_end,time_points=time[i_1:],t0=st,x0=None,stimulus=stimulus,simulator=simulator)
            idx = [tp >= self.s_times[k] for tp in self.time_points].index(True) - 1
            if event==1:   # cell divided
    	         r = self.division(k)
    	         self.n_divisions[idx] += 1
    	         if r >= 0:
    		     self.delete(r)
    		     self.n_removed[idx] +=1
#           if simulator.teststop:   # cell divided
#                self.cells[k].divide(st)
#                r = self.division(k)
#                self.n_divisions[idx] += 1
#                if r >= 0:
#                    self.delete(r)
#                    self.n_removed[idx] +=1
#           elif event==2 and simulator.teststop: # cell divided
#               self.cells[k].divide(st)
#               r = self.division(k)
#               self.n_divisions[idx] += 1
#               if r >= 0:
#                   self.delete(r)
#                   self.n_removed[idx] +=1
	    elif event==2:   # cell death
	        self.delete(k)
	        self.n_deaths [idx] +=1
            if self.s_times != []:
                st=min(self.s_times)   # st - minimal simulated time of all cells, which are simulated
                k=self.s_times.index(st)
            else:
	        st = t_end+1  # population died out

	#cells=self.get_cells


    def get_cells(self):
        """
        list of all cells
        """
        cells=[]
        for i in self.mother_cells:
	    i.get_cells(cells)
	return cells


    def get_numbers(self):
        """
        Number of born, dead, alive and simulated cells
        """
        born = np.zeros(len(self.time_points))      # all cells that were born
        born[1:len(born)] += np.cumsum(self.n_divisions)
        dead = np.zeros(len(self.time_points))      # cells that died
        dead[1:len(dead)] += np.cumsum(self.n_deaths)
        alive = len(self.mother_cells)*np.ones(len(self.time_points))     # living cells = mother_cells + born - dead
        alive += born-dead
        simulated = alive.copy()                       # simulated cells = alive - removed
        simulated[1:len(simulated)] -= np.cumsum(self.n_removed)
        return alive, dead
        #return born, dead, alive, simulated







    def get_species_distribution(self, index, mode='sim_to_end'):
        """
        (error if parameter is not a heterogeneous parameter)
        gets the distributions of species in the population (end of simulation)
        parameter: heterogeneous parameter
        mode: string indicating which cells to consider. 'all': all cells, 'alive': only alive cells. 'not_rem': only cells that were not removed,  'sim_to_end': only cells that were simulated until the end

        returns an array with the sorted values of the species
        """
        c = self.get_cells()
        i=len(c)-1
        if mode == 'alive':
            while i >= 0:
	        if not c[i].alive:   # throw away dead cells
	            c.pop(i)
                i -= 1
	elif mode == 'not_rem':
            while i >= 0:
	        if c[i].removed:   # throw away removed cells
	            c.pop(i)
	        i -= 1
	elif mode == 'sim_to_end':
	    while i >= 0:
	        if c[i].simulated_time < self.time_points[-1]:   # throw away cells that were not simulated until the end
	            c.pop(i)
	        i -= 1
	elif mode == 'deadcells':
	    while i >= 0:
	        if c[i].alive:   # throw away alive cells, not checked!!
	            c.pop(i)
	        i -= 1
        vals = [0]*len(c)
        for i in xrange(len(c)):
	    vals[i] = c[i].state[index]
	vals.sort()
	return vals
	
	
	
	
    def get_death_distribution(self, index, mode='death'):
        """
        (error if parameter is not a heterogeneous parameter)
        gets the distributions of species in the population (end of simulation)
        parameter: heterogeneous parameter

        returns an array with the sorted values of the species
        """
        c = self.get_cells()
        i=len(c)-1
        if mode == 'death':
            while i >= 0:
	        if not c[i].dead:   # throw away alive cells
	            c.pop(i)
                i -= 1
        vals = [0]*len(c)
        for i in xrange(len(c)):
	    vals[i] = c[i].state[index]
	vals.sort()
	return vals
	
	

    def get_unsorted_species_distribution(self, index, mode='sim_to_end'):
        """
        (error if parameter is not a heterogeneous parameter)
        gets the distributions of species in the population (end of simulation)
        parameter: heterogeneous parameter
        mode: string indicating which cells to consider. 'all': all cells, 'alive': only alive cells. 'not_rem': only cells that were not removed,  'sim_to_end': only cells that were simulated until the end

        returns an array with the unsorted values of the species
        """
        c = self.get_cells()
        i=len(c)-1
        if mode == 'alive':
            while i >= 0:
	        if not c[i].alive:   # throw away dead cells
	            c.pop(i)
                i -= 1
	elif mode == 'not_rem':
            while i >= 0:
	        if c[i].removed:   # throw away removed cells
	            c.pop(i)
	        i -= 1
	elif mode == 'sim_to_end':
	    while i >= 0:
	        if c[i].simulated_time < self.time_points[-1]:   # throw away cells that were not simulated until the end
	            c.pop(i)
	        i -= 1
	elif mode == 'sim_not_to_end':
	    while i >= 0:
	        if not c[i].simulated_time < self.time_points[-1]:   # throw away cells that were not simulated until the end
	            c.pop(i)
	        i -= 1
        vals = [0]*len(c)
        for i in xrange(len(c)):
	    vals[i] = c[i].state[index]
	return vals


    def get_parameter_distribution(self, parameter, mode='sim_to_end'):
        """
        error if parameter is not a heterogeneous parameter
        gets the distributions of heterogeneous parameters in the population
        parameter: heterogeneous parameter
        mode: string indicating which cells to consider. 'all': all cells, 'alive': only alive cells. 'not_rem': only cells that were not removed,  'sim_to_end': only cells that were simulated until the end

        returns an array with the sorted values of the parameter
        """
        c = self.get_cells()
        i=len(c)-1
        if mode == 'alive':
            while i >= 0:
	        if not c[i].alive:   # throw away dead cells
	            c.pop(i)
                i -= 1
	elif mode == 'not_rem':
            while i >= 0:
	        if c[i].removed:   # throw away removed cells
	            c.pop(i)
	        i -= 1
	elif mode == 'sim_to_end':
	    while i >= 0:
	        if c[i].simulated_time < self.time_points[-1]:   # throw away cells that were not simulated until the end
	            c.pop(i)
	        i -= 1
        vals = [0]*len(c)
        for i in xrange(len(c)):
	    vals[i] = c[i].hetvals[parameter]
	vals.sort()
	return vals
	
	
    def get_unsorted_parameter_distribution(self, parameter, mode='sim_t_end'):
        """
        error if parameter is not a heterogeneous parameter
        gets the distributions of heterogeneous parameters in the population
        parameter: heterogeneous parameter
        returns an array with the unsorted values of the parameter
        """
        c = self.get_cells()
        i=len(c)-1
        vals = [0]*len(c)
        for i in xrange(len(c)):
	    vals[i] = c[i].hetvals[parameter]
	return vals
	
	



    """
    def get_weighted_numbers(self):
        Number of born, dead, alive and simulated cells, weighted according to the number of removed cells. The simulated cells represent also the cells that were removed.
        sum_removed = np.cumsum(self.n_removed)
        weight = [(float(self.maximal_number)+r)/self.maximal_number for r in sum_removed]
        born = np.zeros(len(self.time_points))      # all cells that were born
        born[1:len(born)] += np.cumsum(np.multiply(self.n_divisions, weight))
        dead = np.zeros(len(self.time_points))      # cells that died
        dead[1:len(dead)] += np.cumsum(self.n_deaths)
        alive = len(self.mother_cells)*np.ones(len(self.time_points))     # living cells = mother_cells + born - dead
        alive[1:len(alive)] += born-alive
        simulated = alive                       # simulated cells = alive - removed
        simulated[1:len(simulated)] -= np.cumsum(self.n_removed)
        return born, dead, alive, simulated








    def get_numbers(self):

    #pop_n - total number of cells that were born before time_points
    #sim_n - number of cells that are simulated at time points
    #cell is counted only if it was not removed immediately after birth

        cells=self.get_cells()
        l=len(self.time_points)
        sim_n=[0]*l
        pop_n=[0]*l
        for i in cells:
	    t=i.time_points
	    if t!=[]:
                i1=[tp >= t[0] for tp in self.time_points].index(True)
                i2=[tp >= t[-1] for tp in self.time_points].index(True)
                for j in xrange(i1,i2+1):
	            sim_n[j]+=1
	        for j in xrange(i1,l):
		    pop_n[j]+=1


        return [pop_n,sim_n]
        """
