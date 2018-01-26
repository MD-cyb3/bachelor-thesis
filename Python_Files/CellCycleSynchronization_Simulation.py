# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 08:49:37 2017

@author: Michi

cell cycle synchronization
"""


'''
import all important packages
'''
import numpy as np
import brn
from popsim import population
from popsim import pdf
from popsim import models
import matplotlib.pyplot as plt
from math import log

import os

import save_values

'''
compute point psi on limit cycle with minimum distance to point z 
in neighbourhood of limit cycle (mapping 1) and corresponding phase on unit
circle (mapping 2) --> reduced phase model
'''
def mapping(z, limit_cycle):
    # array with vectors from point z to points on limit cycle
    sub = np.array(z) - np.array(limit_cycle)
    # get index of point xi on limit cycle with least euclidean distance
    xi_index = np.argmin(np.einsum('ij,ij->i',sub,sub))
    xi = limit_cycle[xi_index]
    # frequency of oscillators (T = len(limit_cycle))
    omega = (2 * np.pi) / len(limit_cycle)
    # corresponding phase on unit circle
    theta_rad = omega * xi_index
    # representation as complex number on unit circle
    r = 1
    j = complex(0, 1)
    theta_unit_circle = r * np.exp(j * theta_rad)
    return theta_rad, theta_unit_circle

''' 
set simulation parameters #3
third simulation with division of cells at stopcondition
'''
# turn interactive mode off
plt.ioff()

number = 3000 # number of mother cells at t=0
t_0 = 0.
t_end = 1 # in hours
n = 2 # simulation steps

# time for one period
T = 2*np.pi
# factor for multiplication
# factor = limit_cycle_period / (2*np.pi)
factor = 937.5
# set growth rate
growth_rate = log(2) / T

'''
initialize cells in steady state distribution
'''
# create initial pdf
initial_ss_pdf = pdf.my_pdf(0., T, growth_rate)
# sample time points for initialzation of cell states 
# from steady state age distribution
ss_time_points = []
for i in xrange(number):
    ss_time_points.append(initial_ss_pdf.sample(num=0))
    
# multiply time points with caclulated factor, to distribute the cells on the
# whole limit cycle (not just in the interval [2, 2 * np.pi])
ss_time = []
for i in ss_time_points:
    ss_time.append(int(round(i * factor)))
    
'''    
# plot steady-state age structured cell distribution
plt.figure()
plt.plot(initial_pdf.varax, initial_pdf.pdf)
plt.xlabel('cell age on unit circle', fontdict=font)
plt.ylabel('steady state distribution', fontdict=font)
plt.xlim(0, T)
plt.ylim(0, 0.24)
'''    
    
'''
initialize cells in normal distribution
'''
# create inital pdf
initial_n_pdf = pdf.normal(np.pi, 1.5, numpoints=100, 
                           varax=np.linspace(0., 2*np.pi, num=100))
# sample time points for initialzation of cell states 
# from normal distribution
n_time_points = []
for i in xrange(number):
    n_time_points.append(initial_n_pdf.sample(num=0))
    
# multiply time points with caclulated factor, to distribute the cells on the
# whole limit cycle (not just in the interval [2, 2 * np.pi]) 
n_time = []
for i in n_time_points:
    n_time.append(int(round(i * factor)))
    
'''  
# plot normal age cell distribution
plt.figure()
plt.plot(initial_n_pdf.varax, initial_n_pdf.pdf)
plt.xlabel('cell age on unit circle', fontdict=font)
plt.ylabel('normal distribution', fontdict=font)
plt.xlim(0, T)
'''

# values of model states at cell division
division_states = min(limit_cycle_trajectory, key = lambda t: t[1])

'''
initial heterogeneous parameters
'''
# center of lognormal distribution is a single value of the various states
# at cell division
pdf_Ma = pdf.lognormal(division_states[0], 0.1, numpoints=400)
pdf_Mb = pdf.lognormal(division_states[1], 0.1, numpoints=400)
pdf_Me = pdf.lognormal(division_states[2], 0.1, numpoints=400)
pdf_E2F = pdf.lognormal(division_states[3], 0.1, numpoints=400)
pdf_Cdc20 =pdf.lognormal(division_states[4], 0.1, numpoints=400)
het_params__2 = {'Ma': pdf_Ma, 'Mb': pdf_Mb, 'Me': pdf_Me, 'E2F': pdf_E2F, 
                 'Cdc20': pdf_Cdc20}

# cell cycle model
net = models.CellCycle_MD()

# number of simulation loops
maximum = 200

# kappa for von Mises distribution
kappa = 70

# values for simulation study
# div_fac_values = [1.005, 1.1, 2., 10., 100.]
# epsilon_values = [0.01, 0.03, 0.1, 0.5, 1.0]
div_fac_values = [10., 100.]
epsilon_values = [0.03, 0.05]

for max_div_fac in div_fac_values:
    print "Div_Fac = " + str(max_div_fac) + "\n"
    
    for epsilon in epsilon_values:
        print "Epsilon = " + str(epsilon) + "\n"
        
        # store moments of cell population
        moments = [] # absolute values
        m1_values = [] # normal values
        
        # store calculated input
        input_values = []        
        
        '''
        set folder and file names to save values
        '''
        folder_name = 'Epsilon_%.3f_DivFac_%.3f' % (epsilon, max_div_fac)
                  
        file_name_1 = 'Time_Dependent_Values'
        file_name_2 = 'Distribution_Values'
        file_name_3 = 'Phase_Dependent_Values'
        file_name_4 = 'Information about this folder'
        
        information_text = 'SimTime_%d, Kappa_%d, NumberOfCells%d' % (maximum-1, kappa, number)
        information_text += os.linesep
        information_text += 'Time_Dependent_Values with input u(t), moments m1(t), and absolut of moments |m1(t)|'
        information_text += os.linesep
        information_text += 'Distribution_Values with distribution n(x,0), n(x,t_end), p(x,0), p(x,t_end)'
        information_text += os.linesep
        information_text += 'Phase_Dependent_Values with phase theta(0), theta(t_end)'
        
        name_of_path = 'csv_files_simulation_SECONDstudy/' + folder_name + '/'
        name_of_path_plot_folder = name_of_path + 'Plots/'
        
        # create directory for plots
        directory_plots = os.path.dirname(name_of_path_plot_folder)
        if not os.path.exists(directory_plots):
            os.makedirs(directory_plots) 
        
        ''' 
        simulate cells #3
        set stopcondition with species V and Mb coupled
        stopcondition = lambda x: x[V_index] >= 1.95 and x[Mb_index] <= Mb_threshold 
        '''
        # start simulation loop
        for loop in range(1, maximum):
            # output for user
            print "Loop = " + str(loop) + "\n"     
        
            '''
            initialization of all required variables
            '''
            # all 6 states (+ volume), just the values of one simulation loop
            Ma = []
            Mb = []
            Md = []
            Me = []
            E2F = []
            Cdc20 = []
            Volume = []
            
            '''
            set growth, Volume sizes and initial states of simulation #3
            '''
            # simulation method: simulator with stopcondition
            method__2 = brn.simulation.VodeSimulatorStop
            
            # pdf for growth of cells
            growth_pdf__2 = None
            
            # set inital states of species
            initial_states__2 = None
            
            # initial volume
            size_pdf__2 = None
            
            # no stimulus required
            stimulus__2 = {}
            
            # maximal number of cells is same number as initial cells
            maximal_number__2 = number
        
            # trajectorie of limit cycle for intialization
            trajectories = limit_trajectories
            
            '''
            inheritance at cell division
            '''
            # (factor = 1 --> mother and daughter identical)
            # (factor = np.inf --> heterogeneous parameters are sampled from complete 
            # distrubtion), factor is here chosen as small as possible to ensure 
            # correct distribution of cells on cell cycle
            max_div_fac__2 = max_div_fac
            
            # set input
            if loop == 1:
                net.set({'input': 0.})
            else:
                net.set({'input': input})
                
            # initialize cell population
            pop = population.Population(net, number, hetparameters = het_params__2,
                                        growth = growth_pdf__2, size_pdf = size_pdf__2, 
                                        max_divergence_factor=max_div_fac__2)
            
            '''
            set initial states of population and initial volume
            '''
            # first simulation loop
            if loop == 1:
                # initialize volume
                Volume_initialization = []
                for i in range(number):
                    Volume_initialization.append(V_trajectory[n_time[i]])
                    
                for i in pop.mother_cells:
                    # initialize volume (size) and all 6 state-variables
                    i.size = Volume_initialization[pop.mother_cells.index(i)]
                    for index in range(len(i.state)-1):
                        i.state[index] = \
                                trajectories[n_time[pop.mother_cells.index(i)]][index]
        
            # other simulation loops
            else:
                for i in pop.cells:
                    # initialize volume
                    i.size = Volume_end[pop.cells.index(i)]
                    
                    # initialize all 6 state-variables
                    for index in range(len(i.state)-1):
                        i.state[index] = initialstates[pop.cells.index(i)][index] 
        
            # start simulating the cell population
            pop.Population_Simulator(t_0 = t_0, t_end = t_end, n = n, 
                                 initial_states = initial_states__2,
                                 maximal_number = maximal_number__2, 
                                 method = method__2, stimulus = stimulus__2, 
                                 stopcondition = lambda x: x[V_index] >= 1.95 and 
                                                           x[Mb_index] <= Mb_threshold)
        
            '''
            create variables of all states (+ volume) and store their values 
            '''
            # size of cell at the end for initialization of next simulation loop
            Volume_end = [] 
            
            for i in pop.cells:
                Ma.append([x[Ma_index] for x in i.trajectory])
                Mb.append([x[Mb_index] for x in i.trajectory])
                Md.append([x[Md_index] for x in i.trajectory])
                Me.append([x[Me_index] for x in i.trajectory])
                E2F.append([x[E2F_index] for x in i.trajectory])
                Cdc20.append([x[Cdc20_index] for x in i.trajectory])
                Volume.append([x[V_index] for x in i.trajectory])
                Volume_end.append(i.size)
        
            number = len(Ma)
        
            '''
            compute solution of ode (last state of cells)
            compute theta on unit circle via mapping function
            comput initialstates for next simulation-loop
            '''
            solution = [0] * number # solution of ODEs
            initialstates = [0] * number # initial states for next loop
            # variables for unit_circle and theta
            theta = [0] * number
            theta_unit_circle = [0] * number
        
            # for initialization of next simulation loop    
            for index in range(number): 
                # compute solution of ODEs
                solution[index] = (Ma[index][-1], Mb[index][-1], Me[index][-1], 
                                   E2F[index][-1], Cdc20[index][-1])
                # phase theta of the oscillators (cells)
                theta[index], theta_unit_circle[index] = \
                            mapping(solution[index], limit_cycle_trajectory)
                # save initialstates as last states of recent loop        
                initialstates[index] = (Ma[index][-1], Mb[index][-1], Md[index][-1],
                                        Me[index][-1], E2F[index][-1], 
                                        Cdc20[index][-1])
                        
            theta = [theta[index] for index in range(len(Ma))]
            
            # split theta_unit_circle in real and imaginary
            Re = []
            Im = []
            for ind in range(len(Ma)):    
                Re.append(theta_unit_circle[ind].real)
                Im.append(theta_unit_circle[ind].imag)
            
            # plot phase unit circle
            if loop == 1:
                theta_begin = theta
#                plt.figure(figsize=(8,7.6))
#                plt.plot(Re, Im, marker = 'o', linestyle = 'None')
#                plt.xlabel('Re', fontdict=font)
#                plt.ylabel('Im', fontdict=font)
#                plt.xlim(-1, 1)
#                plt.ylim(-1, 1)
#                plt.grid(True)  
            
            if (loop == maximum-1):
                theta_end = theta
#                plt.figure(figsize=(8,7.6))
#                plt.plot(Re, Im, marker = 'o', linestyle = 'None')
#                plt.xlabel('Re', fontdict=font)
#                plt.ylabel('Im', fontdict=font)
#                plt.xlim(-1, 1)
#                plt.ylim(-1, 1)
#                plt.grid(True)  
                
            varax = np.linspace(0., 2*np.pi, num=100)
            # kernel density estimation
            theta_pdf = pdf.estimate(theta, method="kde", varax=varax, 
                                     h=0.05, points=100)
                                     
            # von Mises estimation
            theta_pdf_mises = pdf.estimate(theta, method="vonMises", varax=varax, h=kappa,
                                           points=100)      
        
            '''
            import phase response curve (PRC)
            '''
            data = open('Phase_Response_Curve_Z_v2010_startCb.dat')
            data = data.read().split('\t')
            data_list = []
            for index in range(len(data)):
                data_list.append(data[index].split('\n'))
                
            # x in [0, 2* np.pi], z = values of PRC
            x = []
            z = []
            for index in range(len(data_list)-2):
                x.append(float(data_list[index+1][1]))
                z.append(float(data_list[index+2][0]))
                
        #    # plot PRC
        #    if (loop == 1):
        #        plt.figure(figsize=(8,7))
        #        plt.plot(x,z)
        #        plt.xlabel('phase', fontdict=font)
        #        plt.ylabel('phase response curve', fontdict=font)
        #        plt.xlim(0, 2*np.pi)    
                
            '''
            transformation n(x,t) --> p(x,t), see paper       
            '''        
            q_tilde = 2 * growth_rate * np.exp(-growth_rate * np.array(x))
            theta_tilde_pdf = theta_pdf_mises.pdf/q_tilde
            transformed_theta_pdf = \
                theta_tilde_pdf / np.trapz(np.array(theta_tilde_pdf), x=np.array(x))
        
            # plot estimated distribution graph  and save first and last distribution                        
            if (loop == 1):
                theta_pdf_mises_begin = theta_pdf_mises.pdf
                transformed_theta_pdf_begin = transformed_theta_pdf
                fig_density_begin = plt.figure(figsize=(8,7))
                plt.plot(theta_pdf_mises.varax, theta_pdf_mises.pdf)
                plt.xlabel('phase of cells', fontdict=font)
                plt.ylabel('density on unit circle', fontdict=font)
                plt.xlim(0, 2*np.pi)
                plt.ylim(0, 0.60)
                fig_density_begin.savefig('csv_files_simulation_SECONDstudy/' + folder_name + '/Plots/Density_Begin.png')
                plt.close(fig_density_begin)
                
            if (loop == maximum-1):
                theta_pdf_mises_end = theta_pdf_mises.pdf
                transformed_theta_pdf_end = transformed_theta_pdf
                fig_density = plt.figure(figsize=(8,7))
                plt.plot(theta_pdf_mises.varax, theta_pdf_mises.pdf)
                plt.xlabel('phase of cells', fontdict=font)
                plt.ylabel('density on unit circle', fontdict=font)
                plt.xlim(0, 2*np.pi)
                plt.ylim(0, 0.60)    
                fig_density.savefig('csv_files_simulation_SECONDstudy/' + folder_name + '/Plots/Density_End.png')
                plt.close(fig_density)
            
            # complex number
            j = complex(0, 1)
            k = 1
        
            # compute first circular moment
            y = np.exp(j * k * np.array(x)) * transformed_theta_pdf
            m_1 = np.trapz(y, x=np.array(x))
            m_minus1 = np.conj(m_1)
            
            m1_values.append(m_1)
            moments.append(abs(m_1))
        
            '''
            calculate input u
            '''
            d_minus1_complex = (growth_rate - j) * m_1
            d_0_complex = m_1 * m_minus1 * (-2) * growth_rate
            d_1_complex = (growth_rate + j) * m_minus1
            
            d_minus1 = d_minus1_complex * np.exp(j * -1 * np.array(x))
            d_0 = d_0_complex.real
            d_1 = d_1_complex * np.exp(j * 1 * np.array(x))   
            
            d_p = d_minus1 + d_0 + d_1    
            y1 = z * d_p * transformed_theta_pdf
            
            input_complex = (np.trapz(y1, x=np.array(x)))
            input_real = input_complex.real
            input_values.append(input_real)
            
            if input_real < 0:
                input = 0.
            else:
                # epsilon = strength of input
                input = epsilon * input_real 
                
            # plot limit cycle
            if loop == maximum-1:
                fig_limit_cycle = plt.figure(figsize=(8,7))
                for i in range(len(Ma)):
                    plt.plot(Ma[i], Mb[i], marker = 'o', linestyle = 'None')
                plt.xlabel('cyclin A-CDK2', fontdict=font)
                plt.ylabel('cyclin B-CDK1', fontdict=font)
                fig_limit_cycle.savefig('csv_files_simulation_SECONDstudy/' + folder_name + '/Plots/LimitCycle_End.png')
                plt.close(fig_limit_cycle)
        
        '''
        Create plots
        '''
        # linear interpolation of absolute of first circular moment, for 300 loops
        #x = range(len(moments))
        #xvals_list = [2.54, 17.96, 43.0, 64.94, 88.99, 115.025, 137.24, 162.0, 185.2, 
        #              210.24, 232.45, 281.12]
        #xvals = np.array(xvals_list)
        #interp_y = moments
        #interp = np.interp(xvals, x, interp_y)
        
        # plot absolute of first circular moment
        fig_moments = plt.figure(figsize=(8,7))
        plt.plot(range(len(moments)), moments)
        #plt.plot(xvals, interp, linewidth=3, label='linear interpolated')
        plt.xlabel('simulation time in hours', fontdict=font)
        plt.ylabel('absolute moments', fontdict=font)
        plt.legend(loc=4)
        fig_moments.savefig('csv_files_simulation_SECONDstudy/' + folder_name + '/Plots/Moments.png')
        plt.close(fig_moments)       
        
        # plot input
        fig_input = plt.figure(figsize=(8,7))
        plt.plot(range(len(input_values)), input_values)
        plt.xlabel('simulation time in hours', fontdict=font)
        plt.ylabel('input u', fontdict=font)
        fig_input.savefig('csv_files_simulation_SECONDstudy/' + folder_name + '/Plots/Input.png')
        plt.close(fig_input)        
        
        '''
        save all values for parameter study
        '''
        save_values.save_time_dependent_values(folder_name, file_name_1, input_values,
                                               m1_values, moments)
                                            
        save_values.save_distribution_values(folder_name, file_name_2, varax, 
                                             theta_pdf_mises_begin, theta_pdf_mises_end,
                                             transformed_theta_pdf_begin,
                                             transformed_theta_pdf_end)
                                             
        save_values.save_phase_dependent_values(folder_name, file_name_3, 
                                                theta_begin, theta_end)
                                             
        save_values.save_information_file(folder_name, file_name_4, information_text)

'''
-----------other plots if necessary, first: cells = pop.get_cells()------------
'''

'''
# plot Ma
plt.figure()
for i in range(len(Ma)):
    plt.plot(range(len(Ma[i])), Ma[i])
    
# plot Limit Cycle
plt.figure()
for i in range(len(Ma)):
    plt.plot(Ma[i], Mb[i])

# plot Ma
plt.figure()
for i in cells:
    plt.plot(i.time_points, [x[0] for x in i.trajectory])
    plt.title('cyclin A/Cdk2', fontdict=font)

# plot Mb
plt.figure()
    
for i in cells:
    plt.plot(i.time_points, [x[1] for x in i.trajectory])
    plt.title('cyclin B/Cdk1', fontdict=font)

# plot Md
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[2] for x in i.trajectory])
    plt.title('cyclin D/Cdk4-6', fontdict=font)
    
# plot Me
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[3] for x in i.trajectory])
    plt.title('cyclin E/Cdk2', fontdict=font)
    
# plot E2F
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[3] for x in i.trajectory])
    plt.title('E2F', fontdict=font)

# plot Cdc20
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[3] for x in i.trajectory])
    plt.title('Cdc20', fontdict=font)

# plot Ma, Mb, Md, Me
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[0] for x in i.trajectory],
                           i.time_points, [x[1] for x in i.trajectory],
                           i.time_points, [x[2] for x in i.trajectory],
                           i.time_points, [x[3] for x in i.trajectory])
    plt.title('Ma, Mb, Md, Me', fontdict=font)
    
# plot Volume
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[6] for x in i.trajectory])
    plt.title('V', fontdict=font)


# plot Volume-Mb
plt.figure()

for i in cells:
    plt.plot(i.time_points, [x[6]-x[1] for x in i.trajectory])
    plt.title('V-Mb', fontdict=font)

# plot limit cycle
plt.figure()

for i in cells:
    plt.plot([x[0] for x in i.trajectory], [x[1] for x in i.trajectory])
    plt.title('limit cycle, Mb over Ma', fontdict=font)
    
    
plt.show()
'''
