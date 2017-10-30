# -*- coding: utf-8 -*-
"""
brn.test.generic: generic utilities for test code

Copyright (C) 2010 Steffen Waldherr waldherr@ist.uni-stuttgart.de
"""

import numpy as np
import brn


def CellCycle_MD():
    """
    5-variable skeleton model for mammalian cell cycle
    """
    stoich = np.zeros((7,14))
    stoich[0,0] = 1
    stoich[0,4] = -1
    stoich[0,13] = 1
    stoich[1,1] = 1
    stoich[1,5] = -1
    stoich[2,2] = 1
    stoich[2,6] = -1
    stoich[3,3] = 1
    stoich[3,7] = -1
    stoich[4,8] = 1
    stoich[4,10] = -1
    stoich[5,9] = 1
    stoich[5,11] = -1
    stoich[6,12] = 1
    
    species = ['Ma', 'Mb', 'Md', 'Me', 'E2F', 'Cdc20', 'V']
    variables = []
    reactions = ['v0', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', 'v9', 
                 'v10', 'v11', 'v12', 'v13']
    parameters = ['vsa', 'vsb', 'vsd', 'vse', 'Vda', 'Vdb', 'Vdd', 'Vde', 
                  'V1e2f', 'V2e2f', 'V1cdc20','V2cdc20', 'Kda', 'Kdb', 'Kdd', 
                  'Kde', 'K1e2f', 'K2e2f', 'K1cdc20', 'K2cdc20', 'Kgf', 
                  'E2Ftot', 'Cdc20tot', 'GF', 'g', 'Stimulus', 'Vref', 'alpha', 
                  'vu', 'nu', 'ku', 'input']
    net = brn.Brn(stoich,species,reactions,parameters,variables)
    net.set({'v0':'vsa * E2F', 'v1':'vsb * Ma', 'v2':'vsd * (GF / (Kgf + GF))',
             'v3':'vse * E2F', 'v4':'Vda * Cdc20 * (Ma / (Kda + Ma))',
             'v5':'Vdb * Cdc20 * (Mb / (Kdb + Mb))', 
             'v6':'Vdd * (Md / (Kdd + Md))', 
             'v7':'Vde * Ma * (Me / (Kde + Me))',
             'v8':'V1e2f * ((E2Ftot - E2F) / \
                            (K1e2f + E2Ftot - E2F)) * (Md + Me)',
             'v9':'V1cdc20 * Mb * ((Cdc20tot - Cdc20) / \
                                   (K1cdc20 + Cdc20tot - Cdc20))',
             'v10':'V2e2f * (E2F / (K2e2f + E2F)) * Ma', 
             'v11':'V2cdc20 * (Cdc20 / (K2cdc20 + Cdc20))',
             'v12':'alpha * (Vref - V)', 
             'v13': '(vu * (1 - Ma**(nu) / (ku**(nu) + Ma**(nu)))) * input'})
    net.set({'vsa':0.175, 'vsb':0.21, 'vsd':0.175, 'vse':0.21, 'Vda':0.245, 
             'Vdb':0.28, 'Vdd':0.245, 'Vde':0.35, 'V1e2f':0.805, 'V2e2f':0.7, 
             'V1cdc20':0.21, 'V2cdc20':0.35, 'Kda':0.1, 'Kdb':0.005, 'Kdd':0.1, 
             'Kde':0.1, 'K1e2f':0.01, 'K2e2f':0.01, 'K1cdc20':1., 'K2cdc20':1., 
             'Kgf':0.1, 'E2Ftot':3., 'Cdc20tot':5., 'GF':1., 'g':1., 'V':1., 
             'Vref':2., 'alpha':0.25, 'vu':1.6, 'nu':8., 'ku':4., 'input': 0.})
    return net




def apoptose_DI():
    """
    with v11 in equation 2, with stimulus, with volume
    v0 is new with respect to Eissing: input TRAIL
    Turnover adapted - protein number doubled during 17h
    maximal turnover: 10⁵ molecukes per hour (Schwanhaeuser2013)
    """
    stoich=np.zeros((9,8))
    stoich[0,0]=-1
    stoich[1,0]=1
    stoich[0,1]=-1
    stoich[1,1]=1
    stoich[2,2]=-1
    stoich[3,2]=1
    stoich[1,3]=-1
    stoich[6,3]=-1
    stoich[7,3]=1
    stoich[1,4]=1
    stoich[6,4]=1
    stoich[7,4]=-1
    stoich[3,5]=-1
    stoich[4,5]=-1
    stoich[5,5]=1
    stoich[3,6]=1
    stoich[4,6]=1
    stoich[5,6]=-1
    stoich[4,7]=-1

    stoich = np.concatenate((stoich,np.identity(9)*-1),axis=1)    #degradation
    stoich = np.concatenate((stoich,np.identity(9)),axis=1)		#production
    species=['C8','C8a','C3','C3a','IAP','C3a_IAP','CARP','C8a_CARP','V']   # 'V' as species
    variables=[]
    reactions = ['v0','v1','v2','v3','v4','v5','v6','v7']
    degradation =  ['d0','d1','d2','d3','d4','d5','d6','d7','d8']
    production =   ['p0','p1','p2','p3','p4','p5','p6','p7','p8']
    reactions.extend(degradation)
    reactions.extend(production)
    parameters=['k0','k1','k2','k3','k4','k5','k6','k7','TRAIL','g','kp0','kp2','kp4','kp6']
    net = brn.Brn(stoich,species,reactions,parameters,variables)
    net.set({'v0':'k0*TRAIL*C8','v1':'k1*C3a*C8/V','v2':'k2*C8a*C3/V','v3':'k3*C8a*CARP/V','v4':'k4*C8a_CARP','v5':'k5*C3a*IAP/V','v6':'k6*C3a_IAP','v7':'k7*C3a*IAP/V'})
    net.set({'p0':'kp0*V','p1':0.,'p2':'kp2*V','p3':0.,'p4':'kp4*V','p5':0.,'p6':'kp6*V','p7':0.,'p8':'g'})
    net.set({'d0':'6.15e-2*C8','d1':'5.8e-3*C8a','d2':'4.76e-2*C3','d3':'5.8e-3*C3a','d4':'1e-1*IAP','d5':'1.73e-2*C3a_IAP','d6':'7.5e-2*CARP','d7':'1.16e-2*C8a_CARP','d8':0.})
    # old #net.set({'d0':'3.9e-3*C8','d1':'5.8e-3*C8a','d2':'3.9e-3*C3','d3':'5.8e-3*C3a','d4':'1.16e-2*IAP','d5':'1.73e-2*C3a_IAP','d6':'1.0e-3*CARP','d7':'1.16e-2*C8a_CARP','d8':0.})
    net.set({'k0':5.8e-5,'k1':1e-5,'k2':5.8e-5,'k3':5e-4,'k4':2.1e-1,'k5':5e-4,'k6':2.1e-1,'k7':3.0e-4,'g':1.,'V':1.,'TRAIL':0.,'kp0':8e3,'kp2':1e3,'kp4':4e3,'kp6':3e3})
    # old #net.set({'k0':5.8e-5,'k1':1e-5,'k2':5.8e-5,'k3':5e-4,'k4':2.1e-1,'k5':5e-4,'k6':2.1e-1,'k7':3.0e-4,'g':1.,'V':1.,'TRAIL':0.,'kp0':5.07e2,'kp2':8.19e1,'kp4':4.64e2})
    return net




def apoptose_net_1():
    """
    #without v11 in equation 2
    """
    stoich=np.zeros((8,13))
    stoich[0,1]=-1
    stoich[0,8]=-1
    stoich[1,1]=1
    stoich[1,4]=-1
    stoich[2,0]=-1
    stoich[2,9]=-1
    stoich[3,0]=1
    stoich[3,2]=-1
    stoich[3,5]=-1
    stoich[4,2]=-1
    stoich[4,3]=-1
    stoich[4,7]=-1
    stoich[5,2]=1
    stoich[5,6]=-1
    stoich[6,10]=-1
    stoich[6,11]=-1
    stoich[7,10]=1
    stoich[7,12]=-1
    species=['C8','C8_star','C3','C3_star','IAP','C3_star_IAP','BAR','C8_star_BAR']
    reactions = ['v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12','v13']
    parameters=['k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12','k13','k_1','k_2','k_3','k_4','k_5','k_6','k_7','k_8','k_9','k_10','k_11','k_12','k_13']
    net = brn.Brn(stoich,species,reactions,parameters)
    net.set({'v1':'k1*C8_star*C3','v2':'k2*C3_star*C8','v3':'k3*C3_star*IAP-k_3*C3_star_IAP','v4':'k4*C3_star*IAP','v5':'k5*C8_star','v6':'k6*C3_star','v7':'k7*C3_star_IAP','v8':'k8*IAP-k_8','v9':'k9*C8-k_9','v10':'k10*C3-k_10','v11':'k11*C8_star*BAR-k_11*C8_star_BAR','v12':'k12*BAR-k_12','v13':'k13*C8_star_BAR'})
    net.set({'k1':5.8e-5,'k2':1e-5,'k3':5e-4,'k4':3e-4,'k5':5.8e-3,'k6':5.8e-3,'k7':1.73e-2,'k8':1.16e-2,'k9':3.9e-3,'k10':3.9e-3,'k11':5e-4,'k12':1e-3,'k13':1.16e-2,'k_1':0,'k_2':0,'k_3':0.21,'k_4':0,'k_5':0,'k_6':0,'k_7':0,'k_8':464.,'k_9':507.,'k_10':81.9,'k_11':0.21,'k_12':40.,'k_13':0})
    return net




def makesimplenet():
    """
    the chain -> A -> B -> with constant / linear reaction rates
    and a variable C = B
    """
    stoich = np.asarray([[1,-1,0],[0,1,-1]])
    species = ['A','B']
    reactions = ['v1','v2','v3']
    parameters = ['k2','k3']
    variables = ['C']
    net = brn.Brn(stoich,species,reactions,parameters,variables)
    net.set({'v1':1, 'v2':'k2*A', 'v3':'k3*C', 'C':'B', 'k2':1, 'k3':1})
    return net

def makescaledsimplenet():
    """
    the chain -> A -> B -> with constant / linear reaction rates
    and a variable C = B
    """
    stoich = np.asarray([[1,-1,0],[0,1,-1]])
    species = ['A','B']
    reactions = ['v1','v2','v3']
    parameters = ['k2','k3']
    variables = ['C']
    net = brn.Brn(stoich,species,reactions,parameters,variables,scale=2.)
    net.set({'v1':1, 'v2':'k2*A', 'v3':'k3*C', 'C':'B', 'k2':1, 'k3':1})
    return net

def make_conservation_net():
    """
    the network x0 -> x1 -> x2 with linear reaction rates
    """
    stoich = [[-1,0],[1,-1],[0,1]]
    parameters = ['k1','k2']
    #net = brn.Brn(stoich,parameters=parameters)
    net = brn.Brn(stoich,parameters=parameters,scale=2.)
    net.set({'v0':'k1*x0', 'v1':'k2*x1', 'k1':1, 'k2':1})
    return net

def make_nonlinear_scalar_net():
    """
    xdot = -k + x**2
    """
    stoich = [[-1, 1]]
    parameters = ['k']
    net = brn.Brn(stoich,['x'],['v1','v2'],parameters=parameters)
    net.set({'k':1., 'v1':'k', 'v2':'pow(x,2)'})
    return net

def make_scaled_net():
    """
    a simple network where the species are scaled differentially
    x1dot = (1/k1) k1 x1
    x2dot = (1/k2) k2 x2
    """
    stoich = [[1, 0],[0, 1]]
    species = ['x1','x2']
    reactions = ['v1','v2']
    parameters = ['k1','k2']
    net = brn.Brn(stoich,species,reactions,parameters,scale=['1/k1','1/k2'])
    net.set({'v1':'-k1*x1', 'v2':'-k2*x2', 'k1':1., 'k2':2.})
    return net

def adaptnet():
    """
    not yet complete !!!
    network with adaptation
    -> a ->
    -> m ->
    """
    stoich = np.asarray([[1,-1,0,0],[0,0,1,-1]])
    species = ['a','m']
    reactions = ['v1','v2','v3','v4']
    parameters = ['u','k2','k3','k4']
    net = brn.Brn(stoich,species,reactions,parameters,[])

def scalarnet():
    """
    xdot = (-1) * x
    """
    return brn.Brn([[-1]],['x'],['v'],rates={'v':'x'})

def make_mathexpr_net():
    """
    a test network with expressions from math module
    """
    rates = {'v1':'exp(-x1)', 'v2': 'x1**(2*k1)'}
    ics = {'x1':'sin(k1)'}
    parvals = {'k1':1.}
    return brn.Brn([[1,-1]],['x1'],['v1','v2'],['k1'],rates=rates,ics=ics,parvals=parvals)
