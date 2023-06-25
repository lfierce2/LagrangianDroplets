#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 10:09:03 2022

@author: Laura Fierce
"""

import parcels
import run
import process
import numpy as np
import time

# ================ USER INPUTS ===================
number_of_parcels = 10000  # number of trajectories to make
stop_time = 60.0    # total time of trajectories, seconds
dt = 0.5    # time interval for evaluation of particle growth
avg_SS = np.linspace(-0.04,0.05,19)#np.linspace(-0.03, 0.03, 3) #supersaturation offset, will be added to supersaturation trajectory
N = 1e6 # number concentration, #/m^3
dry_diameter = 0.130e-6 # dry diameter, # m
kappa = 1.0 # 0.65 better? Hygroscopiscity
density = 2160. # aerosol density, kg/m^3

parcel_trace_dir = 'parcel_traces/'
all_particle_trace_dir = 'particle_traces/'
all_particle_trace_dir_uniformS = 'particle_traces_uniformS/'

run.split_to_nodes(
    avg_SS,number_of_parcels,stop_time=stop_time,dt=dt,
    parcel_trace_dir=parcel_trace_dir, 
    all_particle_trace_dir=all_particle_trace_dir,
    all_particle_trace_dir_uniformS=all_particle_trace_dir_uniformS,
    N=N,dry_diameter=dry_diameter,kappa=kappa,density=density,odesystem_type='lny')
