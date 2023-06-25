#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 10:08:26 2022

@author: Laura Fierce
"""

import parcels
import run

import numpy as np
import time

# ================ USER INPUTS ===================
number_of_parcels = 50000  # number of trajectories to make
num_start = 9400
stop_time = 60.0    # total time of trajectories, seconds
dt = 0.5    # time interval for evaluation of particle growth
avg_SS = np.linspace(-0.04,0.05,19) # np.linspace(-0.03, 0.03, 3) #supersaturation offset, will be added to supersaturation trajectory
N = 1e6 # number concentration, #/m^3
dry_diameter = 0.130e-6 # dry diameter, # m
kappa = 1.0 # 0.65 better? Hygroscopiscity

# ================================================
LES_dir ='LES/'
parcel_trace_dir = 'parcel_traces/'
all_particle_trace_dir = 'particle_traces/'
plot_dir = 'figs/'

start_time = time.time()
timesteps = np.linspace(90000,180000,3601)
run.make_directories(
    avg_SS,
    parcel_trace_dir=parcel_trace_dir,all_particle_trace_dir=all_particle_trace_dir,plot_dir=plot_dir)
parcels.trajectories(number_of_parcels, stop_time, num_start=num_start, LES_dir=LES_dir,parcel_trace_dir=parcel_trace_dir, timesteps=timesteps)
