#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:53:54 2022

@author: fier887
"""
import process
import numpy as np

# output_dir = '/pic/projects/sooty2/fierce/output_LagrangeDroplets/'
output_dir = '/Users/fier887/OneDrive - PNNL/Documents/students/FY22_Jesse-Andersen/output_LagrangeDroplets/'
parcel_trace_dir = output_dir + 'parcel_traces/'
# process.plot_spatial_trajectory(parcel_trace_dir=parcel_trace_dir)

# parcel_nums = np.array([1,2,4])
parcel_nums = np.array([18,1,4])
process.plot_some_trajectories(parcel_nums,output_dir,avg_SS=0.,Slims=[-6,6],Dlims=[0.,15.],add_spatial_trajectories=False,plot_dir = 'figures/')
avg_SS = np.linspace(-0.04,0.03,15)
process.plot_activated_fraction_mean_SS(avg_SS, output_dir,plot_dir = output_dir+'figures/')