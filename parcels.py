#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 19:06:25 2022

@author:  Jesse Anderson and Laura Fierce

"""

import numpy as np
import pickle
import os
from scipy.integrate import solve_ivp
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator
import time, random, sys
import matplotlib.pyplot as plt

def trajectories(num_trajectories, stop_time=0, num_start=0,LES_dir='OUT_3D/',parcel_trace_dir='parcel_traces/',timesteps = range(90000,180000,3601)):
    
    time_grid = np.zeros(len(timesteps))
    for tt,timestep in enumerate(timesteps):
        filename = LES_dir+'PiChamber_huji_19K_trj_32_' + str(int(timestep)).zfill(10)+'.nc'
        
        LES_array = Dataset(filename)

        U_t = LES_array['U'][:]     # X Wind Component, m/s, dimensions=(t,z,y,x)
        V_t = LES_array['V'][:]     # Y Wind Component, m/s, dimensions=(t,z,y,x)
        W_t = LES_array['W'][:]     # Z Wind Component, m/s, dimensions=(t,z,y,x)
        T_t = LES_array['TABS'][:]      # Temperature, Kelvin, dimensions=(t,z,y,x)
        Qv_t = LES_array['QV'][:]/1000  # Water Vapor, kg water/kg air,  dimensions=(t,z,y,x)
        
        
        X_grid = LES_array['x'][:]
        Y_grid = LES_array['y'][:]
        Z_grid = LES_array['z'][:]
        time_grid[tt] = LES_array['time'][:]
        
        if tt == 0:
            U_array = U_t
            V_array = V_t
            W_array = W_t
            Tabs = T_t
            Qv = Qv_t   
        else:    
            U_array = np.concatenate((U_array,U_t), axis=0)  # Join the arrays long the time axis
            V_array = np.concatenate((V_array,V_t), axis=0)
            W_array = np.concatenate((W_array,W_t), axis=0)
            Tabs = np.concatenate((Tabs,T_t), axis=0)
            Qv = np.concatenate((Qv,Qv_t), axis=0)
        
    
    print(Qv.shape,X_grid.shape,Y_grid.shape,Z_grid.shape,time_grid.shape)
    # ## Function to interpolate from the virtual sensors. 
    interpolate_U = RegularGridInterpolator((X_grid, Y_grid, Z_grid, time_grid), U_array.T) # .T takes the transpose of the array, so that the dimensions match (x,y,z,t)
    interpolate_V = RegularGridInterpolator((X_grid, Y_grid, Z_grid, time_grid), V_array.T)
    interpolate_W = RegularGridInterpolator((X_grid, Y_grid, Z_grid, time_grid), W_array.T)
    interpolate_T = RegularGridInterpolator((X_grid, Y_grid, Z_grid, time_grid), Tabs.T)
    interpolate_Qv = RegularGridInterpolator((X_grid, Y_grid, Z_grid, time_grid), Qv.T)
    
    
    dt = (time_grid[1] - time_grid[0])/(num_trajectories-1)
    
    # # Solve the ODE, t_span is the timespan of the simulation
    for aa in range(num_start,num_trajectories):
        
        
        start_time = time.time() # start timer for this ensemble member
        
        # initial position of the parcel
        x0 = 2*random.random()
        y0 = 2*random.random()
        z0 = random.random()
        position = np.array([x0,y0,z0])
        
        # dPdt = change in parcel position per change in time (not change in pressure/change in time)
        dXdt = lambda t, position: find_Velocity(t, position, X_grid, Y_grid, Z_grid, interpolate_U, interpolate_V, interpolate_W)
        # t_span = [t_start, t_stop]
        t_start = dt*aa
        t_stop = time_grid[-1]
        t_span = [t_start,t_stop]
        sol = solve_ivp(dXdt, t_span, position)#, max_step=(0.1))
        
        ## Get the particle traces from the solver
        
        time_trace = np.array(list(sol.t))
        x_trace = np.array(list(sol.y[0]))
        y_trace = np.array(list(sol.y[1]))
        z_trace = np.array(list(sol.y[2]))
        temperature_trace = np.zeros((len(time_trace)))
        mixing_ratio_trace = np.zeros((len(time_trace)))
        
        index=0
        for one_time in time_trace:
            ##  Interpolate r and T along the trace  
            #t_trace will be a function of n when the LES field is changing
            temperature_trace[index] = interpolate_T([x_trace[index], y_trace[index], z_trace[index], one_time])
            mixing_ratio_trace[index] = interpolate_Qv([x_trace[index], y_trace[index], z_trace[index], one_time])
            index += 1
                    
        SuperSat_trace = SuperSaturation(mixing_ratio_trace, temperature_trace)
        print('S_mean',np.mean(SuperSat_trace))
        elapsed_time = time.time() - start_time
        print('Trace Number '+str(aa+1)+', elapsed time:', elapsed_time)
        
        trajectory = {}
        trajectory['x'] = x_trace
        trajectory['y'] = y_trace
        trajectory['z'] = z_trace
        trajectory['time'] = time_trace
        trajectory['mixing ratio'] = mixing_ratio_trace
        trajectory['T'] = temperature_trace
        trajectory['SS'] = SuperSat_trace
        # path = os.getcwd()+'/parcel_traces'
        file_name = parcel_trace_dir+'parcel_traces_' + str(aa).zfill(6)+'.pkl'
        pickle.dump(trajectory, open(file_name,'wb'))
        
    return


# ## Function that handles periodic boundary conditions
def RBC(position, X_grid, Y_grid, Z_grid):
    x = position[0]
    y = position[1]
    z = position[2]
    
    if x < np.min(X_grid):
        d = abs(x-np.min(X_grid))
        position[0] = np.min(X_grid)+d
    elif x > np.max(X_grid):
        d = abs(x-np.max(X_grid))
        position[0] = np.max(X_grid)-d
    if y < np.min(Y_grid):
        d = abs(y-np.min(Y_grid))
        position[1] = np.min(Y_grid)+d
    elif y > np.max(Y_grid):
        d = abs(y-np.max(Y_grid))
        position[1] = np.max(Y_grid)-d
    if z < np.min(Z_grid):
        d = abs(z-np.min(Z_grid))
        position[2] = np.min(Z_grid)+d
    elif z > np.max(Z_grid):
        d = abs(z-np.max(Z_grid))
        position[2] = np.max(Z_grid)-d
        
    return position
    


# ## Function that interpolates the velocity vectors 
def find_Velocity(t, position, X_grid, Y_grid, Z_grid, interpolate_U, interpolate_V, interpolate_W): 
    
    ## Reflective boundary conditions, keep the parcel in the box
    position = RBC(position, X_grid, Y_grid, Z_grid)
    
    x = position[0]
    y = position[1]
    z = position[2]

    dxdt = interpolate_U([x, y, z, t])  # interpolates the x wind component at x, y, z, t
    dydt = interpolate_V([x, y, z, t])  # interpolates the y wind component at x, y, z, t
    dzdt = interpolate_W([x, y, z, t])  # interpolates the z wind component at x, y, z, t
    velocity = np.array([dxdt, dydt, dzdt]) # velocity of parcel
    velocity = np.reshape(velocity,3)
    return velocity


# # Define the Calculation of the supersatration, (Magnus Approximation)
def SuperSaturation(mixing_ratio,temperature):
    P_atm = 1000 # Atmospheric pressure (hPa)
    T = temperature-273.15 #Convert to degrees C
    #Saturation vapor pressure (T)
    saturation_pressure = 6.1094*np.exp(17.625*T/(243.04+T)) #(hPa)
    #Actual Vapor Pressure
    vapor_pressure = (mixing_ratio)*P_atm/(0.622+(mixing_ratio)) # (hPa)
    #Supersaturation
    SS = vapor_pressure/saturation_pressure-1
    return SS
