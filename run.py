#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 19:08:29 2022

@author: Laura Fierce

"""

import microphysics
import pickle, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from pyrcel.thermo import kohler_crit
from scipy.signal import savgol_filter
import warnings


def make_directories(avg_SS,parcel_trace_dir='parcel_traces/',all_particle_trace_dir='particle_traces/',all_particle_trace_dir_uniformS='particle_traces_uniformS/',plot_dir='figures/'):
    
    if not os.path.exists(parcel_trace_dir):
        os.mkdir(parcel_trace_dir)
    
    if not os.path.exists(all_particle_trace_dir):
        os.mkdir(all_particle_trace_dir)
        
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    
    
    for SS in avg_SS:
        particle_trace_dir = all_particle_trace_dir + '/avgSS_' + str(int(SS*1000)).zfill(6) +'/'
        if not os.path.exists(particle_trace_dir):
            os.mkdir(particle_trace_dir)
            
        particle_trace_dir_uniformS = all_particle_trace_dir_uniformS + '/avgSS_' + str(int(SS*1000)).zfill(6) +'/'
        
        if not os.path.exists(particle_trace_dir_uniformS):
            os.mkdir(particle_trace_dir_uniformS)    
    


def main(avg_SS, dry_diameter, kappa, density, N, dt, parcel_trace_dir='parcel_traces/',all_particle_trace_dir = 'particle_traces/', constant_SS = False, odesystem_type='lny', ignore_Tp=True, odesolver='BDF'):
    
    print('\n')
    print('Making particle traces...')
    print('\n')
    
    all_particle_traces = []
    for aa, SS in enumerate(avg_SS):
        print('SS offet = ', SS)
        particle_traces = many_trajectories(
                dry_diameter, kappa, density, N,
                odesystem_type=odesystem_type,
                avg_SS = SS,
                P0=101325., accom=0.3, dt=dt,
                parcel_trace_dir = parcel_trace_dir, 
                ignore_Tp=ignore_Tp, smooth=False,
                all_particle_trace_dir=all_particle_trace_dir,
                run_as_avg = constant_SS, odesolver=odesolver)
        all_particle_traces.append(particle_traces)
        print('\n')
        
    return all_particle_traces

def split_to_nodes(
        avg_SS,number_of_parcels,stop_time=60.,dt=0.5,
        parcel_trace_dir = 'parcel_traces/', all_particle_trace_dir='particle_traces/',
        all_particle_trace_dir_uniformS = 'particle_traces_uniformS/',
        N=1e6,dry_diameter=0.130e-6,kappa=1.0,density=2160.,odesystem_type='lny'):
#    output_dir = '/Users/fier887/OneDrive - PNNL/Documents/shared_files/CloudMC_output/library_03/'
    for ii,avg_SS_val in enumerate(avg_SS):
        main_lines = [
            'import parcels',
            'import run',
            'import process',
            'import numpy as np',
            'import time',
            'number_of_parcels = ' + str(number_of_parcels),
            'stop_time =' + str(stop_time),
            'dt = ' + str(dt), 
            'avg_SS = np.array([' + str(avg_SS_val) + '])', #np.linspace(-0.04,0.05,19)
            'N = ' + str(N),
            'odesystem_type = \'' + odesystem_type + '\'',
            'dry_diameter = ' + str(0.130e-6),
            'kappa = ' + str(kappa),
            'density = ' + str(density),
            'parcel_trace_dir = ' + '\'' +  parcel_trace_dir + '\'',
            'all_particle_trace_dir = ' + '\'' + all_particle_trace_dir + '\'',
            'all_particle_trace_dir_uniformS = ' + '\'' + all_particle_trace_dir_uniformS + '\'',
            
            # 'all_particle_traces = run.main(',
            # '    avg_SS, dry_diameter, kappa, density, N, dt,',
            # '    parcel_trace_dir=parcel_trace_dir,',
            # '    all_particle_trace_dir = all_particle_trace_dir,',
            # '    constant_SS=False, odesystem_type=odesystem_type, ignore_Tp=True)', 
            
            'all_particle_traces_uniform = run.main(',
            '    avg_SS, dry_diameter, kappa, density, N, dt,',
            '    parcel_trace_dir=parcel_trace_dir,',
            '    all_particle_trace_dir = all_particle_trace_dir_uniformS,',
            '    constant_SS=True, odesystem_type=odesystem_type, ignore_Tp=True)']
        
        main_filename = 'main_node' + str(ii).zfill(4) + '.py'
        
        with open(main_filename, 'w') as f:
            for line in main_lines:
                f.write(line)
                f.write('\n')
        
        run_lines = [
                '#!/bin/tcsh',
                '#SBATCH -A sooty2',
                '#SBATCH -p slurm',
                '#SBATCH -t 2:00:00',
                '#SBATCH -N 1',
                '#SBATCH -o node' + str(ii).zfill(4) + '.out',
                '#SBATCH -e node' + str(ii).zfill(4) + '.err',
                '#SBATCH -J node' + str(ii).zfill(4) + '.exe.default',
                # 'source /share/apps/python/anaconda3.6/etc/profile.d/conda.sh',
                # 'conda activate py37_pyrcel2',                
                'module load gcc/5.2.0',
                'module load netcdf',
                'unset PYTHONHOME',
                'module load python/3.7.2',
                'python main_node' + str(ii).zfill(4) + '.py']
        
        run_filename = 'run_node' + str(ii).zfill(4) + '.sh'
        with open(run_filename, 'w') as f:
            for line in run_lines:
                f.write(line)
                f.write('\n')
        os.system('chmod +x ' + run_filename)
        os.system('sbatch ' + run_filename)
        

def many_trajectories(dry_diameter, kappa, density, N,
                       avg_SS = 0.,
                       P0 = 101325., accom=0.3, part_ii=None,
                       dt = 1., ignore_Tp=False, smooth=True,
                       parcel_trace_dir = 'parcel_traces/',
                       all_particle_trace_dir = 'particle_traces/', 
                       odesystem_type='linear',
                       run_as_avg = False,
                       add_ignore_Tp=True,
                       odesolver='BDF'):
    
    particle_trace_dir = all_particle_trace_dir + '/avgSS_' + str(int(avg_SS*1000)).zfill(6) +'/'
    parcel_trace_nums = np.sort(np.array([int(onedir[-10:-4]) for onedir in os.listdir(parcel_trace_dir) if onedir.startswith('parcel_traces')]))

    particle_traces = {}
    particle_trace = {'particle_props':{},'traces':{}}
    print('parcel number . . . ')
    print(parcel_trace_nums)
    
    for parcel_num in parcel_trace_nums:
        
        if run_as_avg and parcel_num != min(parcel_trace_nums):
            particle_traces[parcel_num] = particle_trace
        else:
            parcel_filename = get_parcel_filename(parcel_num, parcel_trace_dir = parcel_trace_dir)
            particle_filename = get_particle_filename(parcel_num, particle_trace_dir)
            one_particle_trace = one_particle(
                    dry_diameter, kappa, density, N, parcel_filename, avg_SS=avg_SS, P0=P0, accom=accom, dt=dt, ignore_Tp=ignore_Tp, odesystem_type=odesystem_type, smooth=smooth, run_as_avg = run_as_avg, odesolver=odesolver)
    
            particle_trace['traces'] = one_particle_trace
            
            if add_ignore_Tp:
                one_particle_trace2 = one_particle(
                        dry_diameter, kappa, density, N, parcel_filename, avg_SS=avg_SS, P0=P0, accom=accom, dt=dt, ignore_Tp=True, odesystem_type=odesystem_type, smooth=smooth, run_as_avg = run_as_avg)
                particle_trace['traces_ignoreTp'] = one_particle_trace2
            
            T0 = np.mean(one_particle_trace['T'].ravel())
            r_crit, SS_crit = kohler_crit(T0, dry_diameter/2., kappa)
            particle_trace['props'] = {'Ddry':dry_diameter, 'kappa':kappa, 'density':density, 'N':N, 'Dp_crit':r_crit*2., 'SS_crit':SS_crit}
            f = open(particle_filename,'wb')
            pickle.dump(particle_trace,f)
            f.close()
            particle_traces[parcel_num] = particle_trace
        print(parcel_num)
        
    return particle_traces


def one_particle(dry_diameter, kappa, density, N, parcel_filename, avg_SS=0., P0=101325., accom=0.3, dt=1., ignore_Tp=True, odesystem_type='linear', smooth=True, run_as_avg=False, tau0=False, odesolver='BDF'):
    
    if run_as_avg:
        mixing_ratio_fun, temp_fun, SS_fun, t_end = environmental_trace_funs_constant_SS(parcel_filename, P0=P0, avg_SS=avg_SS, smooth=smooth)
    else:
        mixing_ratio_fun, temp_fun, SS_fun, t_end = environmental_trace_funs(parcel_filename, P0=P0, avg_SS=avg_SS, smooth=smooth)
    
    
    temp_0 = temp_fun(0.) # intitial temperature
    SS_0 = SS_fun(0.) # initial supersaturation
    if SS_0 > 0:
        SS_0 = 0.
    diameter_0, mixing_ratio_0, water_content_0 = microphysics.equilibrate_h2o(np.array([dry_diameter]),np.array([kappa]),np.array([N]),temp_0,SS_0,P0)
    if SS_0 > 0:
        r_crit, S_crit = kohler_crit(temp_0, dry_diameter/2., kappa)
        diameter_0 = r_crit*2.
    diameter_0 = diameter_0[0]
    
    if odesystem_type == 'linear':
        rhs_fcn = lambda t, y: microphysics.one_particle(y, t, mixing_ratio_fun(t), SS_fun(t), temp_fun(t), P0, dry_diameter, kappa, N, accom, density, ignore_Tp=ignore_Tp)
        q_0 = np.array([diameter_0, temp_0])
    elif odesystem_type == 'lny':
        rhs_fcn = lambda t, lny: microphysics.one_particle_lny(lny, t, mixing_ratio_fun(t), SS_fun(t), temp_fun(t), P0, dry_diameter, kappa, N, accom, density, ignore_Tp=ignore_Tp)            
        q_0 = np.log(np.array([diameter_0, temp_0]))
    elif odesystem_type == 'lnDp':
        rhs_fcn = lambda t, lnDp_and_T: microphysics.one_particle_lnDp(lnDp_and_T, t, mixing_ratio_fun(t), SS_fun(t), temp_fun(t), P0, dry_diameter, kappa, N, accom, density, ignore_Tp=ignore_Tp)
        q_0 = np.array([np.log(diameter_0), temp_0])
    
    t_span = [0.,t_end] # 
    t_eval = np.linspace(0.,t_end,int(t_end/dt+1))
    
    try:
        soln = solve_ivp(rhs_fcn, t_span, q_0, method=odesolver, max_step=1E-2, t_eval=t_eval)
        one_particle_trace = {
            't':soln['t'],'Dp':soln['y'][0],'Tp':soln['y'][1],'mixing_ratio':mixing_ratio_fun(soln['t']),'T':temp_fun(soln['t']),'SS':SS_fun(soln['t'])}        
    except:
        try:
            soln = solve_ivp(rhs_fcn, t_span, q_0, method=odesolver,t_eval=t_eval,max_step=1E-3)#, max_step=0.05)#, t_eval = t_eval)
            one_particle_trace = {
                't':soln['t'],'Dp':soln['y'][0],'Tp':soln['y'][1],'mixing_ratio':mixing_ratio_fun(soln['t']),'T':temp_fun(soln['t']),'SS':SS_fun(soln['t'])}    
        except:
            try:
                soln = solve_ivp(rhs_fcn, t_span, q_0, method=odesolver,t_eval=t_eval, rtol=1e-10, max_step=1E-3)#, max_step=0.05)#, t_eval = t_eval)
                one_particle_trace = {
                    't':soln['t'],'Dp':soln['y'][0],'Tp':soln['y'][1],'mixing_ratio':mixing_ratio_fun(soln['t']),'T':temp_fun(soln['t']),'SS':SS_fun(soln['t'])}
            except:
                one_particle_trace = {
                    't':np.array([0.]),'Dp':np.array([diameter_0]),'Tp':np.array([temp_0]),'mixing_ratio':np.array([mixing_ratio_0]),'T':np.array([temp_0]),'SS':np.array([SS_0])}
    
    if (len(soln['y'][0])!=len(t_eval)):
        warnings.warn("ERROR IN PARTICLE TRACE ODEs.")
        sys.exit()

    return one_particle_trace




def unravel_particle_traces(particle_traces):
    ts = np.array([0.])
    for ii in particle_traces.keys():
        ts_new = particle_traces[ii]['traces']['t']
        if len(ts_new)>len(ts):
            ts = ts_new
    print(ts)
        
    Dp = np.zeros([len(particle_traces),len(ts)])
    Tp = np.zeros([len(particle_traces),len(ts)])    
    SS = np.zeros([len(particle_traces),len(ts)])
    mixing_ratio = np.zeros([len(particle_traces),len(ts)])    
    T = np.zeros([len(particle_traces),len(ts)])
    Ddry = np.zeros(len(particle_traces))
    kappa = np.zeros(len(particle_traces))
    SS_crit = np.zeros(len(particle_traces))
    Dp_crit = np.zeros(len(particle_traces))
    jj = 0    
    for ii in particle_traces.keys():
        if True:#len(particle_traces[ii]['traces']['Dp'])==Dp.shape[1]:
            N_ts = len(particle_traces[ii]['traces']['t'])
            Dp[jj,:N_ts] = particle_traces[ii]['traces']['Dp']
            Tp[jj,:N_ts] = particle_traces[ii]['traces']['Tp']            
            T[jj,:N_ts] = particle_traces[ii]['traces']['T']
            SS[jj,:N_ts] = particle_traces[ii]['traces']['SS']
            mixing_ratio[jj,:N_ts] = particle_traces[ii]['traces']['mixing_ratio']            
            Ddry[jj] = particle_traces[ii]['props']['Ddry']
            kappa[jj] = particle_traces[ii]['props']['kappa']
            SS_crit[jj] = particle_traces[ii]['props']['SS_crit']
            Dp_crit[jj] = particle_traces[ii]['props']['Dp_crit'] 
            jj+=1
    return ts, Dp, Tp, T, mixing_ratio, SS, SS_crit, Dp_crit, Ddry, kappa

def unravel_particle_trace(particle_trace):
    ts = np.array([0.])
    ts_new = particle_trace['traces']['t']
    if len(ts_new)>len(ts):
            ts = ts_new
        
    Dp = np.zeros([len(ts)])
    Tp = np.zeros([len(ts)])
    S = np.zeros([len(ts)])
    wv = np.zeros([len(ts)])
    T = np.zeros([len(ts)])
    
    N_ts = len(particle_trace['traces']['t'])
    Dp[:N_ts] = particle_trace['traces']['Dp']
    Tp[:N_ts] = particle_trace['traces']['Tp']            
    T[:N_ts] = particle_trace['traces']['T']
    S[:N_ts] = particle_trace['traces']['SS']
    wv[:N_ts] = particle_trace['traces']['mixing_ratio']            
    Ddry = particle_trace['props']['Ddry']
    kappa = particle_trace['props']['kappa']
    S_crit = particle_trace['props']['SS_crit']
    D_crit = particle_trace['props']['Dp_crit']            
    return ts,Dp,Tp,T,wv,S,S_crit,D_crit,Ddry,kappa

#import matplotlib.pyplot as plt

def get_particle_filename(ii,particle_trace_dir):
    particle_filename = particle_trace_dir + 'particle_' + str(ii).zfill(9) + '.pkl'
    return particle_filename
   
def read_particle_traces(particle_trace_dir):
    trace_files = [onedir for onedir in os.listdir(particle_trace_dir) if not onedir.startswith('.')]
    particle_traces = {}
    for part_id_str in trace_files:
        ii = int(part_id_str[-9:-4])
        particle_traces[ii] = read_particle_trace(ii,particle_trace_dir)
    return particle_traces
    
def read_particle_trace(ii,particle_trace_dir):
    particle_filename = get_particle_filename(ii,particle_trace_dir)    
    f = open(particle_filename,'rb')
    particle_trace = pickle.load(f)
    f.close()
    return particle_trace

def get_parcel_filename(parcel_num,parcel_trace_dir = 'parcel_traces/'):
    parcel_filename = parcel_trace_dir + 'parcel_traces_' + str(parcel_num).zfill(6)+'.pkl'
    return parcel_filename
    
def read_parcel_trace(parcel_num,parcel_trace_dir = 'parcel_traces/'):
    parcel_filename = get_parcel_filename(parcel_num,parcel_trace_dir=parcel_trace_dir)#parcel_trace_dir + '/parcel_traces' + str(parcel_num).zfill(4)
    f = open(parcel_filename,'rb')
    trajectory = pickle.load(f)
    return trajectory

def environmental_trace_funs(parcel_filename, P0=101325., avg_SS=0., smooth=True):
    
    if smooth==True:
        fluctuating=False
    else:
        fluctuating=True
    
    f = open(parcel_filename,'rb')
    trajectory = pickle.load(f)
    f.close()
 
    time_trace = trajectory['time'] # seconds or minutes
    mixing_ratio_trace_orig = trajectory['mixing ratio']    #kg water per kg air
    temp_trace = trajectory['T']  # kelvin
    
    if smooth:
        mixing_ratio_trace_orig = savgol_filter(mixing_ratio_trace_orig, 151, 2)
        temp_trace = savgol_filter(temp_trace, 151, 2)
    if fluctuating:
        error = np.random.normal(np.zeros(len(temp_trace)), 0.01*temp_trace) # 10% error?
        temp_trace = temp_trace - error
        error = np.random.normal(np.zeros(len(mixing_ratio_trace_orig)), 0.01*mixing_ratio_trace_orig)
        mixing_ratio_trace_orig = mixing_ratio_trace_orig - error
    

    SS_trace = trajectory['SS'] + avg_SS - np.mean(trajectory['SS'])
    mixing_ratio_trace = microphysics.MixingRatio(SS_trace, temp_trace, P0)
    
    if any(mixing_ratio_trace < 0):
        print('ERROR: ', parcel_filename,' mixing ratio: ',mixing_ratio_trace,', SS_trace: ', SS_trace)
    
    mixing_ratio_fun = lambda t: np.interp(t, time_trace, mixing_ratio_trace)
    temp_fun = lambda t: np.interp(t, time_trace, temp_trace)
    SS_fun = lambda t: np.interp(t, time_trace, SS_trace)
    t_end = max(time_trace)
    
    return mixing_ratio_fun, temp_fun, SS_fun, t_end

def environmental_trace_funs_constant_SS(parcel_filename, P0=101325., avg_SS=0., smooth=True):
    
    f = open(parcel_filename,'rb')
    trajectory = pickle.load(f)
    f.close()
 
    time_trace = trajectory['time'] # seconds or minutes
    SS_trace = np.linspace(avg_SS, avg_SS, len(time_trace))
    
    mean_temp = np.mean(trajectory['T'])
    temp_trace = np.linspace(mean_temp, mean_temp, len(time_trace))
    
    mixing_ratio_trace = microphysics.MixingRatio(SS_trace, temp_trace, P0)  # kg water per kg air
    
    if any(mixing_ratio_trace < 0):
        print('ERROR: ', parcel_filename,' mixing ratio: ',mixing_ratio_trace,', SS_trace: ', SS_trace)
    
    mixing_ratio_fun = lambda t: np.interp(t, time_trace, mixing_ratio_trace)
    temp_fun = lambda t: np.interp(t, time_trace, temp_trace)
    SS_fun = lambda t: np.interp(t, time_trace, SS_trace)
    t_end = max(time_trace)
    
    return mixing_ratio_fun, temp_fun, SS_fun, t_end
