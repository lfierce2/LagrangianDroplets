#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 09:37:27 2022

@author: fier887
"""
import matplotlib.pyplot as plt
import pickle
import numpy as np
import run
import parcels
from netCDF4 import Dataset
from scipy.optimize import fsolve
import scipy.constants as c
import random


# def timeseries():
#     fig, ax = plt.subplots(3)
#     fig.set_size_inches(3.5,6)
#     return hfigs, haxs1, haxs2
M_dry_air = 28.9647/1000. # kg/mol
M_h2o = 18./1000.
M_co2 = 44./1000.
M_nacl = 58.44/1000.



def plot_some_trajectories(parcel_nums,output_dir,avg_SS=0.,Slims=[-10,10],Dlims=[0.,25.],add_spatial_trajectories=False,plot_dir = 'figures/'):
    fig, axs1 = plt.subplots(len(parcel_nums))
    axs2 = np.array([])
    # axs3 = np.array([])
    fig.set_size_inches(3.5,len(parcel_nums)*1.8)
    ylab_idx = np.floor(len(axs1)/2)
    particle_trace_dir= output_dir + 'particle_traces/avgSS_' + str(int(avg_SS*100.)).zfill(6) + '/'
    parcel_trace_dir= output_dir + 'parcel_traces/'
    for pp,(parcel_num,ax1) in enumerate(zip(parcel_nums,axs1)):
        particle_trace = pickle.load(open(particle_trace_dir + 'particle_' + str(parcel_num).zfill(9) + '.pkl','rb'))
        ts,lnDp,Tp,T,wv,S,S_crit,D_crit,Ddry,kappa = run.unravel_particle_trace(particle_trace)  
        Dp=np.exp(lnDp)
        
        # fig, ax1 = plt.subplots(figsize=(4.5,2.5))
        ax2 = ax1.twinx()
        np.hstack([axs2,ax2])
        
        hln1, = ax1.plot(ts,S*100,linewidth=2.)
        hln2, = ax2.plot(ts,Dp*1e6,linewidth=2.,color='C3')
        
        if pp == ylab_idx:
            ax1.set_ylabel('supersaturation [%]')
            ax2.set_ylabel('droplet diameter [$\mu$m]',rotation=270,labelpad=18)
        if pp < (len(parcel_nums) - 1):
            ax1.set_xticklabels('')
            ax2.set_xticklabels('')    
        else:
            ax1.set_xlabel('time [s]')
            
        ax1.set_xlim([min(ts),max(ts)])
        ax2.set_xlim([min(ts),max(ts)])    
        
        ax1.yaxis.label.set_color(hln1.get_color())
        ax2.yaxis.label.set_color(hln2.get_color())
        
        ax1.spines["right"].set_edgecolor(hln1.get_color())
        ax2.spines["right"].set_edgecolor(hln2.get_color())
        
        ax1.tick_params(axis='y', colors=hln1.get_color())
        ax2.tick_params(axis='y', colors=hln2.get_color())   
        
        ax2.hlines(D_crit*1e6,0.,max(ts),linestyle=':',color=hln2.get_color())
        ax1.hlines(S_crit*100.,0.,max(ts),linestyle=':',color=hln1.get_color())
        print(S_crit*100.)
        ax1.set_ylim(Slims)
        ax2.set_ylim(Dlims)
        if add_spatial_trajectories:
            big_ax_position = ax1.get_position()
            dy = big_ax_position.y1 - big_ax_position.y0
            dx = big_ax_position.x1 - big_ax_position.x0
            little_ax_position = [big_ax_position.x0-dx*0.6, big_ax_position.y0 + dy*0.4, dx/3.,dy/1.8]
            
            ax3 = fig.add_axes(little_ax_position,projection='3d')
            ax3.set_xticks([],labels='')
            ax3.set_yticks([],labels='')
            ax3.set_zticks([],labels='')
            ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax3.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax3.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax3.yaxis.set_color((1.0, 1.0, 1.0, 0.0))
            # ax3.zaxis.set_color((1.0, 1.0, 1.0, 0.0))

            # ax3.patch.set_facecolor('none')
            # ax3.patch.set_edgecolor('none')            
            
            ax3 = plot_spatial_trajectory(parcel_num,parcel_trace_dir,ax3)
            # ax3.zaxis.set_visible(False)
            
            # ax3.patch.set_facecolor('white')
            # ax3.patch.set_alpha(1.0)
            # print(little_ax_position)
            
#    ax2.set_yscale('log')
#    ax.yaxis.grid(color='lightgray', linestyle='dashed')
    plt.tight_layout()
    fig.savefig(plot_dir + 'trajectories.png',dpi=1000)
    # plt.show()
    if add_spatial_trajectories:
        return fig,ax1,ax2,ax3
    else:
        return fig,ax1,ax2

def plot_spatial_trajectory(parcel_num,parcel_trace_dir,ax):
    # print('parcel_num',parcel_num)
#for part_id in particle_traces.keys():
#    parcel_num = particle_traces[part_id]['trajectory_id']
    trajectory = run.read_parcel_trace(parcel_num,parcel_trace_dir = parcel_trace_dir)
#    particle_trace = run.read_particle_trace(ii,particle_trace_dir)
    # Data for a three-dimensional line
    xdata = trajectory['x']
    ydata = trajectory['y']
    zdata = trajectory['z']
    # fig = plt.figure(figsize=[2.5,2.5])
    # ax = plt.axes(projection='3d')    
    ax.plot3D(xdata, ydata, zdata, 'gray')
    # if not cvar == None:
    #     cdata = trajectory[cvar]
    #     ax.scatter3D(xdata, ydata, zdata, c=cdata);
    
#    plt.title('particle '+str(part_id)+ ', trajectory '+ str(parcel_num))
#    plt.title('trajectory '+ str(parcel_num))
    ax.set_xlim([0.,2.])
    ax.set_ylim([0.,2.])        
    ax.set_zlim([0.,1.])
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Z [m]')
    # plt.show()
    # filename = plot_dir + '/activated_fraction_mean_SS.png'
    # plt.tight_layout()
    # fig.savefig(filename,dpi=1000)
    
    # return fig,ax

def plot_activated_fraction_mean_SS(avg_SS, output_dir,plot_dir = 'figures/'):
    all_particle_trace_dir = output_dir+'/particle_traces'
    
    N_tot = []
    N_ccn = []
    N_ccn_kohler = []
    N_ccn_tau0 = []
    mean_SS = []
    
    
    for aa, SS in enumerate(avg_SS):
        SS_thresh = SS
        particle_trace_dir = all_particle_trace_dir + '/avgSS_' + str(int(SS*1000)).zfill(6) +'/'
        particle_traces = run.read_particle_traces(particle_trace_dir)
#        particle_traces = all_particle_traces[aa]
        time, Dp, Tp, T, mixing_ratio, SS_env, SS_crit, Dp_crit, Ddry, kappa = run.unravel_particle_traces(particle_traces)
        Dp=np.exp(Dp)
        Ntot_t = np.zeros([len(time)])
        Nccn_t = np.zeros([len(time)])
        Nccn_t_tau0 = np.zeros([len(time)])
        Nccn_kohler_t = np.zeros([len(time)])
        mean_SS.append(np.mean(SS))
        
        
        # for tt in range(len(time)): 
        #     Ntot_t[tt] = len(Dp[:,tt])#sum(Dp[:,tt] < 1.)
        #     Nccn_t[tt] = sum(Dp[:,tt]>=Dp_crit) #(Dp[:,tt] <1.) & (Dp[:,tt] >= Dp_crit))
        #     Nccn_kohler_t[tt] = sum(SS_crit <= SS_thresh)#sum((Dp[:,tt] < 1.) & (SS_crit <= SS_thresh))
        #     Nccn_t_tau0[tt] = sum(SS_crit <= (SS_env[:,tt]+SS))
        # N_tot.append(Ntot_t)
        # N_ccn.append(Nccn_t)
        # N_ccn_kohler.append(Nccn_kohler_t)
        # N_ccn_tau0.append(Nccn_t_tau0)
        
        N0 = 1.
        Ns = np.zeros(Dp.shape)
        wps = np.zeros(Dp.shape)
        for ii in range(Dp.shape[0]):
            N_vs_t,wp_vs_t = N_vs_t__settling(N0,time,Dp[ii,:],Tp[ii,:],SS_env[ii,:]+1.,rho_p=1000.,ht=1.)
            Ns[ii,:] = N_vs_t
            wps[ii,:] = wp_vs_t
            
        Dp_bins = np.logspace(-7,-4,100);
        plt.hist(Dp[:,-1],bins=Dp_bins,weights=Ns[:,-1]); plt.xscale('log'); plt.title(SS); plt.show();
        print(aa,SS)
        for tt in range(len(time)): 
            Ntot_t[tt] = sum(Ns[:,tt]) #len(Dp[:,tt])#sum(Dp[:,tt] < 1.)
            Nccn_t[tt] = sum(Ns[:,tt]*(Dp[:,tt]>=Dp_crit)) #sum(Dp[:,tt]>=Dp_crit) #(Dp[:,tt] <1.) & (Dp[:,tt] >= Dp_crit))
            Nccn_kohler_t[tt] = sum(Ns[:,tt]*(SS_crit <= SS_thresh))#sum((Dp[:,tt] < 1.) & (SS_crit <= SS_thresh))
            Nccn_t_tau0[tt] = sum(Ns[:,tt]*(SS_crit <= (SS_env[:,tt]+SS)))
        N_tot.append(Ntot_t)
        N_ccn.append(Nccn_t)
        N_ccn_kohler.append(Nccn_kohler_t)
        N_ccn_tau0.append(Nccn_t_tau0)
                

    
    N_tot = np.array(N_tot)
    N_ccn = np.array(N_ccn)
    N_ccn_kohler = np.array(N_ccn_kohler)
    mean_SS = np.array(mean_SS)*100
    
    # SS_mids_avg0,s_pdf = get_s_pdf(avg_SS=0.,LES_dir='/Users/fier887/OneDrive - PNNL/Documents/shared_files/ICLASS/LES/OUT_3D/')
    # grid = plt.GridSpec(9,3,hspace=0.4)
    # fig = plt.figure()
    # ax = fig.add_subplot(grid[3:,:])
    # ax1 = fig.add_subplot(grid[0,:])
    # ax2 = fig.add_subplot(grid[1,:])
    # ax3 = fig.add_subplot(grid[2,:])
    # fig.set_size_inches(3.5,3.5)
    
    # ax1.plot(SS_mids_avg0-2.,s_pdf)
    # ax2.plot(SS_mids_avg0,s_pdf)
    # ax3.plot(SS_mids_avg0+2.,s_pdf)
    
    #plt.legend(loc='upper left')    
    fig, ax = plt.subplots(1,1)
    hln_uniform, = ax.plot(np.array([min(mean_SS),np.mean(SS_crit*100.),np.mean(SS_crit*100.),max(mean_SS)]),np.array([0.,0.,1.,1.]),color='C1',linewidth=3.)    
    hln_variable, = ax.plot(mean_SS,(N_ccn/N_tot)[:,len(time)-1],linewidth=2.,color='C2')
    hln_tau0, = ax.plot(mean_SS,(N_ccn_tau0/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color='C2')
    # ax.plot(mean_SS,(N_ccn_kohler/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color='C1')
    
    # ax.plot(mean_SS,(N_ccn_kohler/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color='C1')
    
    fig.set_size_inches(3.5,2.5)
    
    ax.set_xlabel('mean supersaturation (%)')
    ax.set_ylabel('activated fraction')
    
    ax.set_xlim([-4.,3.])
    ax.set_ylim([0.,1.001])
    ax.legend([hln_variable,hln_tau0,hln_uniform],['turbulent','$\\tau_{\mathrm{evap}}=0$','uniform'],fontsize=10,loc='upper left',handletextpad=0.5)
    # ax1.set_xlim([-4.,3.])
    # ax2.set_xlim([-4.,3.])
    # ax3.set_xlim([-4.,3.])
    
    filename = plot_dir + '/activated_fraction_mean_SS.png'
    plt.tight_layout()
    fig.savefig(filename,dpi=1000)
    
    return fig,ax,[hln_variable,hln_tau0,hln_uniform]

def get_s_pdf(avg_SS=0.,LES_dir='/Users/fier887/OneDrive - PNNL/Documents/shared_files/ICLASS/LES/OUT_3D/'):    
    SS_all = np.array([])
    for ii in range(0,3001,25):
        ## Load the variable from the LES output
        filename = LES_dir+'PiChamber_huji_19K_surfmod_trj_32_000009'+str(ii).zfill(4)+'.nc'
        LES_array = Dataset(filename)
        T_t = LES_array['TABS'][:]      # Temperature, Kelvin, dimensions=(t,z,y,x)
        Qv_t = LES_array['QV'][:]/1000  # Water Vapor, kg water/kg air,  dimensions=(t,z,y,x)
        SS_t = parcels.SuperSaturation(Qv_t, T_t)*100.
        SS_all = np.hstack([SS_all,SS_t.ravel()])
    SS_true_mean = np.mean(SS_all)
    N,SS_edges,patches = plt.hist(SS_all-SS_true_mean+avg_SS,100)
    # SS_mids = output[1]
    dSS = SS_edges[1] - SS_edges[0]
    SS_mids = SS_edges[:-1]+dSS/2.
    pdf = N/(sum(N)*dSS)
    return SS_mids,pdf
    
def get_nonoverlapping_xy(N_samples,r_closest=0.1):
    
    xs = np.array([np.random.uniform()])
    ys = np.array([np.random.uniform()])
    while len(xs) < N_samples:
        x = np.random.uniform()
        y = np.random.uniform()
        # r = np.hstack([xs/2.,ys/2.,(1-xs)/2.,(1-ys)/2.,np.sqrt((xs-x)**2 + (ys-y)**2)])
        r = np.hstack([np.sqrt((xs-x)**2 + (ys-y)**2)])
        if all(r>r_closest):
            xs = np.append(xs,x)
            ys = np.append(ys,y)        
    return xs, ys

def plot_activated_fraction_mean_SS__with_particles(avg_SS,output_dir,plot_dir = 'figures/',ss_plot = [-0.03,0],tt_plot=-1,N_particles=30,r_closest=0.1,droplets_above=True):
    # all_particle_trace_dir = output_dir+'/particle_traces'
    all_particle_trace_dir_uniform = output_dir+'/particle_traces_uniformS'
    
    all_particle_trace_dir = output_dir + '../output_LagrangeDroplets_old/particle_traces'
    
    fig = plt.figure(constrained_layout=True)
    # fig.set_size_inches(7.,2.5)
    fig.set_size_inches(3.5,5.)
    plt.show()
    nrows = 3
    if droplets_above:
        gs = fig.add_gridspec(nrows=nrows, ncols=len(ss_plot),height_ratios=[1,1,3],wspace=0.1)
    else:
        gs = fig.add_gridspec(nrows=nrows, ncols=len(ss_plot),height_ratios=[3,1,1])
    Dp_uniform = np.zeros(len(ss_plot))
    frac_aerosol_uniform = np.zeros(len(ss_plot))
    ax_particles = []
    for ii in range(2):
        ax_ii = []
        for jj in range(len(ss_plot)):
            if droplets_above:
                ax_ii.append(fig.add_subplot(gs[ii,jj]))
            else:
                ax_ii.append(fig.add_subplot(gs[nrows-2+ii,jj]))
        ax_particles.append(ax_ii)
    ax_particles = np.array(ax_particles)
    if droplets_above:
        ax = fig.add_subplot(gs[2:,:])
        # ax_leg = fig.add_subplot(gs[:2,-1])
    else:
        ax = fig.add_subplot(gs[:-2,:])
        # ax_leg = fig.add_subplot(gs[-2:,-1])



    
    N_tot = []
    N_ccn = []
    N_ccn_kohler = []
    N_ccn_tau0 = []
    mean_SS = []
    # Dp_bins = np.logspace(-7,-4,100);
    jj = 0
    col_inactive = np.array([106, 106, 106])/255.#np.array([212, 213, 214])/255.#np.array([181, 5, 2])/255.#
    edge_col_inactive = np.array([0.,0.,0.])#np.array([59, 60, 61])/255.
    col_active = np.array([92, 154, 255])/255.
    edge_col_active =np.array([50, 94, 166])/255.#np.array([0.,0.,0.])#
    
    multiplier = 1.8e6
    
    leg_strs = ['0.2~$\mu m$','1~$\mu m$','5~$\mu m$','CCN-active','CCN-inactive']
    marker_sizes = np.array([0.2e-7,1.e-6,5e-6,1e6,1e6])*multiplier
    
    marker_cols = np.zeros([5,3])
    marker_cols[:2,:] = 1.
    marker_cols[3,:] = col_active
    marker_cols[4,:] = col_inactive    
    
    marker_edgecols = np.zeros([5,3])    
    marker_cols[3,:] = edge_col_active
    marker_cols[4,:] = edge_col_inactive    
    
    linecol_uniform = 'C1'
    linecol_turbulent = 'C2'
    
    cols = np.zeros([N_particles,3])
    cols_uniform = np.zeros([N_particles,3])    
    edge_cols = np.zeros([N_particles,3])
    edge_cols_uniform = np.zeros([N_particles,3])    
    xs,ys = get_nonoverlapping_xy(N_particles,r_closest=r_closest);
    for aa, SS in enumerate(avg_SS):
        SS = np.round(SS,8)
        particle_trace_dir = all_particle_trace_dir + '/avgSS_' + str(int(avg_SS[aa]*1000)).zfill(6) +'/'
        particle_traces = run.read_particle_traces(particle_trace_dir)
        random_index = random.sample(range(len(particle_traces)),N_particles)
        these_particles = random_index[:N_particles]
        
        
        if SS in ss_plot:
            particle_trace_dir_uniform = all_particle_trace_dir_uniform + '/avgSS_' + str(int(avg_SS[aa]*1000)).zfill(6) +'/'
            particle_traces_uniform = run.read_particle_traces(particle_trace_dir_uniform)
            time, Dp_uniform_all, Tp, T, mixing_ratio, SS_env, SS_crit, Dp_crit_uniform, Ddry_uniform, kappa = run.unravel_particle_traces(particle_traces_uniform)
            if len(Dp_uniform_all) == 1:
                Dp_uniform_all=np.exp(Dp_uniform_all[0])
            else:
                Dp_uniform_all=np.exp(Dp_uniform_all)
            Dp_uniform[jj] = Dp_uniform_all[tt_plot]
            frac_aerosol_uniform[jj] = Ddry_uniform**3/Dp_uniform_all[tt_plot]**3
            # Dp_uniform[jj] = get_Dwet(particle_traces_uniform[0]['props']['Ddry'], particle_traces_uniform[0]['props']['kappa'], (100+SS)/100., np.mean(particle_traces_uniform[0]['traces']['T']))
            print('aa',aa,'SS',SS,'Dp_uniform[jj]',Dp_uniform[jj])
            # if SS == min(ss_plot):
            #     Dwet_low = get_Dwet(particle_traces_uniform[0]['props']['Ddry'], particle_traces_uniform[0]['props']['kappa'], (100+SS)/100., np.mean(particle_traces_uniform[0]['traces']['T']))
            #     print('aa',aa,'SS',SS,'Dwet_low',Dwet_low)
            # else:
            #     Dwet_crit = get_Dwet(particle_traces_uniform[0]['props']['Ddry'], particle_traces_uniform[0]['props']['kappa'], (100+SS)/100., np.mean(particle_traces_uniform[0]['traces']['T']))#particle_traces[0]['props']['Dp_crit']
            #     print('aa',aa,'SS',SS,'Dwet_crit',Dwet_crit)
#        particle_traces = all_particle_traces[aa]
        time, Dp, Tp, T, mixing_ratio, SS_env, SS_crit, Dp_crit, Ddry, kappa = run.unravel_particle_traces(particle_traces)
        Dp=np.exp(Dp)
        Ntot_t = np.zeros([len(time)])
        Nccn_t = np.zeros([len(time)])
        Nccn_t_tau0 = np.zeros([len(time)])
        Nccn_kohler_t = np.zeros([len(time)])
        mean_SS.append(np.mean(SS))
        
        for tt in range(len(time)): 
            Ntot_t[tt] = len(Dp[:,tt])#sum(Dp[:,tt] < 1.)
            Nccn_t[tt] = sum(Dp[:,tt]>=Dp_crit) #(Dp[:,tt] <1.) & (Dp[:,tt] >= Dp_crit))
            Nccn_kohler_t[tt] = sum(SS_crit <= SS)#sum((Dp[:,tt] < 1.) & (SS_crit <= SS_thresh))
            Nccn_t_tau0[tt] = sum(SS_crit <= (SS_env[:,tt]+SS))
        N_tot.append(Ntot_t)
        N_ccn.append(Nccn_t)
        N_ccn_kohler.append(Nccn_kohler_t)
        N_ccn_tau0.append(Nccn_t_tau0)
        
        # run.one_particle(dry_diameter, kappa, density, N, parcel_filename, 
                         # avg_SS=avg_SS, ignore_Tp=True, odesystem_type='linear', smooth=True, run_as_avg=False, tau0=False)
        # print(aa,SS)
        
        print('SS',SS,'SS in ss_plot',SS in ss_plot)
        if SS in ss_plot:
            
            for ii in range(N_particles):
                if Dp_uniform[jj]>=Dp_crit_uniform[0]:
                    cols_uniform[ii,:] = col_active
                    edge_cols_uniform[ii,:] = edge_col_active

                else:
                    cols_uniform[ii,:] = col_inactive
                    edge_cols_uniform[ii,:] = edge_col_inactive

            for ii in range(N_particles):
                if Dp[these_particles[ii],tt_plot]>=Dp_crit[these_particles[ii]]:
                    cols[ii,:] = col_active
                    edge_cols[ii,:] = edge_col_active
                else:
                    cols[ii,:] = col_inactive
                    edge_cols[ii,:] = edge_col_inactive
            hscat_uniform = ax_particles[0,jj].scatter(xs,ys,np.ones(N_particles)*Dp_uniform[jj]*multiplier,cols_uniform,edgecolors=edge_cols_uniform,linewidth=0.25)
            # if SS == min(ss_plot):
            #     ax_particles[0,jj].scatter(xs,ys,np.ones(N_particles)*Dwet_low*multiplier)
            # else:
            #     ax_particles[0,jj].scatter(xs,ys,np.ones(N_particles)*Dwet_crit*multiplier)
            ax_particles[1,jj].scatter(xs,ys,Dp[these_particles,tt_plot]*multiplier,cols,edgecolors=edge_cols,linewidth=0.2)
            print('aa',aa,'SS',SS,'mean Dp',np.mean(Dp[:,tt_plot]),'gm Dp',np.exp(np.mean(np.log(Dp[:,tt_plot]))))
            ax_particles[0,jj].set_xticks([])
            ax_particles[0,jj].set_yticks([])            
            ax_particles[1,jj].set_xticks([])
            ax_particles[1,jj].set_yticks([])
            
            ax_particles[0,jj].spines['left'].set_color(linecol_uniform)
            ax_particles[0,jj].spines['right'].set_color(linecol_uniform)
            ax_particles[0,jj].spines['top'].set_color(linecol_uniform)
            ax_particles[0,jj].spines['bottom'].set_color(linecol_uniform)
            
            ax_particles[1,jj].spines['left'].set_color(linecol_turbulent)
            ax_particles[1,jj].spines['right'].set_color(linecol_turbulent)
            ax_particles[1,jj].spines['top'].set_color(linecol_turbulent)
            ax_particles[1,jj].spines['bottom'].set_color(linecol_turbulent)     
            
            ax_particles[0,jj].set_title('$s=$' + str(int(SS*100))+'%',fontsize=10)
            
            ax_particles[0,jj].patch.set_facecolor('w')
            ax_particles[1,jj].patch.set_facecolor('w')
            

            # ax_particles[0,jj].set_title('r$\bar{s}=$' + str(int(SS*100))+'%',fontsize=10)
            jj += 1
            
            
            # ax_particles[0,jj].hist(Dp[:,-1],bins=Dp_bins); plt.xscale('log'); plt.title(SS);# plt.show();
            
            
    N_tot = np.array(N_tot)
    N_ccn = np.array(N_ccn)
    N_ccn_kohler = np.array(N_ccn_kohler)
    mean_SS = np.array(mean_SS)*100
    
    # SS_mids_avg0,s_pdf = get_s_pdf(avg_SS=0.,LES_dir='/Users/fier887/OneDrive - PNNL/Documents/shared_files/ICLASS/LES/OUT_3D/')
    # grid = plt.GridSpec(9,3,hspace=0.4)
    # fig = plt.figure()
    # ax = fig.add_subplot(grid[3:,:])
    # ax1 = fig.add_subplot(grid[0,:])
    # ax2 = fig.add_subplot(grid[1,:])
    # ax3 = fig.add_subplot(grid[2,:])
    # fig.set_size_inches(3.5,3.5)
    
    # ax1.plot(SS_mids_avg0-2.,s_pdf)
    # ax2.plot(SS_mids_avg0,s_pdf)
    # ax3.plot(SS_mids_avg0+2.,s_pdf)
    
    #plt.legend(loc='upper left')    
    # fig, ax = plt.subplots(1,1)
    hln_uniform, = ax.plot(np.array([min(mean_SS),np.mean(SS_crit*100.),np.mean(SS_crit*100.),max(mean_SS)]),np.array([0.,0.,1.,1.]),color=linecol_uniform,linewidth=3.)    
    hln_variable, = ax.plot(mean_SS,(N_ccn/N_tot)[:,len(time)-1],linewidth=2.,color=linecol_turbulent)
    hln_tau0, = ax.plot(mean_SS,(N_ccn_tau0/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color=linecol_turbulent)
    # ax.plot(mean_SS,(N_ccn_kohler/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color='C1')
    
    # ax.plot(mean_SS,(N_ccn_kohler/N_tot)[:,len(time)-1],linestyle='--',linewidth=2.,color='C1')
    
    
    ax.set_xlabel('mean supersaturation, $s$')#', ' + r'$\bar{s}$')
    ax.set_ylabel('activated fraction')
    
    ax.set_xlim([-4.,3.])
    ax.set_ylim([0.,1.001])
    
    
    ss_ticks = np.arange(-4.,3.)
    ax.set_xticks([int(ss_tick) for ss_tick in ss_ticks])
    ax.set_xticklabels([str(int(ss_tick))+'%' for ss_tick in ss_ticks])    
    # ax.set_xticklabels([str(int(s_percent)) + '%' for s_percent in ss_plot*100.])
    ax.legend(
        [hln_uniform,hln_variable,hln_tau0],
        ['uniform $s$','with turbulent\nfluctuations in $s$','with fluctuations\nbut $\\tau_{\mathrm{evap}}=0$'],
        frameon=False,fontsize=9,loc='upper left',handletextpad=0.5)
    # ax1.set_xlim([-4.,3.])
    # ax2.set_xlim([-4.,3.])
    # ax3.set_xlim([-4.,3.])
    
    # lgnd = ax_leg.legend([*hscat_uniform.legend_elements("sizes", num=3),*hscat_uniform.legend_elements("colors",num=2)],leg_strs,loc="upper left", fontsize=10)
    # lgnd = ax_leg.legend([hscat_uniform,hscat_uniform,hscat_uniform,hscat_uniform,hscat_uniform],leg_strs,loc="upper left", fontsize=10)
    
    # # lgnd = ax_leg.legend([hscat_uniform]*len(leg_strs),leg_strs,loc="upper left", fontsize=10)
    
    # for (handle,maker_size,marker_col,marker_edgecol) in zip([*hscat_uniform.legend_elements("sizes", num=3),*hscat_uniform.legend_elements("colors",num=2)],marker_sizes,marker_cols,marker_edgecols):
    #     print(handle,maker_size,marker_col,marker_edgecol)
    #     handle.set_sizes(maker_size)
    #     handle.set_markerfacecolor(marker_col)
    #     handle.set_markeredgecolor(marker_edgecol)        
    # ax_leg.get_xaxis().set_visible(False)
    # ax_leg.get_yaxis().set_visible(False)
    filename = plot_dir + '/activated_fraction_mean_SS.png'
    plt.tight_layout()
    fig.patch.set_alpha(0.)
    fig.savefig(filename,dpi=1000)
    # fig.savefig(filename,dpi=1000,transparent=True)
    
    return fig,ax,[hln_variable,hln_tau0,hln_uniform]

def N_vs_t__settling(N0,ts,Dp_vs_t,Tv_vs_t,Sv_vs_t,rho_p=1000.,ht=1.):
    dt = ts[1] - ts[0]
    N_vs_t = np.zeros(Dp_vs_t.shape)
    N_vs_t[0] = N0
    wp_vs_t = np.zeros(Dp_vs_t.shape)
    for tt in range(0,len(Dp_vs_t)):
        wp_terminal = get_terminal_velocity(Dp_vs_t[tt],rho_p,Tv_vs_t[tt],Sv_vs_t[tt],x_co2=410e-6,p=101325.)
        wp_vs_t[tt] = wp_terminal
        if tt > 0:
            N_vs_t[tt] = N_vs_t[tt-1]*np.exp(-wp_terminal*dt/ht)
    return N_vs_t,  wp_vs_t
        
def get_terminal_velocity(Dp,rho_p,Tv_inf,S_inf,x_co2=410e-6,p=101325.):
    # eqn. 9.49 from Seinfeld and Pandis    
    rho_g = get_air_density(Tv_inf,S_inf,x_co2=x_co2,p=p)
    
    # break out function -- make it jit
    wp_terminal = fsolve(lambda wp: wp - np.sqrt(c.g*(4*Dp*rho_p)/(3*rho_g*get_drag_coefficient(wp-0.,Dp,Tv_inf,S_inf))),1.)
    return wp_terminal

def get_drag_coefficient(velocity_diff,Dp,Tv,S):
    Re = get_particle_Re(velocity_diff,Dp,Tv,S)
    if Re == 0:
        Cd = 0.
    elif Re<=1.:
        Cd = 24/Re
    else:
        Cd = (1 + 0.15*Re**0.687)*24/Re
    return Cd

def get_particle_Re(velocity_diff,Dp,Tv,S,x_co2=410e-6,p=101325.):
    rho_air = get_air_density(Tv,S,x_co2=x_co2,p=p)
    mu = get_dynamic_viscosity_of_air(Tv)
    nu = mu/rho_air
    
    Re = abs(velocity_diff)*Dp/nu
    return Re

def get_dynamic_viscosity_of_air(Tv):
    # http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html
    
    mu0 = 17.15e-6
    T0 = 273.15
    mu = mu0*(Tv/T0)**0.7
    return mu

def get_air_density(Tv,S,x_co2=410e-6,p=101325.):
    x_h2o = get_h2o_mixing_ratio_from_Sv(S,Tv,p=p)
    rho_air = (M_dry_air*(1-x_h2o-x_co2) + x_h2o*M_h2o + x_co2*M_co2)*p/(c.R*Tv)#*(1+x_h2o+x_co2))
    
    return rho_air

def get_h2o_mixing_ratio_from_Sv(S,Tv,p=101325.):
    p_v = get_vapor_pressure(S,Tv)
    mixing_ratio = p_v/p#(p - p_v)
    return mixing_ratio

def get_vapor_pressure(S,Tv):
    p_sat = get_saturation_vapor_pressure(Tv)
    p_v = S*p_sat
    return p_v # Pa

def get_saturation_vapor_pressure(Tv_K):
    Tv = Tv_K-273.15
    p_sat = 611.21*np.exp((18.678 - Tv/234.5)*(Tv/(257.14+Tv)))
    return p_sat # Pa


def get_Dwet(Ddry, kappa, RH, temp):
    import numpy as np
    from scipy.optimize import brentq
    if RH>0 and kappa>0:
        sigma_w = 0.072; rho_w = 1000; M_w = 18/1e3; R=8.314;
        A = 4*sigma_w*M_w/(R*temp*rho_w)
        zero_this = lambda gf: RH/np.exp(A/(Ddry*gf))-(gf**3-1)/(gf**3-(1-kappa))
        return Ddry*brentq(zero_this,1,10000)
    else:
        return Ddry


