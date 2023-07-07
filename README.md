# Particle-based model of turbulent microphysics in the MTU Pi Chamber
[![DOI](https://zenodo.org/badge/542801260.svg)](https://zenodo.org/badge/latestdoi/542801260)

  * Models fluctuations in temperature and supersaturation due to turbulent flow within clouds.
  * Calculates droplet growth and cloud droplet activation of aerosol particles exposed to turbulent flow.

Anderson, J.C., Beeler, P., Ovchinnikov, M., Cantrell, W., Krueger, S., Shaw, R.A., Yang, F. and Fierce, L., 2022. Enhancements in cloud condensation nuclei concentrations from turbulent fluctuations in supersaturation. arXiv preprint arXiv:2210.15766.


## Dependencies

  * numpy
  * scipy
  * netCDF4

## Instructions

1. Download large eddy simulation output, which is provided [here](https://portal.nersc.gov/project/m1657/LagrangianDroplets/).
2. Clone repository
3. Open main.py
4. Ensure that the working directory contains main.py, microphysics.py, parcels.py, process.py, and run.py 
4. Move "OUT_3D" folder to the working directory
5. Enter desired input parameters (lines 13-20)
6. Run

## Function Descriptions

  * run.make_directories
  	 * Creates folders called "parcel_traces", "particle_traces", and figures in current working directory.
  	 * Creates folders for each average supersaturation within "particle_traces".
  	   * example: Particle traces run with 3% average supersaturation will be saved in "/particle_traces/avgSS_000030".
  * parcels.trajectories
    * Reads output from large eddy simulation and calculates particle trajectories.
    * Parcel x, y, z positions, as well as temperature and supersaturation are saved in "/parcel_traces".
      * example: Data for first parcel is saved in "parcel_traces/parcel_traces_dynamic000000.pkl".
   * run.main
     * Calculates particle/droplet properties based on parcel trajectories.
     * Returns a dict with size = (average supersaturations, parcels). Contains the following data and units:
       * traces
         * t - time, seconds
         * Dp - particle diameter, meters
         * Tp - particle temperature, Kelvin
         * mixing_ratio - water vapor mixing ratio, kg water per kg air
         * T - air temperature, Kelvin
         * SS - supersaturation, unitless
       * props
         * Ddry - dry diameter of particles, meters
         * kappa - hygroscopiscity of particles, unitless
         * density - particle density, kg per cubic meter
         * N - total number concentration, number per cubic meter
         * Dp_crit - critical diameter for droplet activation, meters
         * SS_crit - critical supersaturation for droplet activation, unitless
   * process.plot_Nccn_timeseries
     * Plots fraction of activated particles as a function of time for each average supersaturation.
     * Plot is saved in "/figures"
   * process.plot_S_Dp_timeseries
     * Plots timeseries of supersaturation, droplet diameter, and droplet diameter/critical diameter of one parcel and one average supersaturation.
     * By default the first parcel and the first average supersaturation are plotted.
     * Plot is saved in "/figures"
   * process.plot_distribution_SS_Dp
     * Plots histogram of supersaturation and droplet diameter for each average supersaturation.
     * Plot is saved in "/figures"
   * process.plot_spatial_trajectory
     * Plots trajectory of one parcel
     * By default the first air parcel is plotted.
     * Plot is saved in "/figures"
   * process.plot_activated_fraction_mean_SS
     * Plots activated fraction as a function of average supersaturation.
     * Plot is saved in "/figures"

## License
See the [LICENSE](https://github.com/lfierce2/LagrangianDroplets/blob/main/LICENSE.md) file for license rights and limitations (GNU).
