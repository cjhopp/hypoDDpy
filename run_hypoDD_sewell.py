#!/usr/bin/env python

""" This program runs hypoDDpy."""

import glob
import sys, os
# sys.path.insert('/home/chet/hypoDDpy')

from hypoddpy import HypoDDRelocator

# Location of input QML catalog
cat_file = '/Volumes/GeoPhysics_07/users-data/hoppche/detections/det_cat_mcc0.4_shift0.2_ALL_LOCATED_uncert0.05.xml'
# time slice directory
wav_dir = '/Volumes/GeoPhysics_07/users-data/hoppche/stefan_sac/all_MF_dets'

# station files dir
sta_file = '/Users/home/hoppche/data/stations/Mercury_Network_staxml.xml'

# working dir
work_dir = '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD'

# output dir and catalog file
out_dir = work_dir
out_file = 'det_cat_mcc0.4_shift0.2_ALL_LOCATED_HypoDD.xml'
cc_plot_dir = out_dir + '/cc_plots/'

# number of cores for parallel cross-correlation processing
ncores = 10

### ph2dt Settings
ph2dt_sets = {
    # min pick weight. Leave as 0, don't know why this would change?
    'MINWGHT' : 0,  
    # maximum distance (in km) between event pair and station. 
    # Set large to include all stations
    'MAXDIST'  : 50,
    # Maximum hypocentral separation between event pairs in kms. Set to ensure 
    # events within same spatial cluster are considered together whilst excluding
    # events that are obviously in different area
    'MAXSEP' : 10, 
    # Maximum number of neighbours per event. Should be high to allow all possible
    # events within geographic cluster defined by other parameters
    'MAXNGH' : 1000,
    # Minimum number of links required to define a neighbour
    'MINLNK' : 8,
    # Minimum number of links per pair 
    'MINOBS' : 1,
    # Max number of links per pair. 
    # Should set to total number of stations to consider all phase pairs 
    # within geographic cluster
    'MAXOBS' : 40
            }

### HypoDD Settings
hypodd_sets = {
    # 1 = cross-corr only, 2 = absolute (cat) data only, 3 = cross-corr and catalog
    'IDAT' : 3,
    # 1 = P-phase, 2 = S-phase, 3 = P & S phase
    'IPHASE' : 1,
    # Maximum distance between centroid of event cluster and stations
    'DIST' : 50,
    # Minimum number of x-corr or catalog links per event pair to form a continuous cluster
    # Set to 0 to disable clustering within hypoDD
    'OBSCC' : 0,
    'OBSCT' : 0,
    #  min and max distance between event pairs and stations
    # Set to -999 to disable
    'MINDS' : -999,
    'MAXDS' : -999,
    # Maximum azimuthal gap between individual events pairs and stations
    # Set to -999 to disable
    'MAXGAP' : -999,
    # Initial locations. 1 = start from cluster centroid, 2 = start from catalog locations
    'ISTART' : 2,
    # Iterations
    # List each iteration as a string in following order... 
    # NITER = number of iterations for this set of parameters
    # WTCCP, WTCCS = Weight for cross-corr (P then S) -999 to not use this data
    # WTCTP, WTCTS = Weight for catalog (P then S) -999 to not use this data
    # WRCC, WRCT = Cutoffs to remove outliers (cross-corr then cat). 
        # 0-1 = static (number of secs) >=1 = dynamic, multiple of standard dev
    # WDCC, WDCT = Max event separation distance (cross-corr then cat) -999 to deactivate
    # DAMP = damping. Aim for condition numbers between about 40-80
                    #   Cross-corr Data   #    Catalog Data    #  
             # NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP 
    'iters' : ["   5  0.01  0.01 -999 -999   1.0 0.005 -999 -999   20",
               "   5  0.01  0.01 -999 -999   1.0 0.005 -999    4   20",
               "   5  1.00  0.01 -999    2  0.01 0.005 -999    4   20",
               "   5  1.00  0.01 -999    2  0.01 0.005 -999    4   20",
               "   5  1.00  0.01 -999  0.5  0.01 0.005 -999    4   20"]
            }

### Cross-correlation Plotting
# Event pair interval to output plots of cross-correlations (use 0 to turn off)
cc_plot_int = 0


""" 
    Relocator parameters
    :param working_dir: The working directory where all temporary and final
        files will be placed.
    :param cc_time_before: Time to start cross correlation before pick time
        in seconds.
    :param cc_time_after: Time to start cross correlation after pick time
        in seconds.
    :param cc_maxlag: Maximum lag time tested during cross correlation.
    :param cc_filter_min_freq: Lower corner frequency for the Butterworth
        bandpass filter to be applied during cross correlation.
    :param cc_filter_max_freq: Upper corner frequency for the Butterworth
        bandpass filter to be applied during cross correlation.
    :param cc_p_phase_weighting: The cross correlation travel time
        differences can be calculated on several channels. This dict
        specified which channels to calculate it for and how to weight the
        channels to determine the final cross correlated traveltime between
        two events. This assumes the waveform data adheres to the SEED
        naming convention.
        This dict applies to all P phase picks.
        Examples:
            {"Z": 1.0} - Only use the vertical channel.
            {"E": 1.0, "N": 1.0} - Use east and north channel and weight
                                   them equally.
    :param cc_s_phase_weighting: See cc_p_phase_weighting. Just for S
        phases.
    :param cc_min_allowed_cross_corr_coeff: The minimum allowed
        cross-correlation coefficient for a differential travel time to be
        accepted.
    """

relocator = HypoDDRelocator(
    working_dir=work_dir,
    cc_time_before=0.1,
    cc_time_after=0.2,
    cc_maxlag=0.2,
    cc_filter_min_freq=3.0,
    cc_filter_max_freq=20.0,
    cc_p_phase_weighting={"Z": 1.0},
    cc_s_phase_weighting= {},
    cc_min_allowed_cross_corr_coeff=0.5,
    ph2dt_sets=ph2dt_sets, 
    hypodd_sets=hypodd_sets,
    cc_plot_int=cc_plot_int,
    cc_plot_dir=cc_plot_dir)

# Setup the velocity model. This is just a constant velocity model.
relocator.setup_velocity_model(
    model_type="layered_p_velocity_with_constant_vp_vs_ratio",
    layer_tops=[(-0.6, 1.9),
                (0.2, 2.6),
                (1.0, 3.5),
                (1.5, 3.6),
                (2.0, 3.9),
                (3.0, 4.9),
                (5.0, 5.4),
                (10.0, 5.8),
                (20.0, 6.9),
                (50.0, 7.4)],
    vp_vs_ratio=1.7)

# Add the necessary files. Call a function multiple times if necessary.
relocator.add_event_files(cat_file)

data_file_list = []
for root, dirs, files in os.walk(wav_dir):
    for file in files:
        data_file_list.append(os.path.join(root, file))
relocator.add_waveform_files(data_file_list)
relocator.add_station_files(glob.glob(sta_file))

# Start the relocation with the desired output file.
relocator.start_relocation(output_event_file=os.path.join(out_dir, out_file),
                           ncores=ncores)