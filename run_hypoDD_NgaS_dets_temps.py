#!/usr/bin/env python

""" This program runs hypoDDpy."""

from glob import glob
from obspy import read_events
import csv
import os
import gc

from hypoddpy.hypodd_relocator import HypoDDRelocator

# Location of input QML catalog
# cat_file = '/Volumes/GeoPhysics_07/users-data/hoppche/detections/det_cat_mcc0.4_shift0.2_ALL_LOCATED_uncert0.05.xml'
cat_files_S = [
    '/Volumes/GeoPhysics_07/users-data/hoppche/catalogs/qml/year_long/hypoDD/det_cat_mcc0.5_shift0.2_NGAS_HypoDD_cc0.7_newlocs.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2012/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2012a_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2012/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2012b_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2012/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2012e_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2012/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2012f_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013a_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013b_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013c_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013d_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013e_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2013f_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014a_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014b_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014c_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014d_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014e_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2014f_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015a_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015b_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015c_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015d_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015e_LOCATED_justlocs_uncert0.05.xml',
    '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/location_cats/NgaS/det_cat_mcc0.4_shift0.2_NgaS_2015f_LOCATED_justlocs_uncert0.05.xml',
    ]
out_dirs_S = [
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2012a',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2012b',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2012e',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2012f',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013a',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013b',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013c',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013d',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013e',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2013f',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014a',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014b',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014c',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014d',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014e',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2014f',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015a',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015b',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015c',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015d',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015e',
    '/Volumes/GeoPhysics_07/users-data/hoppche/hypoDD/all_dets_NgaS_2015f',
    ]
# time slice directory
wav_dir = '/Volumes/GeoPhysics_07/users-data/hoppche/stefan_sac/all_MF_dets'

# station files dir
sta_file = '/Users/home/hoppche/data/stations/Mercury_Network_staxml.xml'

# number of cores for parallel cross-correlation processing
ncores = 20

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
    'MAXSEP' : 2,
    # Maximum number of neighbours per event. Should be high to allow all possible
    # events within geographic cluster defined by other parameters
    'MAXNGH' : 1000,
    # Minimum number of links required to define a neighbour
    'MINLNK' : 6,
    # Minimum number of links per pair
    'MINOBS' : 1,
    # Max number of links per pair.
    # Should set to total number of stations to consider all phase pairs
    # within geographic cluster
    'MAXOBS' : 45
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
    'iters' : ["   5  0.10  -999    2    1   1.0  -999    2    1   150",
               "   5  0.30  -999    2    1   0.7  -999    2    1   100",
               "   5  0.50  -999    2    1  0.50  -999    2    1    70",
               "   5  0.70  -999    2  0.75 0.30  -999    2  0.75   50",
               "   5  1.00  -999    2  0.5  0.01  -999    2  0.5    30"]
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
for f_list, outdirs in zip([cat_files_S],
                           [out_dirs_S]):
    for cat_file, outdir in zip(f_list[1:], outdirs):
        print('Running file {}'.format(cat_file))
        relocator = HypoDDRelocator(
            working_dir=outdir,
            cc_time_before=0.1,
            cc_time_after=0.5,
            cc_maxlag=0.3,
            cc_filter_min_freq=3.0,
            cc_filter_max_freq=20.0,
            cc_p_phase_weighting={"Z": 1.0},
            cc_s_phase_weighting= {},
            cc_min_allowed_cross_corr_coeff=0.7,
            ph2dt_sets=ph2dt_sets,
            hypodd_sets=hypodd_sets)
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
        print('Adding event files to relocator object')
        relocator.add_event_files([cat_file, f_list[0]])
        print('Creating list of all self_detection wav directories.')
        self_files = [
            '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2012/selfs_rotnga_2012.txt',
            '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2013/selfs_rotnga_2013.txt',
            '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2014/selfs_rotnga_2014.txt',
            '/Volumes/GeoPhysics_07/users-data/hoppche/detections/2015/selfs_rotnga_2015.txt']
        selfs = []
        for self_file in self_files:
            with open(self_file, 'r') as f:
                rdr = csv.reader(f)
                for row in rdr:
                    selfs.append(str(row[0]))
        data_file_list = []
        for self in selfs:
            print('Adding template: %s to file list' % os.path.join(wav_dir,
                                                                    self, '*'))
            data_file_list.extend(glob(os.path.join(wav_dir, self, '*')))
        temp_wav_len = len(data_file_list)
        print('Finding only event wavs in catalog')
        print('Reading catalog...')
        cat = read_events(cat_file)
        names = [str(ev.resource_id).split('/')[-1] for ev in cat]
        del cat
        gc.collect()
        for name in names:
            print('Adding detection: %s to file list' % os.path.join(wav_dir,
                                                                     name, '*'))
            data_file_list.extend(glob(os.path.join(wav_dir, name, '*')))
        print('Length of wav list with just temps: {}'.format(temp_wav_len))
        print('Length of waveforms list after dets: {}'.format(
            len(data_file_list)))
        print('Adding wavs to relocator object')
        relocator.add_waveform_files(data_file_list)
        print('Adding station files to relocator object')
        relocator.add_station_files(glob(sta_file))
        print('Starting relocation run')
        # Start the relocation with the desired output file.
        relocator.start_relocation(output_event_file=os.path.join(
            outdir, '{}_hypoDD_0.7.xml'.format(
                '_'.join(outdir.split('_')[-2:]))),
                                   ncores=ncores)