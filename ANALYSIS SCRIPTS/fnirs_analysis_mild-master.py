# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 19:47:20 2024

@author: benri
"""

import numpy as np
import mne
import math
# matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import os
import pandas as pd
from collections import defaultdict

# mne.set_config('MNE_BROWSER_BACKEND', 'qt')
#from nirx_movement import mark_aux_movement_bad
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm, statsmodels_to_results
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)
#import pandas as pd
from nilearn.plotting import plot_design_matrix
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.io.snirf import read_snirf_aux_data
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.visualisation import plot_glm_group_topo, plot_glm_surface_projection

# Import StatsModels
import statsmodels.formula.api as smf
#from scipy import signal
from scipy import stats


# Ben preprocessing specific imports
from mne.preprocessing.nirs import optical_density, beer_lambert_law, temporal_derivative_distribution_repair, scalp_coupling_index
from itertools import compress
from mne_nirs.signal_enhancement import short_channel_regression
from mne_modified_beer_lambert_law import mne_modified_beer_lambert_law
from mne.io import read_raw_nirx
from mpl_toolkits.mplot3d import Axes3D


from plot_nirs import plot_nirs_evoked_error
# ---------------------------------------------------------------
# -----------------          Data Parameters            ---------
# ---------------------------------------------------------------
wdir = os.path.dirname(__file__)

# Define Subject Files
# Define Subject Files
root = ''
user = 'Home'
if user == 'Laptop':
    data_root = 'C:/Users/benri/Downloads/'

elif user == 'Desktop':
    data_root = '/home/apclab/Downloads/'
    mild_master_root  = '/home/apclab/Documents/GitHub/MILD-Master'
elif user == 'Home':
    data_root = '/home/ben/Downloads/'
    mild_master_root = '/home/ben/Documents/GitHub/MILD-Master'

all_fnirs_data_folders = [data_root + "2025-01-16/2025-01-16_001",
data_root + "2025-01-20/2025-01-20_001",
data_root + "2025-01-20/2025-01-20_002",
data_root + "2025-01-20/2025-01-20_003",
data_root + "2025-01-21/2025-01-21_001",
data_root + "2025-01-21/2025-01-21_002",
data_root + "2025-01-22/2025-01-22_001",
data_root + "2025-01-22/2025-01-22_002",
data_root + "2025-01-23/2025-01-23_001",
data_root + "2025-01-24/2025-01-24_001",
data_root + "2025-01-24/2025-01-24_002",
data_root + "2025-01-24/2025-01-24_003",
data_root + "2025-01-27/2025-01-27_003",
data_root + "2025-01-29/2025-01-29_001",
data_root + "2025-01-29/2025-01-29_001",
data_root + "2025-01-30/2025-01-30_001",
data_root + "2025-01-30/2025-01-30_002",
data_root + "2025-02-04/2025-02-04_001",
data_root + "2025-02-05/2025-02-05_001",
data_root + "2025-02-05/2025-02-05_002",
data_root + "2025-02-06/2025-02-06_001",
data_root + "2025-02-06/2025-02-06_002",
data_root + "2025-02-10/2025-02-10_001",
data_root + "2025-02-12/2025-02-12_001",
data_root + "2025-02-13/2025-02-13_001",
data_root + "2025-02-13/2025-02-13_002",
data_root + "2025-02-18/2025-02-18_001",                                                                                                            
data_root + "2025-02-19/2025-02-19_001",
data_root + "2025-02-19/2025-02-19_002",
data_root + "2025-02-20/2025-02-20_001",
data_root + "2025-02-20/2025-02-20_002",
data_root + "2025-02-24/2025-02-24_001",
data_root + "2025-02-24/2025-02-24_002"]

# All subject IDs
subject_ID = ['mild_master_1',
'mild_master_3',
'mild_master_4',
'mild_master_5',
'mild_master_7',
'mild_master_6',
'mild_master_8',
'mild_master_9',
'mild_master_10',
'mild_master_11',
'mild_master_12',
'mild_master_13',
'mild_master_14',
'mild_master_15',
'mild_master_16',
'mild_master_17',
'mild_master_18',
'mild_master_19',
'mild_master_20',
'mild_master_21',
'mild_master_22',
'mild_master_23',
'mild_master_24',
'mild_master_25',
'mild_master_26',
'mild_master_27',
'mild_master_28',
'mild_master_29',
'mild_master_30',
'mild_master_31',
'mild_master_32',
'mild_master_33',
'mild_master_34']


# The subjects we would like to run right now
curr_subject_ID = ['mild_master_1',
'mild_master_3',
'mild_master_4',
'mild_master_5']

# 'mild_master_6',
# 'mild_master_8',
# 'mild_master_9',
# 'mild_master_10',
# 'mild_master_11',
# 'mild_master_12',
# 'mild_master_14',
# 'mild_master_15',
# 'mild_master_16',
# 'mild_master_17',
# 'mild_master_18',
# 'mild_master_19',
# 'mild_master_20',
# 'mild_master_22',
# 'mild_master_23',
# 'mild_master_24',
# 'mild_master_25',
# 'mild_master_26',
# 'mild_master_27',
# 'mild_master_28',
# 'mild_master_29',
# 'mild_master_30',
# 'mild_master_31',
# 'mild_master_32',
# 'mild_master_33',
# 'mild_master_34']

curr_folder_indices = [index for index, element in enumerate(subject_ID) if np.isin(element,curr_subject_ID)]
curr_fnirs_data_folders = [all_fnirs_data_folders[i] for i in curr_folder_indices]

glm_dur = 5

n_subjects = len(curr_subject_ID)

n_long_channels = 101
fs = 6.8
tmin, tmax = -2, 16
n_timepoints = math.ceil((tmax - tmin)*fs)
task_type = 'Ben_SvN'


all_subjects_bad_channels = []
group_df = pd.DataFrame()  # To store channel level results

left_hem_channels = [[1,8],[1,7],[2,8],[2,7],[2,6],[3,7],[3,6],[3,5],[7,8],[7,7],[7,18],[7,17],[8,8],[8,7],[8,6],[8,18],
[8,17],[8,16],[9,7],[9,6],[9,5],[9,17],[9,16],[9,15],[10,6],[10,5],[10,16],[10,15],[10,14],[10,21],[15,18],[15,17],[15,22],
[16,18],[16,17],[16,16],[16,22],[16,23],[17,17],[17,16],[17,15],[17,21],[17,22],[17,23],[18,23],[18,22],[18,16],[18,15],[18,21],[19,15],[19,14],[19,21]]

left_hem_channel_names = ["S" + str(value[0]) + "_D" + str(value[1]) + " hbo" for idx, value in enumerate(left_hem_channels)]
# left_hem_channel_names = ["S1_D8 hbo",
#                           "S1_D7 hbo",
#                           "S2_D8 hbo",
#                           "S2_D7 hbo",
#                           "S2_D6 hbo",
#                           "S3_D7 hbo",
#                           "S3_D6 hbo",
#                           "S3_D5 hbo",
#                           "S7_D8 hbo",
#                           "S7_D7 hbo",
#                           "S7_D18 hbo",
#                           "S7_D17 hbo",
#                           "S8_D8 hbo",
#                           "S8_D7 hbo",
#                           "S8_D6 hbo",
#                           "S8_D18 hbo",
#                           "S8_D17 hbo",
#                           "S8_D16 hbo",
#                           "S9_D7 hbo",
#                           "S9_D6 hbo",
#                           "S9_D5 hbo",
#                           "S9_D17 hbo",
#                           "S9_D16 hbo",
#                           "S9_D15 hbo",
#                           "S10_D6 hbo",
#                           "S10_D5 hbo",
#                           "S10_D16 hbo",
#                           "S10_D15 hbo",
#                           "S10_D14 hbo",
#                           "S10_D21 hbo",
#                           "S15_D18 hbo",
#                           "S15_D17 hbo",
#                           "S15_D22 hbo",
#                           "S16_D18 hbo",
#                           "S16_D17 hbo",
#                           "S16_D16 hbo",
#                           "S16_D22 hbo",
#                           "S16_D23 hbo",
#                           "S17_D17 hbo",
#                           "S17_D16 hbo",
#                           "S17_D15 hbo",
#                           "S17_D21 hbo",
#                           "S17_D22 hbo",
#                           "S17_D23 hbo",
#                           "S18_D23 hbo",
#                           "S18_D22 hbo",
#                           "S18_D16 hbo",
#                           "S18_D15 hbo",
#                           "S18_D21 hbo",
#                           "S19_D15 hbo",
#                           "S19_D14 hbo",
#                           "S19_D21 hbo"]


right_hem_channels = [[4,4],[4,3],[4,2],[5,3],[5,2],[5,1],[6,2],[6,1],[11,13],[11,4],[11,3],[11,12],[11,11],[12,4],[12,3],
    [12,2],[12,12],[12,11],[12,10],[13,11],[13,10],[13,9],[13,3],[13,2],[13,1],[14,2],[14,1],[14,10],[14,9],[20,13],[20,12],
    [20,20],[21,20],[21,12],[21,11],[21,19],[22,20],[22,12],[22,11],[22,10],[22,19],[23,11],[23,10],[23,9],[23,19],[24,10],
    [24,9],[24,19]]
right_hem_channel_names = ["S" + str(value[0]) + "_D" + str(value[1]) + " hbo" for idx, value in enumerate(right_hem_channels)]


# set up the arrays to hold all subject data
subject_data_itd5 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd15 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild15 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd5_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd15_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild15_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd5_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd15_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild15_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd5_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd15_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild15_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)


subject_data_itd5_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd15_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild15_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd5_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd15_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild15_GLM = np.full((n_subjects, n_long_channels), np.nan)



# put into a larger array with all subjects data!
subject_data_itd5_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd15_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild15_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)

num_channels_removed = np.full(n_subjects, np.nan)
age = np.full(n_subjects, np.nan)
sex = np.full(n_subjects, np.nan)

range_BH_response = np.zeros((n_subjects, n_long_channels))

all_evokeds = defaultdict(list)
subject_info = []
# loop through all subjects and all sessions (takes a while)

all_epochs = []
for ii, subject_num in enumerate(range(n_subjects)):
    
    
    
    os.environ["OMP_NUM_THREADS"] = "1"
    
    subject = curr_subject_ID[ii]
    task_type = 'mild-master'
    save_dir = "C:/Users/benri/Documents/GitHub/MILD-Master/ANALYSIS SCRIPTS/RESULTS DATA/fNIRS_Data/"
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    plot_dir = "C:/Users/benri/Documents/GitHub/MILD-Master/ANALYSIS SCRIPTS/RESULTS DATA/fNIRS_Plots/"
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)

    plot_steps = False

    # ---------------------------------------------------------------
    # -----------------      Load the Data        ---------
    # ---------------------------------------------------------------
    data = mne.io.read_raw_nirx(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}_config.hdr",
                                    verbose=False, preload=True)
    
    data_snirf = mne.io.read_raw_snirf(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}.snirf",
                                    optode_frame="mri", preload=True)
    
    aux_snirf = read_snirf_aux_data(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}.snirf",
                                    data_snirf)

    # ---------------------------------------------------------------
    # -----------------      Preprocess the Data            ---------
    # ---------------------------------------------------------------  
    
    if subject != "mild_master_1":
        data.annotations.rename({'5.0':	'az_itd=0_az=5',
    '6.0':	'az_itd=0_az=15',
    '7.0':	'az_itd=5_az=0',
    '8.0':	'az_itd=0_az=5',
    '9.0':	'az_itd=0_az=15',
    '10.0':	'az_itd=15_az=0',
    '11.0':	'az_itd=5_az=0',
    '12.0':	'az_itd=15_az=0'})
        data_snirf.annotations.rename({'5':	'az_itd=0_az=5',
    '6':	'az_itd=0_az=15',
    '7':	'az_itd=5_az=0',
    '8':	'az_itd=0_az=5',
    '9':	'az_itd=0_az=15',
    '10':	'az_itd=15_az=0',
    '11':	'az_itd=5_az=0',
    '12':	'az_itd=15_az=0'})
    else:
        data.annotations.rename({'4.0':	'az_itd=0_az=5',
    '5.0':	'az_itd=0_az=15',
    '6.0':	'az_itd=5_az=0',
    '7.0':	'az_itd=0_az=5',
    '8.0':	'az_itd=0_az=15',
    '9.0':	'az_itd=15_az=0',
    '10.0':	'az_itd=5_az=0',
    '11.0':	'az_itd=15_az=0'})
        data_snirf.annotations.rename({'4':	'az_itd=0_az=5',
    '5':	'az_itd=0_az=15',
    '6':	'az_itd=5_az=0',
    '7':	'az_itd=0_az=5',
    '8':	'az_itd=0_az=15',
    '9':	'az_itd=15_az=0',
    '10':	'az_itd=5_az=0',
    '11':	'az_itd=15_az=0'})
        
    
    # Trying out shifting the trigger to the stim onset (rather than cue)
    #data_snirf.annotations.onset = data_snirf.annotations.onset + 0.702
    
    # ---------------------------------------------------------------
    # -------------               Preprocessing             ---------
    # ---------------------------------------------------------------
    
    events, event_dict = mne.events_from_annotations(data, verbose=False)

    if subject != "mild_master_5":
        this_sub_short_regression = True
    else:
        this_sub_short_regression = False

    raw_haemo_temp, null = preprocess_NIRX(data, data_snirf, event_dict,
                                           save=False,
                                           savename=save_dir + f'{subject}_{task_type}_preproc_nirs.fif',
                                           plot_steps=False,
                                           crop=False, crop_low=0, crop_high=0,
                                           events_modification=False, reject=True,
                                           short_regression=this_sub_short_regression, events_from_snirf=False,
                                           drop_short=False, negative_enhancement=False,
                                           snr_thres=1.5, sci_thres=0.6, filter_type='iir', filter_limits=[0.01,0.3])


    if subject != "mild_master_5":
        raw_haemo_short = get_short_channels(raw_haemo_temp)

    raw_haemo_filt = get_long_channels(raw_haemo_temp)

    # extra_regressors = aux_snirf.reset_index(drop=True)
    #
    # order = 4  # You can adjust the order as needed
    # [b, a] = signal.iirfilter(N=order, Wn=0.01 / (0.5 * raw_haemo_filt.info['sfreq']), btype='high', ftype='butter')
    # filtered_signals = extra_regressors.iloc[:, 1:].apply(lambda col: signal.filtfilt(b, a, col), axis=0)


    # ---------------------------------------------------------------
    # -------------               Epoching                  ---------
    # ---------------------------------------------------------------
    reject_criteria = None #dict(hbo=20e-6)#5e-6
    #flat_criteria = dict(hbo=0.05e-6)
    

    epochs = mne.Epochs(raw_haemo_filt, events,  # events_block,
                        event_id=event_dict,  # event_dict_total,
                        tmin=tmin, tmax=tmax,
                        baseline= (tmin, 0),
                        reject = reject_criteria,
                       # flat = flat_criteria,
                        preload=True, detrend=None, verbose=True,
                        on_missing='warn')
    #epochs.plot_drop_log()
    #plt.show()
    #epochs.drop_bad()
    
    all_epochs.append(epochs)
    
    
    n_conditions = 4
    conditions = ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']


    n_conditions = len(conditions)
    evoked_hbo = np.zeros((n_conditions,), dtype=object)
    evoked_hbo_error = np.zeros((n_conditions,), dtype=object)
    evoked_hbr = np.zeros((n_conditions,), dtype=object)
    evoked_hbr_error = np.zeros((n_conditions,), dtype=object)
    vlim = 0.2
    n_conditions = 1
    fig = plt.figure(figsize=(8, 5), dpi=200)
    
    epochs_colors = ['r','g','b','y']

    for i, cond in enumerate(conditions):
        evoked_hbo[i] = epochs[cond].copy().average(picks='hbo')
        evoked_hbo_error[i] = epochs[cond].copy().standard_error(picks='hbo')
        evoked_hbr[i] = epochs[cond].copy().average(picks='hbr')
        evoked_hbr_error[i] = epochs[cond].copy().standard_error(picks='hbr')
        # fig = plot_nirs_evoked_error(fig, evoked_hbo[i], evoked_hbo_error[i],
        #                               colour=epochs_colors[i], ylim=[-vlim, vlim])
        # fig = plot_nirs_evoked_error(fig, evoked_hbr[i], evoked_hbr_error[i],
        #                      colour='b', ylim=[-vlim, vlim], add_legend=True)
        
    # save to all subject array
    for i, cond in enumerate(conditions):
        all_evokeds[cond].append(epochs[cond].average())

    # mark where the bad channels are
    chan_hbo = epochs.copy().pick('hbo').info['ch_names']
    chan_hbr = epochs.copy().pick('hbr').info['ch_names']
    chan_hbo_bad = list(epochs.copy().pick('hbo').info['bads'])
    chan_hbr_bad = list(epochs.copy().pick('hbr').info['bads'])

    chan_indices_bad_hbo = [i for i in range(len(chan_hbo)) if chan_hbo[i] in chan_hbo_bad]
    chan_indices_good_hbo = [i for i in range(len(chan_hbo)) if chan_hbo[i] not in chan_hbo_bad]

    chan_indices_bad_hbr = [i for i in range(len(chan_hbr)) if chan_hbr[i] in chan_hbr_bad]
    chan_indices_good_hbr = [i for i in range(len(chan_hbr)) if chan_hbr[i] not in chan_hbr_bad]

    all_chan = epochs.copy().info['ch_names']
    chan_bad = list(epochs.copy().info['bads'])
    chan_indices_good_all = [i for i in range(len(all_chan)) if all_chan[i] not in chan_bad]
    # split the data into speech, noise, control,
    data_itd5 = epochs["az_itd=5_az=0"].get_data(picks='hbo')
    data_itd15 = epochs["az_itd=15_az=0"].get_data(picks='hbo')
    data_ild5 = epochs["az_itd=0_az=5"].get_data(picks='hbo')
    data_ild15 = epochs["az_itd=0_az=15"].get_data(picks='hbo')
    
    data_itd5_hbr = epochs["az_itd=5_az=0"].get_data(picks='hbr')
    data_itd15_hbr = epochs["az_itd=15_az=0"].get_data(picks='hbr')
    data_ild5_hbr = epochs["az_itd=0_az=5"].get_data(picks='hbr')
    data_ild15_hbr = epochs["az_itd=0_az=15"].get_data(picks='hbr')
    
    # ---------------------------------------------------------------
    # -----------------    Baselining and Averaging         ---------
    # ---------------------------------------------------------------
    data_itd5_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_itd15_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_ild5_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_ild15_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)

    data_itd5_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_itd15_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_ild5_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_ild15_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)

    for ichannel in range(len(chan_indices_good_hbo)):
        
        data_itd5_avg[ichannel,:] = np.nanmean(data_itd5[:,ichannel,:], axis=0)
        data_itd15_avg[ichannel,:] = np.nanmean(data_itd15[:,ichannel,:], axis=0)
        data_ild5_avg[ichannel,:] = np.nanmean(data_ild5[:,ichannel,:], axis=0)
        data_ild15_avg[ichannel,:] = np.nanmean(data_ild15[:,ichannel,:], axis=0)
        
    for ichannel in range(len(chan_indices_good_hbr)):
        data_itd5_avg_hbr[ichannel,:] = np.nanmean(data_itd5_hbr[:,ichannel,:], axis=0)
        data_itd15_avg_hbr[ichannel,:] = np.nanmean(data_itd15_hbr[:,ichannel,:], axis=0)
        data_ild5_avg_hbr[ichannel,:] = np.nanmean(data_ild5_hbr[:,ichannel,:], axis=0)
        data_ild15_avg_hbr[ichannel,:] = np.nanmean(data_ild15_hbr[:,ichannel,:], axis=0)
    
    # need to mark the indices where the good channels are!

    # put into a larger array with all subjects data!
    subject_data_itd5[ii, chan_indices_good_hbo, :] = 1e6*data_itd5_avg.copy()
    subject_data_itd15[ii, chan_indices_good_hbo, :] = 1e6*data_itd15_avg.copy()
    subject_data_ild5[ii, chan_indices_good_hbo, :] = 1e6*data_ild5_avg.copy()
    subject_data_ild15[ii, chan_indices_good_hbo, :] = 1e6*data_ild15_avg.copy()
  
    subject_data_itd5_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_itd5_avg_hbr.copy()
    subject_data_itd15_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_itd15_avg_hbr.copy()
    subject_data_ild5_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_ild5_avg_hbr.copy()
    subject_data_ild15_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_ild15_avg_hbr.copy()
    
    subject_data_itd5_baselined = subject_data_itd5.copy()
    subject_data_itd15_baselined = subject_data_itd15.copy()
    subject_data_ild5_baselined = subject_data_ild5.copy()
    subject_data_ild15_baselined = subject_data_ild15.copy()
    
    subject_data_itd5_hbr_baselined = subject_data_itd5_hbr.copy()
    subject_data_itd15_hbr_baselined = subject_data_itd15_hbr.copy()
    subject_data_ild5_hbr_baselined = subject_data_ild5_hbr.copy()
    subject_data_ild15_hbr_baselined = subject_data_ild15_hbr.copy()

    # ---------------------------------------------------------------
    # -----------------     GLM                             ---------
    # ---------------------------------------------------------------

    # try to remove some of the conditions in the raw_haemo_filt annotations

    # raw_haemo_temp_crop = raw_haemo.crop(tmin=raw_haemo.annotations.onset[45] - 15)

    # raw_haemo_short_crop = get_short_channels(raw_haemo_temp_crop)
    # raw_haemo_filt_crop = get_long_channels(raw_haemo_temp_crop)

    # # try cropping
    # raw_haemo_filt_crop.resample(5)
    # raw_haemo_short_crop.resample(5)
    
    raw_haemo_filt_for_glm = get_long_channels(raw_haemo_filt).copy()
    
    # drop bad channels from glm
    all_subjects_bad_channels.append([epochs.copy().info['bads']])
    
    raw_haemo_filt_for_glm.drop_channels(epochs.copy().info['bads'])

    raw_haemo_filt_for_glm.annotations.set_durations(glm_dur)
    
    

    design_matrix = make_first_level_design_matrix(raw_haemo_filt_for_glm,
                                                        drift_model=None,
                                                        high_pass=0.01,  # Must be specified per experiment
                                                        hrf_model='spm',
                                                        stim_dur=glm_dur)
    # # add_regs=filtered_signals)

    design_matrix["Linear"] = np.arange(0, np.shape(design_matrix)[0])
    #design_matrix_hbo["ShortHbO"] = np.mean(raw_haemo_short.copy().pick(picks="hbo").get_data(), axis=0)

    #design_matrix_hbr["ShortHbR"] = np.mean(raw_haemo_short.copy().pick(picks="hbr").get_data(), axis=0)
    min_max_scaler = MinMaxScaler()
    X_minmax = min_max_scaler.fit_transform(design_matrix)
    design_matrix_min_max = pd.DataFrame(X_minmax, columns=design_matrix.columns.tolist())

    # pre-whiten
    raw_haemo_filt_for_glm._data = np.subtract(raw_haemo_filt_for_glm._data,
                                            np.nanmean(raw_haemo_filt_for_glm._data, axis=1)[:, np.newaxis])
    glm_est = run_glm(raw_haemo_filt_for_glm, design_matrix, noise_model='ar1')

    # record the glm est for each condition, for each subject
    # will adjust the beta values by the BH correction method

    glm_est_df = glm_est.pick(picks='data').to_dataframe()

    # Append stats df to a group results
    individual_results = glm_est.to_dataframe()
    individual_results["ID"] = subject
    # Convert to uM for nicer plotting below.
    individual_results["theta"] = [t * 1.e6 for t in individual_results["theta"]]

    group_df = pd.concat([group_df, individual_results], ignore_index=True)

    # GLM Topoplot just this participant

    glm_hbo = glm_est.copy().pick(picks="hbo")
    conditions_to_plot = ['az_itd=0_az=5','az_itd=0_az=15']

    this_sub_left_hem = [idx for idx, value in enumerate(glm_hbo.ch_names) if value in left_hem_channel_names]
    this_sub_right_hem = [idx for idx, value in enumerate(glm_hbo.ch_names) if value in right_hem_channel_names]

    fig_topo, ax_topo = plt.subplots(nrows=1, ncols=len(conditions_to_plot))
    for icond, condition in enumerate(conditions_to_plot):
        glm_hbo.copy().pick(this_sub_left_hem).plot_topo(conditions=condition, axes=ax_topo[icond], colorbar=False, vlim=(-0.2, 0.2),sensors=False)
        glm_hbo.copy().pick(this_sub_right_hem).plot_topo(conditions=condition, axes=ax_topo[icond], colorbar=False,vlim=(-0.2, 0.2),sensors=False)
    plt.savefig(mild_master_root + "/CASUAL FIGURES/" + subject + "_individual_topo.png")

##############################
## Take mean during stim ####
############################

# Take Means
index_stim_start = int(-tmin*fs)
index_stim_end = int((-tmin + 11.6)*fs)
# itd5
mean_during_stim_itd5 = np.nanmean(subject_data_itd5_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd5_hbr = np.nanmean(subject_data_itd5_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# itd15
mean_during_stim_itd15 = np.nanmean(subject_data_itd15_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd15_hbr = np.nanmean(subject_data_itd15_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# ild5
mean_during_stim_ild5 = np.nanmean(subject_data_ild5_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild5_hbr = np.nanmean(subject_data_ild5_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

#ild15
mean_during_stim_ild15 = np.nanmean(subject_data_ild15_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild15_hbr = np.nanmean(subject_data_ild15_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

## Save breath uncorrected and corrected GLM data

# block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd5_baselined.shape], names = names)
itd5_df = pd.DataFrame({'subject_data_itd5':subject_data_itd5_baselined.flatten()},index=index)['subject_data_itd5']
itd15_df = pd.DataFrame({'subject_data_itd15':subject_data_itd15_baselined.flatten()},index=index)['subject_data_itd15']
ild5_df = pd.DataFrame({'subject_data_ild5':subject_data_ild5_baselined.flatten()},index=index)['subject_data_ild5']
ild15_df = pd.DataFrame({'subject_data_ild15':subject_data_ild15_baselined.flatten()},index=index)['subject_data_ild15']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_uncorr_block_average_hbo.csv',index=True)


# block average means HBO
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd5.shape], names = names)
itd5_df = pd.DataFrame({'mean_during_stim_itd5':mean_during_stim_itd5.flatten()},index=index)['mean_during_stim_itd5']
itd15_df = pd.DataFrame({'mean_during_stim_itd15':mean_during_stim_itd15.flatten()},index=index)['mean_during_stim_itd15']
ild5_df = pd.DataFrame({'mean_during_stim_ild5':mean_during_stim_ild5.flatten()},index=index)['mean_during_stim_ild5']
ild15_df = pd.DataFrame({'mean_during_stim_ild15':mean_during_stim_ild15.flatten()},index=index)['mean_during_stim_ild15']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_hbo.csv',index=True)


# Uncorrected block average means HBR
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd5_hbr.shape], names = names)
itd5_df = pd.DataFrame({'mean_during_stim_itd5_hbr':mean_during_stim_itd5_hbr.flatten()},index=index)['mean_during_stim_itd5_hbr']
itd15_df = pd.DataFrame({'mean_during_stim_itd15_hbr':mean_during_stim_itd15_hbr.flatten()},index=index)['mean_during_stim_itd15_hbr']
ild5_df = pd.DataFrame({'mean_during_stim_ild5_hbr':mean_during_stim_ild5_hbr.flatten()},index=index)['mean_during_stim_ild5_hbr']
ild15_df = pd.DataFrame({'mean_during_stim_ild15_hbr':mean_during_stim_ild15_hbr.flatten()},index=index)['mean_during_stim_ild15_hbr']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_hbr.csv',index=True)

# ---------------------------------------------------------------
# -----------------     Subject Averaging                ---------
# ---------------------------------------------------------------
# # for each subject, take the beta values and compute a mean and standard error
subject_data_itd5_GLM_mean = np.nanmean(1e6*subject_data_itd5_GLM, axis=0)
subject_data_itd15_GLM_mean = np.nanmean(1e6*subject_data_itd15_GLM, axis=0)
subject_data_ild5_GLM_mean = np.nanmean(1e6*subject_data_ild5_GLM, axis=0)
subject_data_ild15_GLM_mean = np.nanmean(1e6*subject_data_ild15_GLM, axis=0)

subject_data_itd5_GLM_std = np.nanstd(1e6*subject_data_itd5_GLM, axis=0) / (np.sqrt(n_subjects)-1)
subject_data_itd15_GLM_std = np.nanstd(1e6*subject_data_itd15_GLM, axis=0) / (np.sqrt(n_subjects)-1)
subject_data_ild5_GLM_std = np.nanstd(1e6*subject_data_ild5_GLM, axis=0) / (np.sqrt(n_subjects)-1)
subject_data_ild15_GLM_std = np.nanstd(1e6*subject_data_ild15_GLM, axis=0) / (np.sqrt(n_subjects)-1)


# ---------------------------------------------------------------
# -----------------     PLotting GLM Averages           ---------
# ---------------------------------------------------------------

caxis_lim = 0.2

fig, ax_topo = plt.subplots(1, 4)
groups_single_chroma = dict(
    Left_Hemisphere=picks_pair_to_idx(raw_haemo_filt.copy().pick(picks='hbo'), left_hem_channels,
                                      on_missing='warning'),
    Right_Hemisphere=picks_pair_to_idx(raw_haemo_filt.copy().pick(picks='hbo'), right_hem_channels,
                                       on_missing='warning'))
# Run group level model and convert to dataframe
group_results = group_df.query("Condition in ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']")

import seaborn as sns
sns.catplot(x="Condition",y="theta",col="ID",hue="Chroma",data=group_results,col_wrap=5,errorbar=None,palette="muted",height=4, s=10)
plt.savefig(mild_master_root + "/CASUAL FIGURES/beta_values_by_participant.png")


group_theta_for_catplot = group_results.groupby(by=['Chroma','ch_name','ID','Condition'],as_index=False)['theta'].mean()
sns.catplot(x="Condition",y="theta",hue = "ID", data=group_theta_for_catplot.query("Chroma == 'hbo'"), errorbar=None, height=7, s=10, legend=False)
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_catplot.png")

group_theta_for_topoplot = group_results.query("Chroma in ['hbo']").groupby(by=['ch_name','Condition'],as_index=False)['theta'].mean()
fig, topo_axes = plt.subplots(nrows=1, ncols=1,figsize=(18,36))
this_info = raw_haemo_filt.copy().pick(picks="hbo")
this_info.drop_channels([i for i in this_info.ch_names if i not in np.unique(group_theta_for_topoplot['ch_name'])])
this_info = this_info.info
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']")['theta'],this_info,sensors=False, axes = topo_axes)
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot.png")

# NOW SPLIT IT UP BY CHANNEL

# ch_model = smf.mixedlm("theta ~ -1 + ch_name:Chroma:Condition",group_results, groups=group_results["ID"]).fit(method='nm')
# ch_model_df = statsmodels_to_results(ch_model)

# sns.catplot(x="Condition",y="Coef.",data=ch_model_df.query("Chroma == 'hbo'"), errorbar=None, palette="muted", height=4, s=10)
# plt.savefig(mild_master_root + "/CASUAL FIGURES/beta_values_by_condition.png")
### Plot group results
# plot_glm_group_topo(raw_haemo_filt.copy().pick(picks="hbo"),
#                     ch_model_df.query("Condition in ['az_itd=5_az=0']"),                ####Change here####
#                     colorbar=False, axes=ax_topo[0],
#                     vlim=(-caxis_lim, caxis_lim), cmap='summer')

# plot_glm_group_topo(raw_haemo_filt.copy().pick(picks="hbo").pick(groups_single_chroma['Right_Hemisphere']),
#                     ch_model_df.query("Condition in ['az_itd=5_az=0']"),              ####Change here#####
#                     colorbar=True, axes=ax_topo[0],
#                     vlim=(-caxis_lim, caxis_lim), cmap='summer')
# plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topo.png")

caxis_lim = 0.1

fig, (ax1,ax2,ax3,ax4) = plt.subplots(1, 4)
im, _ = mne.viz.plot_topomap(subject_data_itd5_GLM_mean, all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='cubic',
                              vlim=(-caxis_lim, caxis_lim), cmap ='summer', axes=ax1, show=True)
#cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label('Beta (a.u.)')
ax1.set_title('GLM Beta: itd5')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_itd15_GLM_mean, all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='cubic',
                              vlim=(-caxis_lim, caxis_lim), cmap ='summer', axes=ax2, show=False)
#cbar = fig.colorbar(im, ax=ax2)
#cbar.set_label('Beta (a.u.)')
ax2.set_title('GLM Beta: itd15')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_ild5_GLM_mean, all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='cubic',
                              vlim=(-caxis_lim, caxis_lim), cmap ='summer', axes=ax3, show=False)
#cbar = fig.colorbar(im, ax=ax3)
#cbar.set_label('Beta (a.u.)')
ax3.set_title('GLM Beta: ild5')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_ild15_GLM_mean, all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='cubic',
                              vlim=(-caxis_lim, caxis_lim), cmap ='summer', axes=ax4, show=False)
cbar = fig.colorbar(im, ax=ax4)
cbar.set_label('Beta (a.u.)')
ax4.set_title('GLM Beta: ild15')





# ---------------------------------------------------------------
# -----------------     Topomap of Mean HbO             ---------
caxis_min = -0.1
caxis_max = 0.1

all_epochs[0].pick('hbo').info['bads'] = []

fig, (ax1,ax2,ax3,ax4) = plt.subplots(1, 4)
im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_itd5, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='nearest',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax1, show=False)
#cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label('Beta (a.u.)')
ax1.set_title('Mean Delta HbO: itd5')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_itd15, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='nearest',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax2, show=False)
#cbar = fig.colorbar(im, ax=ax2)
#cbar.set_label('Beta (a.u.)')
ax2.set_title('Mean Delta HbO: itd15')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_ild5, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='nearest',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax3, show=False)
#cbar = fig.colorbar(im, ax=ax3)
#cbar.set_label('Beta (a.u.)')
ax3.set_title('Mean Delta HbO: ild5')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_ild15, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='nearest',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax4, show=False)
cbar = fig.colorbar(im, ax=ax4)
cbar.set_label('Mean DeltaHbO')
ax4.set_title('Mean Delta HbO: ild15')
plt.show()

caxis_min = -0.05
caxis_max = 0.05

mean_during_stim_itd5[np.isnan(mean_during_stim_itd5)] = 0
mean_during_stim_itd15[np.isnan(mean_during_stim_itd15)] = 0
mean_during_stim_ild5[np.isnan(mean_during_stim_ild5)] = 0
mean_during_stim_ild15[np.isnan(mean_during_stim_ild15)] = 0

all_epochs[0].pick('hbo').info['bads'] = []

fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1, 5)
im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_itd5, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='cubic',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax1, show=False)
#cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label('Beta (a.u.)')
ax1.set_title('Mean Delta HbO: itd5')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_itd15, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local',  image_interp='cubic',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax2, show=False)
#cbar = fig.colorbar(im, ax=ax2)
#cbar.set_label('Beta (a.u.)')
ax2.set_title('Mean Delta HbO: itd15')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_ild5, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='cubic',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax3, show=False)
#cbar = fig.colorbar(im, ax=ax3)
#cbar.set_label('Beta (a.u.)')
ax3.set_title('Mean Delta HbO: ild5')
plt.show()

im, _ = mne.viz.plot_topomap(np.nanmean(mean_during_stim_ild15, axis=0), all_epochs[0].pick('hbo').info,
                      extrapolate='local', image_interp='cubic',
                              vlim=(caxis_min, caxis_max), cmap ='summer', axes=ax4, show=False)
cbar = fig.colorbar(im, ax=ax4, cax=ax5)
cbar.set_label('Mean DeltaHbO')
ax4.set_title('Mean Delta HbO: ild15')
plt.show()

# ---------------------------------------------------------------
# -----------------     Scatterplot of Mean HbO         ---------
# ---------------------------------------------------------------
# Build Data Frame
mean_hbo_all_conditions = np.stack((mean_during_stim_itd5,mean_during_stim_itd15,
                                          mean_during_stim_ild5,mean_during_stim_ild15),axis=0)

mean_hbo_df = pd.DataFrame(columns=['MeanHbO','S','Condition'])

mean_hbo_df['MeanHbO'] = pd.Series(np.ravel(mean_hbo_all_conditions))

mean_hbo_df['S'] = pd.Series(np.tile(np.array([np.repeat(isub,np.size(mean_hbo_all_conditions,axis=2)) for isub in range(len(curr_subject_ID))]).ravel(),np.size(mean_hbo_all_conditions,axis=0)))
mean_hbo_df['Condition'] = pd.Series(np.array([np.repeat(cond,np.size(mean_hbo_all_conditions,axis=1)*np.size(mean_hbo_all_conditions,axis=2)) for cond in conditions]).ravel())

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 5))
import seaborn as sns
sns.lineplot(x="Condition", y="MeanHbO", hue="S", data=mean_hbo_df)

# ---------------------------------------------------------------
# -----------------     Topomap of Group Block Averages ---------
# ---------------------------------------------------------------
conditions = ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']
n_conditions = len(conditions)
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(17, 5))
lims = dict(hbo=[-0.2, 0.2], hbr=[-0.2, 0.2])

for pick, color in zip(["hbo", "hbr"], ["r", "b"]):
    for idx, cond in enumerate(conditions):
        mne.viz.plot_compare_evokeds(
            {cond: all_evokeds[cond]},
            combine="mean",
            picks=pick,
            axes=axes[idx],
            show=False,
            colors=[color],
            legend=False,
            ylim=lims,
            ci=0.95,
        )
        axes[idx].set_title(f"{cond}")
axes[0].legend(["Oxyhaemoglobin", "Deoxyhaemoglobin"])

