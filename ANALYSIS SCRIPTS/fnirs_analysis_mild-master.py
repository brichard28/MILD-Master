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

mne.set_config('MNE_BROWSER_BACKEND', 'qt')
#from nirx_movement import mark_aux_movement_bad
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)
#import pandas as pd
from nilearn.plotting import plot_design_matrix
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.io.snirf import read_snirf_aux_data
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.visualisation import plot_glm_group_topo

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
data_root = 'C:/Users/benri/Downloads/' 

# Define Subject Files
root = ''
user = 'Noptop'
if user == 'Laptop':
    data_root = 'C:/Users/benri/Downloads/'

else:
    data_root = 'D:\\Downloads\\'

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
'mild_master_5',
'mild_master_6',
'mild_master_8',
'mild_master_9',
'mild_master_10',
'mild_master_11',
'mild_master_12',
'mild_master_14',
'mild_master_15',
'mild_master_16',
'mild_master_17',
'mild_master_18',
'mild_master_19',
'mild_master_20',
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

curr_folder_indices = [index for index, element in enumerate(subject_ID) if np.isin(element,curr_subject_ID)]
curr_fnirs_data_folders = [all_fnirs_data_folders[i] for i in curr_folder_indices]

masker_type = 'speech' # type of masker to analyze on this run
glm_dur = 5

n_subjects = len(curr_subject_ID)

n_long_channels = 101
fs = 6.8
tmin, tmax = -2, 16
n_timepoints = math.ceil((tmax - tmin)*fs)
task_type = 'Ben_SvN'


all_subjects_bad_channels = []

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
preprocessing_type = "Eli"

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
    if preprocessing_type == "Eli":
        data = mne.io.read_raw_nirx(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}_config.hdr",
                                    verbose=False, preload=True)
    
        data_snirf = mne.io.read_raw_snirf(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}.snirf",
                                           optode_frame="mri", preload=True)
    
        aux_snirf = read_snirf_aux_data(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}.snirf",
                                        data_snirf)
    elif preprocessing_type == "Ben":
        data = read_raw_nirx(f"{curr_fnirs_data_folders[ii]}/", verbose=False, preload=True)
        data_snirf = mne.io.read_raw_snirf(f"{curr_fnirs_data_folders[ii]}/{curr_fnirs_data_folders[ii][-14:]}.snirf",
                                           optode_frame="mri", preload=True)
    
    
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
    #cue_dur = 1.8
    #data.annotations.onset = data.annotations.onset + cue_dur
    #data_snirf.annotations.onset = data_snirf.annotations.onset + cue_dur
    
    # ---------------------------------------------------------------
    # -------------               Preprocessing             ---------
    # ---------------------------------------------------------------
    
    if preprocessing_type == "Eli":
        events, event_dict = mne.events_from_annotations(data, verbose=False)
        
        if subject != "mild_master_5":
            this_sub_short_regression = True
        else:
            this_sub_short_regression = False
    
        raw_haemo_temp, null = preprocess_NIRX(data, data_snirf, event_dict,
                                               save=True,
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

    elif preprocessing_type == "Ben":
        events, event_dict = mne.events_from_annotations(data, verbose=False)

       # del event_dict['control']
        # Convert to optical density
        raw_od = optical_density(data_snirf)
         
        # Scalp Coupling Index, label bad channels
        sci = scalp_coupling_index(raw_od)

        # Add 'bads' to info
         
        raw_od.info['bads'] = list(compress(raw_od.ch_names,sci < 0.8))
        
        

        # raw_od = short_channel_regression(raw_od, max_dist=0.01)
        
        # Apply TDDR
        raw_od = temporal_derivative_distribution_repair(raw_od, verbose=False)
         
        # Resample to 3 Hz
        #raw_od.resample(3) # 10
         
        # Create separate object for block averages (will run short channel on these, but use short channels as a regressor in the GLM for betas)
        raw_od_regressed = short_channel_regression(raw_od.copy(), max_dist=0.01)
         
        #raw_haemo = beer_lambert_law(raw_od, ppf=0.1)
        raw_haemo = mne_modified_beer_lambert_law(raw_od_regressed) # TRYING ELIS FUNCTION
         
        # Filter data
        iir_params = dict({"order":2,"ftype":"butter","padlen":10000}) # 3
        raw_haemo = raw_haemo.filter(0.01, 0.1, iir_params=iir_params, method='iir', verbose=False) #0.01,0.3
        #raw_haemo.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005) # 0.01, 0.3, 0.2, 0.005
        
        raw_haemo_short = get_short_channels(raw_haemo)
        raw_haemo_filt = get_long_channels(raw_haemo)
        
        
        num_channels_removed[ii] = len(list(raw_haemo_filt.info['bads']))/2
        age[ii] = 2025 - data.info['subject_info']['birthday'][0]
        sex[ii] = data.info['subject_info']['sex']


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
    

    # for ichannel in range(len(chan_indices_good)):
        
    #     data_itd5_avg[ichannel,:] = np.nanmean(data_itd5[:,ichannel,:]  - np.nanmean(data_itd5[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_itd15_avg[ichannel,:] = np.nanmean(data_itd15[:,ichannel,:]  - np.nanmean(data_itd15[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5_avg[ichannel,:] = np.nanmean(data_ild5[:,ichannel,:]  -  np.nanmean(data_ild5[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild15_avg[ichannel,:] = np.nanmean(data_ild15[:,ichannel,:]  -  np.nanmean(data_ild15[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
        
    #     data_itd5_avg_hbr[ichannel,:] = np.nanmean(data_itd5_hbr[:,ichannel,:]  -  np.nanmean(data_itd5_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_itd15_avg_hbr[ichannel,:] = np.nanmean(data_itd15_hbr[:,ichannel,:]  -  np.nanmean(data_itd15_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5_avg_hbr[ichannel,:] = np.nanmean(data_ild5_hbr[:,ichannel,:]  -  np.nanmean(data_ild5_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild15_avg_hbr[ichannel,:] = np.nanmean(data_ild15_hbr[:,ichannel,:]  -  np.nanmean(data_ild15_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    
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
    

    #plot for this subject
    # cmin = -0.05
    # cmax = 0.15
    # time = np.linspace(tmin,tmax,num=n_timepoints)
    # fig, axes = plt.subplots(4, 2)
    # curr_ax = axes[0, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_itd5_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 50 HbO')
    
    # curr_ax = axes[1, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_itd15_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 500 HbO')
    
    # curr_ax = axes[2, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg HbO')
    
    # curr_ax = axes[3, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_ild15_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg + Mag HbO')
    
    # curr_ax = axes[0, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_itd5_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 50 HbR')
    
    # curr_ax = axes[1, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_itd15_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 500 HbR')
    
    # curr_ax = axes[2, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg HbR')
    
    # curr_ax = axes[3, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_ild15_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg + Mag HbR')
    
    # fig.suptitle(subject)
    # run a GLM to extract beta values for each condition

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
    
    

    design_matrix_hbo = make_first_level_design_matrix(raw_haemo_filt_for_glm.pick(picks='hbo'),
                                                        drift_model=None,
                                                        high_pass=0.01,  # Must be specified per experiment
                                                        hrf_model='spm',
                                                        stim_dur=glm_dur)
    # # add_regs=filtered_signals)

    design_matrix_hbo["Linear"] = np.arange(0, np.shape(design_matrix_hbo)[0])
    #design_matrix_hbo["ShortHbO"] = np.mean(raw_haemo_short.copy().pick(picks="hbo").get_data(), axis=0)

    #design_matrix_hbr["ShortHbR"] = np.mean(raw_haemo_short.copy().pick(picks="hbr").get_data(), axis=0)
    min_max_scaler = MinMaxScaler()
    X_minmax = min_max_scaler.fit_transform(design_matrix_hbo)
    design_matrix_min_max = pd.DataFrame(X_minmax, columns=design_matrix_hbo.columns.tolist())
    # if False:
    # # plotting optional
    #     fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    #     ax_img = plot_design_matrix(design_matrix_min_max, ax=ax1)
    #     plt.show()
    
    #     s = mne_nirs.experimental_design.create_boxcar(raw_haemo_filt, stim_dur=glm_dur)
    #     fig, ax2 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    #     plt.plot(raw_haemo_filt.times[0:2000], s[0:2000, 3])
    #     plt.plot(design_matrix_hbo['itd=500_az=0_mag=0'])
    #     plt.legend(["Stimulus", "Expected Response"])
    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Amplitude")
    #     plt.show()

    # print(f'running GLM for subject {ii + 1}')

    # pre-whiten
    raw_haemo_filt_for_glm._data = np.subtract(raw_haemo_filt_for_glm._data,
                                            np.nanmean(raw_haemo_filt_for_glm._data, axis=1)[:, np.newaxis])
    glm_est = run_glm(raw_haemo_filt_for_glm, design_matrix_hbo, noise_model='ar1')

    # record the glm est for each condition, for each subject
    # will adjust the beta values by the BH correction method

    glm_est_df = glm_est.pick(picks='data').to_dataframe()

    # # put into a larger array with all subjects data!
    subject_data_itd5_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'az_itd=5_az=0']['theta']
    subject_data_itd15_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'az_itd=15_az=0']['theta']
    subject_data_ild5_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'az_itd=0_az=5']['theta']
    subject_data_ild15_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'az_itd=0_az=15']['theta']


    # Append stats df to a larger dataframe
    

##############################
## Take mean during stim ####
############################

pfc_channels = []
stg_channels = []
# Take Means
index_stim_start = int(2*fs) # 8
index_stim_end = int(13.6*fs)
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

# Uncorrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd5_baselined.shape], names = names)
itd5_df = pd.DataFrame({'subject_data_itd5':subject_data_itd5_baselined.flatten()},index=index)['subject_data_itd5']
itd15_df = pd.DataFrame({'subject_data_itd15':subject_data_itd15_baselined.flatten()},index=index)['subject_data_itd15']
ild5_df = pd.DataFrame({'subject_data_ild5':subject_data_ild5_baselined.flatten()},index=index)['subject_data_ild5']
ild15_df = pd.DataFrame({'subject_data_ild15':subject_data_ild15_baselined.flatten()},index=index)['subject_data_ild15']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_uncorr_block_average_{masker_type}_masker.csv',index=True)

# Corrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd5_bh_corr.shape], names = names)
itd5_df = pd.DataFrame({'subject_data_itd5_bh_corr':subject_data_itd5_bh_corr.flatten()},index=index)['subject_data_itd5_bh_corr']
itd15_df = pd.DataFrame({'subject_data_itd15_bh_corr':subject_data_itd15_bh_corr.flatten()},index=index)['subject_data_itd15_bh_corr']
ild5_df = pd.DataFrame({'subject_data_ild5_bh_corr':subject_data_ild5_bh_corr.flatten()},index=index)['subject_data_ild5_bh_corr']
ild15_df = pd.DataFrame({'subject_data_ild15_bh_corr':subject_data_ild15_bh_corr.flatten()},index=index)['subject_data_ild15_bh_corr']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_bh_corr_block_average_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBO
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd5.shape], names = names)
itd5_df = pd.DataFrame({'mean_during_stim_itd5':mean_during_stim_itd5.flatten()},index=index)['mean_during_stim_itd5']
itd15_df = pd.DataFrame({'mean_during_stim_itd15':mean_during_stim_itd15.flatten()},index=index)['mean_during_stim_itd15']
ild5_df = pd.DataFrame({'mean_during_stim_ild5':mean_during_stim_ild5.flatten()},index=index)['mean_during_stim_ild5']
ild15_df = pd.DataFrame({'mean_during_stim_ild15':mean_during_stim_ild15.flatten()},index=index)['mean_during_stim_ild15']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBR
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd5_hbr.shape], names = names)
itd5_df = pd.DataFrame({'mean_during_stim_itd5_hbr':mean_during_stim_itd5_hbr.flatten()},index=index)['mean_during_stim_itd5_hbr']
itd15_df = pd.DataFrame({'mean_during_stim_itd15_hbr':mean_during_stim_itd15_hbr.flatten()},index=index)['mean_during_stim_itd15_hbr']
ild5_df = pd.DataFrame({'mean_during_stim_ild5_hbr':mean_during_stim_ild5_hbr.flatten()},index=index)['mean_during_stim_ild5_hbr']
ild15_df = pd.DataFrame({'mean_during_stim_ild15_hbr':mean_during_stim_ild15_hbr.flatten()},index=index)['mean_during_stim_ild15_hbr']
z = pd.concat([itd5_df,itd15_df,ild5_df,ild15_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_{masker_type}_masker_hbr.csv',index=True)
# ---------------------------------------------------------------
# -----  Correlation Between BH Coeffient and Task Beta ---------
# ---------------------------------------------------------------
# determine if channels with large task-related activation are the ones with large BH activation
# flat_BH_range = np.ravel(range_BH_response)
# flat_speech_beta = np.ravel(subject_data_speech_GLM)
# flat_noise_beta = np.ravel(subject_data_noise_GLM)
# flat_control_beta = np.ravel(subject_data_control_GLM)
#
# flat_combined = np.array([flat_BH_range, flat_speech_beta])
# cols_with_nan = np.any(np.isnan(flat_combined), axis=0)
# flat_combined_clean = flat_combined[:, ~cols_with_nan]
#
# pearson_corr_bh_speech, p_value_bh_speech = stats.pearsonr(flat_combined_clean[0, :], flat_combined_clean[1, :])
#
# plt.plot(flat_combined_clean[0, :], flat_combined_clean[1, :], 'o'); plt.show()

# ---------------------------------------------------------------
# -----------------     Subject Averaging                ---------
# ---------------------------------------------------------------
# # for each subject, take the beta values and compute a mean and standard error
subject_data_itd5_GLM_mean = np.nanmean(1e6*subject_data_itd5_GLM, axis=0)
subject_data_itd15_GLM_mean = np.nanmean(1e6*subject_data_itd15_GLM, axis=0)
subject_data_ild5_GLM_mean = np.nanmean(1e6*subject_data_ild5_GLM, axis=0)
subject_data_ild15_GLM_mean = np.nanmean(1e6*subject_data_ild15_GLM, axis=0)

subject_data_itd5_GLM_mean[np.isnan(subject_data_itd5_GLM_mean)] = 0
subject_data_itd15_GLM_mean[np.isnan(subject_data_itd15_GLM_mean)] = 0
subject_data_ild5_GLM_mean[np.isnan(subject_data_ild5_GLM_mean)] = 0
subject_data_ild15_GLM_mean[np.isnan(subject_data_ild15_GLM_mean)] = 0


subject_data_itd5_GLM_std = np.nanstd(1e6*subject_data_itd5_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_itd15_GLM_std = np.nanstd(1e6*subject_data_itd15_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_ild5_GLM_std = np.nanstd(1e6*subject_data_ild5_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_ild15_GLM_std = np.nanstd(1e6*subject_data_ild15_GLM, axis=0) / np.sqrt(n_subjects)

# # do the same for breath hold correction now...
# subject_data_itd5_GLM_bh_corr_mean = np.nanmean(subject_data_itd5_GLM_bh_corr, axis=0)
# subject_data_itd15_GLM_bh_corr_mean = np.nanmean(subject_data_itd15_GLM_bh_corr, axis=0)
# subject_data_ild5_GLM_bh_corr_mean = np.nanmean(subject_data_ild5_GLM_bh_corr, axis=0)
# subject_data_ild15_GLM_bh_corr_mean = np.nanmean(subject_data_ild15_GLM_bh_corr, axis=0)

# subject_data_itd5_GLM_bh_corr_std = np.nanstd(subject_data_itd5_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_itd15_GLM_bh_corr_std = np.nanstd(subject_data_itd15_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_ild5_GLM_bh_corr_std = np.nanstd(subject_data_ild5_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_ild15_GLM_bh_corr_std = np.nanstd(subject_data_ild15_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)




# ---------------------------------------------------------------
# -----------------     PLotting GLM Averages           ---------
# ---------------------------------------------------------------
# caxis_lim = 8e-8

# # caxis_lim = np.nanmean([np.nanmean(np.abs(subject_data_itd5_GLM_mean)),
# #             np.nanmean(np.abs(subject_data_itd15_GLM_mean)),
# #             np.nanmean(np.abs(subject_data_ild5_GLM_mean)),
# # #             np.nanmean(np.abs(subject_data_ild15_GLM_mean))])
# fig, axes = plt.subplots(int(n_long_channels/4),4)
# fig.set_figwidth(4)
# fig.set_figheight(8)


# # Plot error bars
# for ichannel, curr_axes in enumerate(axes.reshape(-1)):
#     # itd5
#     curr_axes.errorbar(1, subject_data_itd5_GLM_mean[ichannel],np.nanstd(subject_data_itd5_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_itd5_GLM, axis=0) - 1)),fmt='o')    
#     # itd15
#     curr_axes.errorbar(2, subject_data_itd15_GLM_mean[ichannel], 
#                         np.nanstd(subject_data_itd15_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_itd15_GLM, axis=0)- 1) ),
#                         fmt='o')
#     # ild5
#     curr_axes.errorbar(3, subject_data_ild5_GLM_mean[ichannel],
#                         np.nanstd(subject_data_ild5_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_ild5_GLM, axis=0)- 1) ),
#                         fmt='o')
#     # ild15 
#     curr_axes.errorbar(4, subject_data_ild15_GLM_mean[ichannel], 
#                         np.nanstd(subject_data_ild15_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_ild15_GLM, axis=0)- 1) ),
#                         fmt='o')
    
#     curr_axes.set_ylim((-0.1,0.1)) 
#     curr_axes.set_xlabel('Condition')
#     curr_axes.set_title(channel_names[ichannel])
#     curr_axes.set_xticks([1,2,3,4])
#     if ichannel >= 12:
#         curr_axes.set_xticklabels(["50 us ITD","500 us ITD","5 deg ILD","5 deg ILD + Mag"])
#     if ichannel == 4 or ichannel == 8:
#         curr_axes.set_ylabel('Beta', usetex=False)

# plt.subplots_adjust(top=.9, right=0.999, left= 0.1, bottom = 0.07)
# plt.savefig(f'C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\Plots\\Beta_dur_{glm_dur}_{masker_type}_masker.png')
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
plt.show()


caxis_lim = 0.15


groups_single_chroma = dict(
    Left_Hemisphere =  picks_pair_to_idx(raw_haemo_filt.copy().pick(picks='hbo'), [[1,8]], on_missing='warning'),
    Right_Hemisphere =  picks_pair_to_idx(raw_haemo_filt.copy().pick(picks='hbo'), [[2,7]], on_missing='warning'))

fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1, 5)
plot_glm_group_topo(raw_haemo_filt.copy().pick(picks="hbo").pick(groups_single_chroma['Left_Hemisphere']),
                    subject_data_itd5_GLM_mean,                ####Change here####
                    colorbar=False, axes=ax1,
                    vlim=(-caxis_lim, caxis_lim), cmap='summer')

plt.show()


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


# ---------------------------------------------------------------
# -----------------     Compare Mean HbO and GLM Values ---------
# ---------------------------------------------------------------
# fig, axes = plt.subplots(2, 7)
# fig.set_figwidth(16)
# fig.set_figheight(8)
# caxis_lim = 11e-8
# for ichannel, curr_axes in enumerate(axes.reshape(-1)):
#     # itd5
#     curr_axes.scatter(mean_during_stim_itd5[:,ichannel],subject_data_itd5_GLM[:,ichannel])
    
#     # itd15
#     curr_axes.scatter(mean_during_stim_itd15[:,ichannel],subject_data_itd15_GLM[:,ichannel])
    
#     # ild5
#     curr_axes.scatter(mean_during_stim_ild5[:,ichannel],subject_data_ild5_GLM[:,ichannel])
    
#     # ild15
#     curr_axes.scatter(mean_during_stim_ild15[:,ichannel],subject_data_ild15_GLM[:,ichannel])
    
#     curr_axes.set_xlabel(r'Mean $\Delta$Hb ($\mu$M) during stim.', usetex=False)
#     curr_axes.set_ylabel('Beta')
# plt.savefig(f'C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\Plots\\Beta_vs_HbO_dur_{glm_dur}_{masker_type}_masker.png')



# # Uncorrected GLM
# names = ['S','Channel']
# index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd5_GLM.shape], names = names)
# itd5_df_GLM = pd.DataFrame({'subject_data_itd5_GLM':subject_data_itd5_GLM.flatten()},index=index)['subject_data_itd5_GLM']
# itd15_df_GLM = pd.DataFrame({'subject_data_itd15_GLM':subject_data_itd15_GLM.flatten()},index=index)['subject_data_itd15_GLM']
# ild5_df_GLM = pd.DataFrame({'subject_data_ild5_GLM':subject_data_ild5_GLM.flatten()},index=index)['subject_data_ild5_GLM']
# ild15_df_GLM = pd.DataFrame({'subject_data_ild15_GLM':subject_data_ild15_GLM.flatten()},index=index)['subject_data_ild15_GLM']
# z = pd.concat([itd5_df_GLM,itd15_df_GLM,ild5_df_GLM,ild15_df_GLM], ignore_index=True,axis=1)
# z.to_csv(f'all_subjects_uncorr_GLM_{masker_type}_masker.csv',index=True)


# # Corrected GLM

# names = ['S','Channel']
# index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd5_GLM.shape], names = names)
# itd5_df_GLM_bh_corr = pd.DataFrame({'subject_data_itd5_GLM_bh_corr':subject_data_itd5_GLM_bh_corr.flatten()},index=index)['subject_data_itd5_GLM_bh_corr']
# itd15_df_GLM_bh_corr = pd.DataFrame({'subject_data_itd15_GLM_bh_corr':subject_data_itd15_GLM_bh_corr.flatten()},index=index)['subject_data_itd15_GLM_bh_corr']
# ild5_df_GLM_bh_corr = pd.DataFrame({'subject_data_ild5_GLM_bh_corr':subject_data_ild5_GLM_bh_corr.flatten()},index=index)['subject_data_ild5_GLM_bh_corr']
# ild15_df_GLM_bh_corr = pd.DataFrame({'subject_data_ild15_GLM_bh_corr':subject_data_ild15_GLM_bh_corr.flatten()},index=index)['subject_data_ild15_GLM_bh_corr']
# z = pd.concat([itd5_df_GLM_bh_corr,itd15_df_GLM_bh_corr,ild5_df_GLM_bh_corr,ild15_df_GLM_bh_corr], ignore_index=True,axis=1)
# z.to_csv(f'all_subjects_bh_corr_GLM_{masker_type}_masker.csv',index=True)