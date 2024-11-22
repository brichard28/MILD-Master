# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:34:12 2024

@author: benri
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 19:47:20 2024

@author: benri
"""

import numpy as np
import mne
import pickle
import matplotlib
import math
# matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import os
import pandas as pd
from collections import defaultdict

mne.set_config('MNE_BROWSER_BACKEND', 'qt')
#from nirx_movement import mark_aux_movement_bad
import mne_nirs
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)
#import pandas as pd
from nilearn.plotting import plot_design_matrix
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.io.snirf import read_snirf_aux_data
from sklearn.preprocessing import StandardScaler, MinMaxScaler
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
user = 'Laptop'
if user == 'Laptop':
    data_root = 'C:/Users/benri/Downloads/'

else:
    data_root = '/home/ben/Nextcloud/data/nirs/data/'
    

# all_fnirs_data_folders = [data_root + '2024-10-30/2024-10-30_001',
#                           data_root + '2024-10-31/2024-10-31_001']

all_fnirs_data_folders = [data_root +'2024-10-30/2024-10-30_001',
                          data_root +'2024-10-31/2024-10-31_001',
                          data_root + '2024-11-01/2024-11-01_001',
                          data_root + '2024-11-05/2024-11-05_001',
                          data_root + '2024-11-07/2024-11-07_001',
                          data_root + '2024-11-08/2024-11-08_001',
                          data_root + '2024-11-11/2024-11-11_001',
                          data_root + '2024-11-12/2024-11-12_001',
                          data_root + '2024-11-14/2024-11-14_001']



# All subject IDs
subject_ID = ['dense_nirs_2','dense_nirs_3','dense_nirs_4','dense_nirs_5','dense_nirs_6','dense_nirs_7','dense_nirs_8','dense_nirs_9','dense_nirs_10']
# The subjects we would like to run right now
curr_subject_ID = ['dense_nirs_10']

curr_folder_indices = [index for index, element in enumerate(subject_ID) if np.isin(element,curr_subject_ID)]
curr_fnirs_data_folders = [all_fnirs_data_folders[i] for i in curr_folder_indices]


masker_type = 'speech' # type of masker to analyze on this run

glm_dur = 11.6

n_subjects = len(curr_subject_ID)

n_long_channels = 97
fs = 5.1
tmin, tmax = -5, 23
n_timepoints = math.floor((tmax - tmin)*fs) + 1
task_type = 'Ben_SvN'



# set up the arrays to hold all subject data
subject_data_itd50 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500 = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5deg = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5degMag = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5deg_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5degMag_hbr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5deg_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5degMag_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5deg_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5degMag_hbr_baselined = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)


subject_data_itd50_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_itd500_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5deg_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)
subject_data_ild5degMag_bh_corr = np.full((n_subjects, n_long_channels, n_timepoints), np.nan)

subject_data_itd50_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd500_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5deg_GLM = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5degMag_GLM = np.full((n_subjects, n_long_channels), np.nan)



# put into a larger array with all subjects data!
subject_data_itd50_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_itd500_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5deg_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)
subject_data_ild5degMag_GLM_bh_corr = np.full((n_subjects, n_long_channels), np.nan)

num_channels_removed = np.full(n_subjects, np.nan)
age = np.full(n_subjects, np.nan)
sex = np.full(n_subjects, np.nan)

range_BH_response = np.zeros((n_subjects, n_long_channels))

all_evokeds = defaultdict(list)
subject_info = []
# loop through all subjects and all sessions (takes a while)
preprocessing_type = "Eli"


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
    data.annotations.rename({'4.0':	'itd=500_az=0_mag=0',
                                '5.0':	'itd=50_az=0_mag=0',
                                '6.0':	'itd=0_az=5_mag=0',
                                '7.0':	'itd=0_az=5_mag=1',
                                '8.0':	'itd=50_az=0_mag=0',
                                '9.0':	'itd=500_az=0_mag=0',
                                '10.0':	'itd=0_az=5_mag=1',
                                '11.0':	'itd=0_az=5_mag=0'})
    data_snirf.annotations.rename({'4':	'itd=500_az=0_mag=0',
                                '5':	'itd=50_az=0_mag=0',
                                '6':	'itd=0_az=5_mag=0',
                                '7':	'itd=0_az=5_mag=1',
                                '8':	'itd=50_az=0_mag=0',
                                '9':	'itd=500_az=0_mag=0',
                                '10':	'itd=0_az=5_mag=1',
                                '11':	'itd=0_az=5_mag=0'})
    
    # Trying out shifting the trigger to the stim onset (rather than cue)
    cue_dur = [2, 2] #[1.85, 2.05]
    #data.annotations.onset = data.annotations.onset + cue_dur[ii]
    #data_snirf.annotations.onset = data_snirf.annotations.onset + cue_dur[ii]
    
    # ---------------------------------------------------------------
    # -------------               Preprocessing             ---------
    # ---------------------------------------------------------------
    
    if preprocessing_type == "Eli":
        events, event_dict = mne.events_from_annotations(data, verbose=False)
        
        #fig = mne.viz.plot_events(events, event_id=event_dict, sfreq=data.info["sfreq"])
    
        if subject == 'dense_nirs_10':
            # Sources for short channels, 2,5,7,10,13,16,22,23
            channels_to_drop = ['S1_D22 760',
             'S1_D22 850',
             'S3_D24 760',
             'S3_D24 850',
             'S4_D25 760',
             'S4_D25 850',
             'S6_D27 760',
             'S6_D27 850',
             'S8_D29 760',
             'S8_D29 850',
             'S9_D30 760',
             'S9_D30 850',
             'S11_D32 760',
             'S11_D32 850',
             'S12_D33 760',
             'S12_D33 850',
             'S14_D35 760',
             'S14_D35 850',
             'S15_D36 760',
             'S15_D36 850',
             'S17_D38 760',
             'S17_D38 850',
             'S18_D39 760',
             'S18_D39 850',
             'S19_D40 760',
             'S19_D40 850',
             'S20_D41 760',
             'S20_D41 850',
             'S21_D42 760',
             'S21_D42 850',
             'S24_D45 760',
             'S24_D45 850']
            
            data.drop_channels(channels_to_drop)
            data_snirf.drop_channels(channels_to_drop)
            
        
        raw_haemo_temp, null = preprocess_NIRX(data, data_snirf, event_dict,
                                               save=True,
                                               savename=save_dir + f'{subject}_{task_type}_preproc_nirs.fif',
                                               plot_steps=True,
                                               crop=False, crop_low=0, crop_high=0,
                                               events_modification=False, reject=True,
                                               short_regression=True, events_from_snirf=False,
                                               drop_short=False, negative_enhancement=False,
                                               snr_thres=3, sci_thres=0.8, filter_type='iir', filter_limits=[0.01,0.3])
    
        raw_haemo_short = get_short_channels(raw_haemo_temp)
        raw_haemo_filt = get_long_channels(raw_haemo_temp)
        
        
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
        raw_haemo = raw_haemo.filter(0.01, 0.3, iir_params=iir_params, method='iir', verbose=False) #0.01,0.3
        #raw_haemo.filter(0.01, 0.2, h_trans_bandwidth=0.2, l_trans_bandwidth=0.005) # 0.01, 0.3, 0.2, 0.005
        
        raw_haemo_short = get_short_channels(raw_haemo)
        raw_haemo_filt = get_long_channels(raw_haemo)
        
        
        num_channels_removed[ii] = len(list(raw_haemo_filt.info['bads']))/2
        age[ii] = 2024 - data.info['subject_info']['birthday'][0]
        sex[ii] = data.info['subject_info']['sex']


    # ---------------------------------------------------------------
    # -------------               Epoching                  ---------
    # ---------------------------------------------------------------
    reject_criteria = dict(hbo=0.2)#5e-6
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
    
    # Plot epochs at each sensor
    conditions = ['itd=50_az=0_mag=0','itd=500_az=0_mag=0','itd=0_az=5_mag=0','itd=0_az=5_mag=1']  
    n_conditions = len(conditions)
    evoked_hbo = np.zeros((n_conditions,), dtype=object)
    evoked_hbo_error = np.zeros((n_conditions,), dtype=object)
    evoked_hbr = np.zeros((n_conditions,), dtype=object)
    evoked_hbr_error = np.zeros((n_conditions,), dtype=object)
    vlim = 0.2
    n_conditions = 1
    for i, cond in enumerate(['itd=50_az=0_mag=0']):
        evoked_hbo[i] = epochs['itd=50_az=0_mag=0','itd=500_az=0_mag=0','itd=0_az=5_mag=0','itd=0_az=5_mag=1'].copy().average(picks='hbo')
        evoked_hbo_error[i] = epochs['itd=50_az=0_mag=0','itd=500_az=0_mag=0','itd=0_az=5_mag=0','itd=0_az=5_mag=1'].copy().standard_error(picks='hbo')
        evoked_hbr[i] = epochs['itd=50_az=0_mag=0','itd=500_az=0_mag=0','itd=0_az=5_mag=0','itd=0_az=5_mag=1'].copy().average(picks='hbr')
        evoked_hbr_error[i] = epochs['itd=50_az=0_mag=0','itd=500_az=0_mag=0','itd=0_az=5_mag=0','itd=0_az=5_mag=1'].copy().standard_error(picks='hbr')
        fig = plt.figure(figsize=(8, 5), dpi=200)
        fig = plot_nirs_evoked_error(fig, evoked_hbo[i], evoked_hbo_error[i],
                                     colour='r', ylim=[-vlim, vlim])
        fig = plot_nirs_evoked_error(fig, evoked_hbr[i], evoked_hbr_error[i],
                             colour='b', ylim=[-vlim, vlim], add_legend=True)
    


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
    data_itd50 = epochs["itd=50_az=0_mag=0"].get_data(picks='hbo')
    data_itd500 = epochs["itd=500_az=0_mag=0"].get_data(picks='hbo')
    data_ild5deg = epochs["itd=0_az=5_mag=0"].get_data(picks='hbo')
    data_ild5degMag = epochs["itd=0_az=5_mag=1"].get_data(picks='hbo')
    
    data_itd50_hbr = epochs["itd=50_az=0_mag=0"].get_data(picks='hbr')
    data_itd500_hbr = epochs["itd=500_az=0_mag=0"].get_data(picks='hbr')
    data_ild5deg_hbr = epochs["itd=0_az=5_mag=0"].get_data(picks='hbr')
    data_ild5degMag_hbr = epochs["itd=0_az=5_mag=1"].get_data(picks='hbr')
    
    # ---------------------------------------------------------------
    # -----------------    Baselining and Averaging         ---------
    # ---------------------------------------------------------------
    data_itd50_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_itd500_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_ild5deg_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)
    data_ild5degMag_avg = np.full((len(chan_indices_good_hbo), n_timepoints), np.nan)

    data_itd50_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_itd500_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_ild5deg_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    data_ild5degMag_avg_hbr= np.full((len(chan_indices_good_hbr), n_timepoints), np.nan)
    

    # for ichannel in range(len(chan_indices_good)):
        
    #     data_itd50_avg[ichannel,:] = np.nanmean(data_itd50[:,ichannel,:]  - np.nanmean(data_itd50[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_itd500_avg[ichannel,:] = np.nanmean(data_itd500[:,ichannel,:]  - np.nanmean(data_itd500[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5deg_avg[ichannel,:] = np.nanmean(data_ild5deg[:,ichannel,:]  -  np.nanmean(data_ild5deg[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5degMag_avg[ichannel,:] = np.nanmean(data_ild5degMag[:,ichannel,:]  -  np.nanmean(data_ild5degMag[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
        
    #     data_itd50_avg_hbr[ichannel,:] = np.nanmean(data_itd50_hbr[:,ichannel,:]  -  np.nanmean(data_itd50_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_itd500_avg_hbr[ichannel,:] = np.nanmean(data_itd500_hbr[:,ichannel,:]  -  np.nanmean(data_itd500_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5deg_avg_hbr[ichannel,:] = np.nanmean(data_ild5deg_hbr[:,ichannel,:]  -  np.nanmean(data_ild5deg_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    #     data_ild5degMag_avg_hbr[ichannel,:] = np.nanmean(data_ild5degMag_hbr[:,ichannel,:]  -  np.nanmean(data_ild5degMag_hbr[:,ichannel,0:int(-1*tmin*fs)], axis = (0,1)), axis=0)
    
    for ichannel in range(len(chan_indices_good_hbo)):
        
        data_itd50_avg[ichannel,:] = np.nanmean(data_itd50[:,ichannel,:], axis=0)
        data_itd500_avg[ichannel,:] = np.nanmean(data_itd500[:,ichannel,:], axis=0)
        data_ild5deg_avg[ichannel,:] = np.nanmean(data_ild5deg[:,ichannel,:], axis=0)
        data_ild5degMag_avg[ichannel,:] = np.nanmean(data_ild5degMag[:,ichannel,:], axis=0)
        
    for ichannel in range(len(chan_indices_good_hbr)):
        data_itd50_avg_hbr[ichannel,:] = np.nanmean(data_itd50_hbr[:,ichannel,:], axis=0)
        data_itd500_avg_hbr[ichannel,:] = np.nanmean(data_itd500_hbr[:,ichannel,:], axis=0)
        data_ild5deg_avg_hbr[ichannel,:] = np.nanmean(data_ild5deg_hbr[:,ichannel,:], axis=0)
        data_ild5degMag_avg_hbr[ichannel,:] = np.nanmean(data_ild5degMag_hbr[:,ichannel,:], axis=0)
    
    # need to mark the indices where the good channels are!

    # put into a larger array with all subjects data!
    subject_data_itd50[ii, chan_indices_good_hbo, :] = 1e6*data_itd50_avg.copy()
    subject_data_itd500[ii, chan_indices_good_hbo, :] = 1e6*data_itd500_avg.copy()
    subject_data_ild5deg[ii, chan_indices_good_hbo, :] = 1e6*data_ild5deg_avg.copy()
    subject_data_ild5degMag[ii, chan_indices_good_hbo, :] = 1e6*data_ild5degMag_avg.copy()
  
    subject_data_itd50_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_itd50_avg_hbr.copy()
    subject_data_itd500_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_itd500_avg_hbr.copy()
    subject_data_ild5deg_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_ild5deg_avg_hbr.copy()
    subject_data_ild5degMag_hbr[ii, chan_indices_good_hbr, :] = 1e6*data_ild5degMag_avg_hbr.copy()
    
    subject_data_itd50_baselined = subject_data_itd50.copy()
    subject_data_itd500_baselined = subject_data_itd500.copy()
    subject_data_ild5deg_baselined = subject_data_ild5deg.copy()
    subject_data_ild5degMag_baselined = subject_data_ild5degMag.copy()
    
    subject_data_itd50_hbr_baselined = subject_data_itd50_hbr.copy()
    subject_data_itd500_hbr_baselined = subject_data_itd500_hbr.copy()
    subject_data_ild5deg_hbr_baselined = subject_data_ild5deg_hbr.copy()
    subject_data_ild5degMag_hbr_baselined = subject_data_ild5degMag_hbr.copy()
    

    #plot for this subject
    # cmin = -0.05
    # cmax = 0.15
    # time = np.linspace(tmin,tmax,num=n_timepoints)
    # fig, axes = plt.subplots(4, 2)
    # curr_ax = axes[0, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_itd50_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 50 HbO')
    
    # curr_ax = axes[1, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_itd500_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 500 HbO')
    
    # curr_ax = axes[2, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5deg_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg HbO')
    
    # curr_ax = axes[3, 0]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5degMag_avg[:, :]), 'r')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg + Mag HbO')
    
    # curr_ax = axes[0, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_itd50_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 50 HbR')
    
    # curr_ax = axes[1, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_itd500_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ITD 500 HbR')
    
    # curr_ax = axes[2, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5deg_avg_hbr[:, :]), 'b')
    # curr_ax.set_ylim([cmin,cmax])
    # curr_ax.axvline(x=0,color='k')
    # curr_ax.set_title('ILD 5 deg HbR')
    
    # curr_ax = axes[3, 1]
    # curr_ax.plot(time,np.transpose(1e6*data_ild5degMag_avg_hbr[:, :]), 'b')
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
    raw_haemo_filt_for_glm.drop_channels(epochs.copy().info['bads'])

    raw_haemo_filt_for_glm.annotations.set_durations(glm_dur)
    
    

    design_matrix_hbo = make_first_level_design_matrix(raw_haemo_filt_for_glm.pick(picks='hbo'),
                                                        drift_model=None,
                                                        high_pass=0.01,  # Must be specified per experiment
                                                        hrf_model='spm',
                                                        stim_dur=raw_haemo_filt.annotations.duration)
    # # add_regs=filtered_signals)

    design_matrix_hbo["Linear"] = np.arange(0, np.shape(design_matrix_hbo)[0])
    #design_matrix_hbo["ShortHbO"] = np.mean(raw_haemo_short.copy().pick(picks="hbo").get_data(), axis=0)

    #design_matrix_hbr["ShortHbR"] = np.mean(raw_haemo_short.copy().pick(picks="hbr").get_data(), axis=0)
    min_max_scaler = MinMaxScaler()
    X_minmax = min_max_scaler.fit_transform(design_matrix_hbo)
    design_matrix_min_max = pd.DataFrame(X_minmax, columns=design_matrix_hbo.columns.tolist())
    if False:
    # plotting optional
        fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
        ax_img = plot_design_matrix(design_matrix_min_max, ax=ax1)
        plt.show()
    
        s = mne_nirs.experimental_design.create_boxcar(raw_haemo_filt, stim_dur=glm_dur)
        fig, ax2 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
        plt.plot(raw_haemo_filt.times[0:2000], s[0:2000, 3])
        plt.plot(design_matrix_hbo['itd=500_az=0_mag=0'])
        plt.legend(["Stimulus", "Expected Response"])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.show()

    # print(f'running GLM for subject {ii + 1}')

    # pre-whiten
    raw_haemo_filt_for_glm._data = np.subtract(raw_haemo_filt_for_glm._data,
                                            np.mean(raw_haemo_filt_for_glm._data, axis=1)[:, np.newaxis])
    glm_est = run_glm(raw_haemo_filt_for_glm, design_matrix_hbo, noise_model='ar1')

    # record the glm est for each condition, for each subject
    # will adjust the beta values by the BH correction method

    glm_est_df = glm_est.pick(picks='data', exclude='bads').to_dataframe()

    # # put into a larger array with all subjects data!
    subject_data_itd50_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'itd=50_az=0_mag=0']['theta']
    subject_data_itd500_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'itd=500_az=0_mag=0']['theta']
    subject_data_ild5deg_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'itd=0_az=5_mag=0']['theta']
    subject_data_ild5degMag_GLM[ii, chan_indices_good_hbo] = glm_est_df.loc[glm_est_df['Condition'] == 'itd=0_az=5_mag=1']['theta']



##############################
## Take mean during stim ####
############################

pfc_channels = []
stg_channels = []
# Take Means
index_stim_start = int(2*fs) # 8
index_stim_end = int(6.8*fs)
# ITD50
mean_during_stim_itd50 = np.nanmean(subject_data_itd50_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd50_hbr = np.nanmean(subject_data_itd50_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# ITD500
mean_during_stim_itd500 = np.nanmean(subject_data_itd500_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_itd500_hbr = np.nanmean(subject_data_itd500_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

# ild5deg
mean_during_stim_ild5deg = np.nanmean(subject_data_ild5deg_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild5deg_hbr = np.nanmean(subject_data_ild5deg_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

#ild5degMag
mean_during_stim_ild5degMag = np.nanmean(subject_data_ild5degMag_baselined[:,:,index_stim_start:index_stim_end], axis=2)
mean_during_stim_ild5degMag_hbr = np.nanmean(subject_data_ild5degMag_hbr_baselined[:,:,index_stim_start:index_stim_end], axis=2)

## Save breath uncorrected and corrected GLM data

# Uncorrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd50_baselined.shape], names = names)
itd50_df = pd.DataFrame({'subject_data_itd50':subject_data_itd50_baselined.flatten()},index=index)['subject_data_itd50']
itd500_df = pd.DataFrame({'subject_data_itd500':subject_data_itd500_baselined.flatten()},index=index)['subject_data_itd500']
ild5deg_df = pd.DataFrame({'subject_data_ild5deg':subject_data_ild5deg_baselined.flatten()},index=index)['subject_data_ild5deg']
ild5degMag_df = pd.DataFrame({'subject_data_ild5degMag':subject_data_ild5degMag_baselined.flatten()},index=index)['subject_data_ild5degMag']
z = pd.concat([itd50_df,itd500_df,ild5deg_df,ild5degMag_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_uncorr_block_average_{masker_type}_masker.csv',index=True)

# Corrected block averages
names = ['S','Channel','Time_Index']
index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd50_bh_corr.shape], names = names)
itd50_df = pd.DataFrame({'subject_data_itd50_bh_corr':subject_data_itd50_bh_corr.flatten()},index=index)['subject_data_itd50_bh_corr']
itd500_df = pd.DataFrame({'subject_data_itd500_bh_corr':subject_data_itd500_bh_corr.flatten()},index=index)['subject_data_itd500_bh_corr']
ild5deg_df = pd.DataFrame({'subject_data_ild5deg_bh_corr':subject_data_ild5deg_bh_corr.flatten()},index=index)['subject_data_ild5deg_bh_corr']
ild5degMag_df = pd.DataFrame({'subject_data_ild5degMag_bh_corr':subject_data_ild5degMag_bh_corr.flatten()},index=index)['subject_data_ild5degMag_bh_corr']
z = pd.concat([itd50_df,itd500_df,ild5deg_df,ild5degMag_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_bh_corr_block_average_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBO
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd50.shape], names = names)
itd50_df = pd.DataFrame({'mean_during_stim_itd50':mean_during_stim_itd50.flatten()},index=index)['mean_during_stim_itd50']
itd500_df = pd.DataFrame({'mean_during_stim_itd500':mean_during_stim_itd500.flatten()},index=index)['mean_during_stim_itd500']
ild5deg_df = pd.DataFrame({'mean_during_stim_ild5deg':mean_during_stim_ild5deg.flatten()},index=index)['mean_during_stim_ild5deg']
ild5degMag_df = pd.DataFrame({'mean_during_stim_ild5degMag':mean_during_stim_ild5degMag.flatten()},index=index)['mean_during_stim_ild5degMag']
z = pd.concat([itd50_df,itd500_df,ild5deg_df,ild5degMag_df], ignore_index=True,axis=1)
z.to_csv(f'all_subjects_mean_during_stim_{masker_type}_masker.csv',index=True)


# Uncorrected block average means HBR
names = ['S','Channel']
index = pd.MultiIndex.from_product([range(s) for s in mean_during_stim_itd50_hbr.shape], names = names)
itd50_df = pd.DataFrame({'mean_during_stim_itd50_hbr':mean_during_stim_itd50_hbr.flatten()},index=index)['mean_during_stim_itd50_hbr']
itd500_df = pd.DataFrame({'mean_during_stim_itd500_hbr':mean_during_stim_itd500_hbr.flatten()},index=index)['mean_during_stim_itd500_hbr']
ild5deg_df = pd.DataFrame({'mean_during_stim_ild5deg_hbr':mean_during_stim_ild5deg_hbr.flatten()},index=index)['mean_during_stim_ild5deg_hbr']
ild5degMag_df = pd.DataFrame({'mean_during_stim_ild5degMag_hbr':mean_during_stim_ild5degMag_hbr.flatten()},index=index)['mean_during_stim_ild5degMag_hbr']
z = pd.concat([itd50_df,itd500_df,ild5deg_df,ild5degMag_df], ignore_index=True,axis=1)
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
subject_data_itd50_GLM_mean = np.nanmean(1e6*subject_data_itd50_GLM, axis=0)
subject_data_itd500_GLM_mean = np.nanmean(1e6*subject_data_itd500_GLM, axis=0)
subject_data_ild5deg_GLM_mean = np.nanmean(1e6*subject_data_ild5deg_GLM, axis=0)
subject_data_ild5degMag_GLM_mean = np.nanmean(1e6*subject_data_ild5degMag_GLM, axis=0)

subject_data_itd50_GLM_std = np.nanstd(1e6*subject_data_itd50_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_itd500_GLM_std = np.nanstd(1e6*subject_data_itd500_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_ild5deg_GLM_std = np.nanstd(1e6*subject_data_ild5deg_GLM, axis=0) / np.sqrt(n_subjects)
subject_data_ild5degMag_GLM_std = np.nanstd(1e6*subject_data_ild5degMag_GLM, axis=0) / np.sqrt(n_subjects)

# # do the same for breath hold correction now...
# subject_data_itd50_GLM_bh_corr_mean = np.nanmean(subject_data_itd50_GLM_bh_corr, axis=0)
# subject_data_itd500_GLM_bh_corr_mean = np.nanmean(subject_data_itd500_GLM_bh_corr, axis=0)
# subject_data_ild5deg_GLM_bh_corr_mean = np.nanmean(subject_data_ild5deg_GLM_bh_corr, axis=0)
# subject_data_ild5degMag_GLM_bh_corr_mean = np.nanmean(subject_data_ild5degMag_GLM_bh_corr, axis=0)

# subject_data_itd50_GLM_bh_corr_std = np.nanstd(subject_data_itd50_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_itd500_GLM_bh_corr_std = np.nanstd(subject_data_itd500_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_ild5deg_GLM_bh_corr_std = np.nanstd(subject_data_ild5deg_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)
# subject_data_ild5degMag_GLM_bh_corr_std = np.nanstd(subject_data_ild5degMag_GLM_bh_corr, axis=0) / np.sqrt(n_subjects)



# ---------------------------------------------------------------
# -----------------     PLotting Block Averages (PFC)   ---------
# ---------------------------------------------------------------


channel_names = [this_chan_hbo.replace(' hbo','') for this_chan_hbo in chan_hbo]
# ymin = -5e-2
# ymax = 20e-2
# fig, axes = plt.subplots(n_long_channels,5)
# fig.set_figwidth(16)
# fig.set_figheight(8)
# for ichannel in range(n_long_channels):
      
#     lims = dict(hbo=[-5e-2, 20e-2], hbr=[-5e-2, 20e-2])
time = np.linspace(tmin,tmax,num=n_timepoints)
#     baseline_start_index = 0
#     baseline_end_index = int(5*fs)
    
#     # Plot ITD50
#     # HbO
#     curr_ax = axes[ichannel, 0] 
#     curr_data = subject_data_itd50_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r')
#     # HbR
#     curr_data = subject_data_itd50_hbr_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='b')
#     curr_ax.set_ylim((ymin, ymax))
#     curr_ax.set_xticks(np.linspace(-5,20,num=6))
#     if ichannel == 0:
#         curr_ax.set_title('50 us ITD', fontsize = 24)
#     elif ichannel == 3:
#         curr_ax.set_ylabel(r'$\Delta$Hb ($\mu$M)', usetex=False, fontsize = 24)
#     if ichannel == 5:
#         curr_ax.set_xlabel('Time (s)', fontsize = 24)
#         curr_ax.set_xticklabels(np.linspace(-5,20,num=6))
#     else:
#         curr_ax.set_xticklabels(["","","","","",""])
        
#     # Plot ITD500
#     curr_ax = axes[ichannel,1]
#     # HbO
#     curr_data = subject_data_itd500_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r')
#     # HbR
#     curr_data = subject_data_itd500_hbr_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='b')
#     curr_ax.set_ylim((ymin, ymax))
#     curr_ax.set_xticks(np.linspace(-5,20,num=6))
#     if ichannel == 0:
#         curr_ax.set_title('500 us ITD', fontsize = 24)
#     if ichannel == 5:
#         curr_ax.set_xlabel('Time (s)', fontsize = 24)
#         curr_ax.set_xticklabels(np.linspace(-5,20,num=6))
#     else:
#         curr_ax.set_xticklabels(["","","","","",""])
#     # Plot ild5deg
#     curr_ax = axes[ichannel,2]
#     # HbO
#     curr_data = subject_data_ild5deg_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r')
#     # HbR
#     curr_data = subject_data_ild5deg_hbr_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='b')
#     curr_ax.set_ylim((ymin, ymax))
#     curr_ax.set_xticks(np.linspace(-5,20,num=6))
#     if ichannel == 0:
#         curr_ax.set_title('ILD5deg', fontsize = 24)
#     if ichannel == 5:
#         curr_ax.set_xlabel('Time (s)', fontsize = 24)
#         curr_ax.set_xticklabels(np.linspace(-5,20,num=6))
#     else:
#         curr_ax.set_xticklabels(["","","","","",""])
#     # Plot ild5degMag
#     curr_ax = axes[ichannel,3]
#     # HbO
#     curr_data = subject_data_ild5degMag_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r')
#     curr_ax.set_xticks(np.linspace(-5,20,num=6))
#     if ichannel == 5:
#         curr_ax.set_xlabel('Time (s)', fontsize = 24)
#         curr_ax.set_xticklabels(np.linspace(-5,20,num=6))
#     else:
#         curr_ax.set_xticklabels(["","","","","",""])
#     #HbR
#     curr_data = subject_data_ild5degMag_hbr_baselined[:,ichannel,:]
#     curr_mean = np.nanmean(curr_data, axis=0) 
#     curr_error = np.nanstd(curr_data, axis=0)/np.sqrt(np.size(curr_data, axis=0) - 1)
#     curr_ax.plot(time, curr_mean, 'k-')
#     curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='b')
#     curr_ax.set_ylim((ymin, ymax))
#     if ichannel == 0:
#         curr_ax.set_title('ILD5deg + Mag', fontsize = 24)
    
    
#     # Plot sensor location
#     curr_ax = axes[ichannel,4]
#     epochs.copy().pick(chan_hbo[ichannel]).plot_sensors(axes = curr_ax)
# plt.subplots_adjust(top=.9, right=.98, left= 0.05, bottom = 0.07)
# plt.savefig('C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\fNIRS_Plots\\PFC_Block_Averages.svg', format='svg')




# ---------------------------------------------------------------
# -----------------     Plot Grand Average PFC          ---------
# ---------------------------------------------------------------

# ymin = -5e-2
# ymax = 20e-2
# fig, axes = plt.subplots(1, 2)
# fig.set_figwidth(16)
# fig.set_figheight(8)
# # ITD 50 vs ITD 500
# curr_ax = axes[0]
# curr_data = np.nanmean(subject_data_itd50_baselined[:,:,:], axis = 1) # mean over channel
# curr_mean = np.nanmean(curr_data, axis = 0) # mean over subject
# curr_error = np.nanstd(curr_data, axis = 0)/np.sqrt(np.size(curr_data, axis=0) - 1)
# line1 = curr_ax.plot(time, curr_mean, 'ro', markerfacecolor= 'white', label = '50 us ITD')
# curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r', alpha = 0.1)

# curr_data = np.nanmean(subject_data_itd500_baselined[:,:,:], axis = 1) # mean over channel
# curr_mean = np.nanmean(curr_data, axis = 0) # mean over subject
# curr_error = np.nanstd(curr_data, axis = 0)/np.sqrt(np.size(curr_data, axis=0) - 1)
# line2 = curr_ax.plot(time, curr_mean, 'ro', markerfacecolor= 'r', label = '500 us ITD')
# curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r', alpha = 0.1)
# curr_ax.legend()
# curr_ax.set_xlabel('Time (s)',fontsize=24)
# curr_ax.set_ylabel(r'$\Delta$HbO ($\mu$M)', usetex=False, fontsize = 24)

# # ILD 
# curr_ax = axes[1]
# curr_data = np.nanmean(subject_data_ild5deg_baselined[:,:,:], axis = 1) # mean over channel
# curr_mean = np.nanmean(curr_data, axis = 0) # mean over subject
# curr_error = np.nanstd(curr_data, axis = 0)/np.sqrt(np.size(curr_data, axis=0) - 1)
# line1 = curr_ax.plot(time, curr_mean, 'ro', markerfacecolor= 'white', label = '5 deg ILD')
# curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r', alpha = 0.1)

# curr_data = np.nanmean(subject_data_ild5degMag_baselined[:,:,:], axis = 1) # mean over channel
# curr_mean = np.nanmean(curr_data, axis = 0) # mean over subject
# curr_error = np.nanstd(curr_data, axis = 0)/np.sqrt(np.size(curr_data, axis=0) - 1)
# line2 = curr_ax.plot(time, curr_mean, 'ro', markerfacecolor= 'r', label = '5 deg ILD + Mag')
# curr_ax.fill_between(time, curr_mean- curr_error,  curr_mean+curr_error, color='r', alpha = 0.1)
# curr_ax.legend()
# curr_ax.set_xlabel('Time (s)',fontsize=24)
# plt.savefig('C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\fNIRS_Plots\\PFC_Grand_Average.svg', format='svg')


# ---------------------------------------------------------------
# -----------------     PLotting Block Average Means PFC---------
# ---------------------------------------------------------------
# fig, axes = plt.subplots(n_long_channels, 2)
# fig.set_figwidth(10)
# fig.set_figheight(9)
# caxis_lim = 11e-8


# # Plot error bars
# for ichannel in range(n_long_channels):
#     curr_axes = axes[ichannel, 0]
#     # ITD50
#     curr_axes.errorbar(0.5, np.nanmean(mean_during_stim_itd50[:,ichannel], axis=0), 
#                        np.nanstd(mean_during_stim_itd50[:,ichannel],axis=0)/(np.sqrt(np.size(mean_during_stim_itd50, axis=0) - 1)),
#                        fmt='o', capsize= 5.0)    
#     # ITD500
#     curr_axes.errorbar(1, np.nanmean(mean_during_stim_itd500[:,ichannel], axis=0), 
#                        np.nanstd(mean_during_stim_itd500[:,ichannel],axis=0)/(np.sqrt(np.size(mean_during_stim_itd500, axis=0) - 1) ),
#                        fmt='o', capsize= 5.0)
#     # ild5deg
#     curr_axes.errorbar(1.5, np.nanmean(mean_during_stim_ild5deg[:,ichannel], axis=0),
#                        np.nanstd(mean_during_stim_ild5deg[:,ichannel],axis=0)/(np.sqrt(np.size(mean_during_stim_ild5deg, axis=0)  - 1)),
#                        fmt='o', capsize= 5.0)
#     # ild5degMag 
#     curr_axes.errorbar(2, np.nanmean(mean_during_stim_ild5degMag[:,ichannel], axis=0), 
#                        np.nanstd(mean_during_stim_ild5degMag[:,ichannel],axis=0)/(np.sqrt(np.size(mean_during_stim_ild5degMag, axis=0) - 1) ),
#                        fmt='o', capsize= 5.0)
    
#     curr_axes.set_ylim((0,0.15))
#     if ichannel == 2:
#         curr_axes.set_ylabel(r'Mean $\Delta$Hb ($\mu$M) during stim.', usetex=False, fontsize=24)
#     #curr_axes.set_title(channel_names[ichannel])
#     curr_axes.set_xticks([0.5,1,1.5,2])
#     curr_axes.set_xlim([0.4,2.1])
#     if ichannel == 5:
#         curr_axes.set_xticklabels(["Small\n ITD","Large\n ITD","Natural\n ILD","Broadband\n ILD"], fontsize = 18)
#     else:
#         curr_axes.set_xticklabels(["","","",""])
        
#     curr_axes = axes[ichannel, 1]
#     epochs.copy().pick(chan_hbo[ichannel]).plot_sensors(axes = curr_axes)
    
# fig.tight_layout()
# plt.subplots_adjust(top=.98, right=.999, left= 0.1, bottom = 0.1)
# plt.savefig(f'C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\fNIRS_Plots\\Mean_HbO_PFC.svg', format='svg')



# ---------------------------------------------------------------
# -----------------     PLotting GLM Averages           ---------
# ---------------------------------------------------------------
# caxis_lim = 8e-8

# # caxis_lim = np.nanmean([np.nanmean(np.abs(subject_data_itd50_GLM_mean)),
# #             np.nanmean(np.abs(subject_data_itd500_GLM_mean)),
# #             np.nanmean(np.abs(subject_data_ild5deg_GLM_mean)),
# # #             np.nanmean(np.abs(subject_data_ild5degMag_GLM_mean))])
# fig, axes = plt.subplots(int(n_long_channels/4),4)
# fig.set_figwidth(4)
# fig.set_figheight(8)


# # Plot error bars
# for ichannel, curr_axes in enumerate(axes.reshape(-1)):
#     # ITD50
#     curr_axes.errorbar(1, subject_data_itd50_GLM_mean[ichannel],np.nanstd(subject_data_itd50_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_itd50_GLM, axis=0) - 1)),fmt='o')    
#     # ITD500
#     curr_axes.errorbar(2, subject_data_itd500_GLM_mean[ichannel], 
#                         np.nanstd(subject_data_itd500_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_itd500_GLM, axis=0)- 1) ),
#                         fmt='o')
#     # ild5deg
#     curr_axes.errorbar(3, subject_data_ild5deg_GLM_mean[ichannel],
#                         np.nanstd(subject_data_ild5deg_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_ild5deg_GLM, axis=0)- 1) ),
#                         fmt='o')
#     # ild5degMag 
#     curr_axes.errorbar(4, subject_data_ild5degMag_GLM_mean[ichannel], 
#                         np.nanstd(subject_data_ild5degMag_GLM[:,ichannel],axis=0)/(np.sqrt(np.size(subject_data_ild5degMag_GLM, axis=0)- 1) ),
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
caxis_lim = 0.06

fig, (ax1,ax2,ax3,ax4) = plt.subplots(1, 4)
im, _ = mne.viz.plot_topomap(subject_data_itd50_GLM_mean, epochs.pick('hbo').info,
                      extrapolate='box',  image_interp='nearest',
                              vlim=(-caxis_lim, caxis_lim), cmap ='RdBu_r', axes=ax1, show=False)
#cbar = fig.colorbar(im, ax=ax1)
#cbar.set_label('Beta (a.u.)')
ax1.set_title('GLM Beta: ITD50')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_itd500_GLM_mean, epochs.pick('hbo').info,
                      extrapolate='box',  image_interp='nearest',
                              vlim=(-caxis_lim, caxis_lim), cmap ='RdBu_r', axes=ax2, show=False)
#cbar = fig.colorbar(im, ax=ax2)
#cbar.set_label('Beta (a.u.)')
ax2.set_title('GLM Beta: ITD500')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_ild5deg_GLM_mean, epochs.pick('hbo').info,
                      extrapolate='box', image_interp='nearest',
                              vlim=(-caxis_lim, caxis_lim), cmap ='RdBu_r', axes=ax3, show=False)
#cbar = fig.colorbar(im, ax=ax3)
#cbar.set_label('Beta (a.u.)')
ax3.set_title('GLM Beta: ild5deg')
plt.show()

im, _ = mne.viz.plot_topomap(subject_data_ild5degMag_GLM_mean, epochs.pick('hbo').info,
                      extrapolate='box', image_interp='nearest',
                              vlim=(-caxis_lim, caxis_lim), cmap ='RdBu_r', axes=ax4, show=False)
cbar = fig.colorbar(im, ax=ax4)
cbar.set_label('Beta (a.u.)')
ax4.set_title('GLM Beta: ild5degMag')
plt.show()



# ---------------------------------------------------------------
# -----------------     Compare Mean HbO and GLM Values ---------
# ---------------------------------------------------------------
# fig, axes = plt.subplots(2, 7)
# fig.set_figwidth(16)
# fig.set_figheight(8)
# caxis_lim = 11e-8
# for ichannel, curr_axes in enumerate(axes.reshape(-1)):
#     # ITD50
#     curr_axes.scatter(mean_during_stim_itd50[:,ichannel],subject_data_itd50_GLM[:,ichannel])
    
#     # ITD500
#     curr_axes.scatter(mean_during_stim_itd500[:,ichannel],subject_data_itd500_GLM[:,ichannel])
    
#     # ild5deg
#     curr_axes.scatter(mean_during_stim_ild5deg[:,ichannel],subject_data_ild5deg_GLM[:,ichannel])
    
#     # ild5degMag
#     curr_axes.scatter(mean_during_stim_ild5degMag[:,ichannel],subject_data_ild5degMag_GLM[:,ichannel])
    
#     curr_axes.set_xlabel(r'Mean $\Delta$Hb ($\mu$M) during stim.', usetex=False)
#     curr_axes.set_ylabel('Beta')
# plt.savefig(f'C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\Plots\\Beta_vs_HbO_dur_{glm_dur}_{masker_type}_masker.png')



# # Uncorrected GLM
# names = ['S','Channel']
# index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd50_GLM.shape], names = names)
# itd50_df_GLM = pd.DataFrame({'subject_data_itd50_GLM':subject_data_itd50_GLM.flatten()},index=index)['subject_data_itd50_GLM']
# itd500_df_GLM = pd.DataFrame({'subject_data_itd500_GLM':subject_data_itd500_GLM.flatten()},index=index)['subject_data_itd500_GLM']
# ild5deg_df_GLM = pd.DataFrame({'subject_data_ild5deg_GLM':subject_data_ild5deg_GLM.flatten()},index=index)['subject_data_ild5deg_GLM']
# ild5degMag_df_GLM = pd.DataFrame({'subject_data_ild5degMag_GLM':subject_data_ild5degMag_GLM.flatten()},index=index)['subject_data_ild5degMag_GLM']
# z = pd.concat([itd50_df_GLM,itd500_df_GLM,ild5deg_df_GLM,ild5degMag_df_GLM], ignore_index=True,axis=1)
# z.to_csv(f'all_subjects_uncorr_GLM_{masker_type}_masker.csv',index=True)


# # Corrected GLM

# names = ['S','Channel']
# index = pd.MultiIndex.from_product([range(s) for s in subject_data_itd50_GLM.shape], names = names)
# itd50_df_GLM_bh_corr = pd.DataFrame({'subject_data_itd50_GLM_bh_corr':subject_data_itd50_GLM_bh_corr.flatten()},index=index)['subject_data_itd50_GLM_bh_corr']
# itd500_df_GLM_bh_corr = pd.DataFrame({'subject_data_itd500_GLM_bh_corr':subject_data_itd500_GLM_bh_corr.flatten()},index=index)['subject_data_itd500_GLM_bh_corr']
# ild5deg_df_GLM_bh_corr = pd.DataFrame({'subject_data_ild5deg_GLM_bh_corr':subject_data_ild5deg_GLM_bh_corr.flatten()},index=index)['subject_data_ild5deg_GLM_bh_corr']
# ild5degMag_df_GLM_bh_corr = pd.DataFrame({'subject_data_ild5degMag_GLM_bh_corr':subject_data_ild5degMag_GLM_bh_corr.flatten()},index=index)['subject_data_ild5degMag_GLM_bh_corr']
# z = pd.concat([itd50_df_GLM_bh_corr,itd500_df_GLM_bh_corr,ild5deg_df_GLM_bh_corr,ild5degMag_df_GLM_bh_corr], ignore_index=True,axis=1)
# z.to_csv(f'all_subjects_bh_corr_GLM_{masker_type}_masker.csv',index=True)