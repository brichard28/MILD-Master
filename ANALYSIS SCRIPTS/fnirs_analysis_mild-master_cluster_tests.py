import os
import numpy as np
import mne
import pickle
from scipy.io import loadmat
from concurrent.futures import ProcessPoolExecutor
from mne_cluster_based_permutation_test_cpu import mne_cluster_based_permutation_test_cpu
from mne_compute_channel_positions import mne_compute_channel_positions
from mne_cluster_based_adjacency_nirs import mne_cluster_based_adjacency_nirs
import time
from datetime import date
from mne_modified_beer_lambert_law import mne_modified_beer_lambert_law
import pandas as pd
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)

user = 'Home'
if user == 'Laptop':
    data_root = 'C:/Users/benri/Downloads/'
    mild_master_root = 'C:/Users/benri/Documents/GitHub/MILD-Master'

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
data_root + "2025-02-24/2025-02-24_002",
data_root + "2025-02-26/2025-02-26_002",
data_root + "2025-02-27/2025-02-27_001",
data_root + "2025-02-27/2025-02-27_002",
data_root + "2025-03-03/2025-03-03_001",
data_root + "2025-03-05/2025-03-05_001",
data_root + "2025-03-10/2025-03-10_001",
data_root + "2025-03-12/2025-03-12_001",
data_root + "2025-03-12/2025-03-12_002",
data_root + "2025-03-17/2025-03-17_001",
data_root + "2025-03-18/2025-03-18_001",
data_root + "2025-04-23/2025-04-23_001",
data_root + "2025-04-24/2025-04-24_001",
data_root + "2025-05-08/2025-05-08_001"]

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
'mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40',
'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_45',
'mild_master_46','mild_master_47','mild_master_48']


# The subjects we would like to run right now
curr_subject_ID = ['mild_master_1',
'mild_master_3',
'mild_master_4',
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
'mild_master_30',
'mild_master_31',
'mild_master_32',
'mild_master_33',
'mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40',
'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_46','mild_master_47','mild_master_48']
curr_folder_indices = [index for index, element in enumerate(subject_ID) if np.isin(element,curr_subject_ID)]
curr_fnirs_data_folders = [all_fnirs_data_folders[i] for i in curr_folder_indices]

conditions = ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']

save_dir = mild_master_root + "/RESULTS DATA/"

# Load the data (group_results)

group_results = pd.read_csv(mild_master_root + "/RESULTS DATA/group_results_glm_dur_11.csv")
beta_data = group_results.query("Chroma in ['hbo']").groupby(by=['ch_name','ID','Condition'],as_index=False)['theta'].mean()
# divide data by condition
beta_small_itd_df = beta_data.query("Condition in ['az_itd=5_az=0']")
beta_large_itd_df = beta_data.query("Condition in ['az_itd=15_az=0']")
beta_small_ild_df = beta_data.query("Condition in ['az_itd=0_az=5']")
beta_large_ild_df = beta_data.query("Condition in ['az_itd=0_az=15']")


# Convert to subjects x channels arrays
beta_small_itd = beta_small_itd_df.copy().pivot(index='ch_name',columns='ID',values='theta').to_numpy().T
beta_large_itd = beta_large_itd_df.copy().pivot(index='ch_name',columns='ID',values='theta').to_numpy().T
beta_small_ild = beta_small_ild_df.copy().pivot(index='ch_name',columns='ID',values='theta').to_numpy().T
beta_large_ild = beta_large_ild_df.copy().pivot(index='ch_name',columns='ID',values='theta').to_numpy().T


# Load raw haemo
data = mne.io.read_raw_nirx(f"{curr_fnirs_data_folders[0]}/{curr_fnirs_data_folders[0][-14:]}_config.hdr",
                                    verbose=False, preload=True)
data_snirf = mne.io.read_raw_snirf(f"{curr_fnirs_data_folders[0]}/{curr_fnirs_data_folders[0][-14:]}.snirf",
                                    optode_frame="mri", preload=True)
events, event_dict = mne.events_from_annotations(data, verbose=False)
raw_haemo_for_plotting, null, null = preprocess_NIRX(data, data_snirf, event_dict,
                                           save=False,
                                           savename='non.fif',
                                           plot_steps=False,
                                           crop=False, crop_low=0, crop_high=0,
                                           events_modification=False, reject=True,
                                           short_regression=False, events_from_snirf=False,
                                           drop_short=False, negative_enhancement=False,
                                           snr_thres=3, sci_thres=0.8, filter_type='iir', filter_limits=[0.01,0.1], filter_transition_bandwidths=[0.005, 0.1/2])
raw_haemo_for_plotting.info['bads'] = []
info = get_long_channels(raw_haemo_for_plotting).pick(picks='hbo').info


# Compute adjacency
n_channels = len(info.ch_names)
pos, outlines = mne_compute_channel_positions(info, n_channels)
adjacency_sparse = mne_cluster_based_adjacency_nirs(info=info, pos=pos, threshold_mm=20)

alpha = 0.05

ch_sorted = np.unique(beta_small_itd_df['ch_name'])
ch_ind = [np.where(ch_n == np.array(info.ch_names))[0] for i, ch_n in enumerate(ch_sorted)]

# Put data in format for X
beta_small_itd_reind = np.squeeze(beta_small_itd[:, ch_ind])
beta_large_itd_reind = np.squeeze(beta_large_itd[:, ch_ind])
beta_small_ild_reind = np.squeeze(beta_small_ild[:, ch_ind])
beta_large_ild_reind = np.squeeze(beta_large_ild[:, ch_ind])

num_channels = np.shape(beta_large_ild_reind)[1]
# Generate the null Gaussian distribution
mean = 0  # Center the distribution at 0
std_dev = np.nanstd(beta_data['theta'])  # Standard deviation (adjust for spread)
num_samples = np.shape(beta_large_ild_reind)[0]  # Number of data points to generate

null_data = np.empty(np.shape(beta_small_itd_reind))
for ichannel in range(num_channels):
    gaussian_data = np.random.normal(loc=mean, scale=std_dev, size=num_samples)
    null_data[:, ichannel] = gaussian_data

# Run cluster based permutation tests

T_obs, clusters, cluster_p_values, H0, significant_channels = (
    mne_cluster_based_permutation_test_cpu(adjacency_sparse=adjacency_sparse,
                                           X=[beta_large_ild_reind, null_data]))



