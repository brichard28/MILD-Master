
import numpy as np
import mne
import math
import matplotlib
from matplotlib import pyplot as plt
import os
import pandas as pd
from collections import defaultdict

mne.set_config('MNE_BROWSER_BACKEND', 'qt')
#from nirx_movement import mark_aux_movement_bad
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm, statsmodels_to_results
from mne_nirs.channels import (get_long_channels,
                               get_short_channels)
from nilearn.plotting import plot_design_matrix
from run_preproc_NIRS import preprocess_NIRX
from mne_nirs.io.snirf import read_snirf_aux_data
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from mne_nirs.channels import picks_pair_to_idx
from mne_nirs.visualisation import plot_glm_group_topo, plot_glm_surface_projection
from plot_nirs import plot_nirs_evoked_error

from scipy import signal
# ---------------------------------------------------------------
# -----------------          Data Parameters            ---------
# ---------------------------------------------------------------
wdir = os.path.dirname(__file__)

# Define Subject Files
# Define Subject Files
root = ''
user = 'Desktop'
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
data_root + "2025-02-24/2025-02-24_002",
data_root + "2025-02-26/2025-02-26_002",
data_root + "2025-02-27/2025-02-27_001",
data_root + "2025-02-27/2025-02-27_002",
data_root + "2025-03-03/2025-03-03_001",
data_root + "2025-03-05/2025-03-05_001"]

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
'mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40']


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
'mild_master_34'] #,'mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40'

curr_folder_indices = [index for index, element in enumerate(subject_ID) if np.isin(element,curr_subject_ID)]
curr_fnirs_data_folders = [all_fnirs_data_folders[i] for i in curr_folder_indices]

glm_dur = 4

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


right_hem_channels = [[4,4],[4,3],[4,2],[5,3],[5,2],[5,1],[6,2],[6,1],[11,13],[11,4],[11,3],[11,12],[11,11],[12,4],[12,3],
    [12,2],[12,12],[12,11],[12,10],[13,11],[13,10],[13,9],[13,3],[13,2],[13,1],[14,2],[14,1],[14,10],[14,9],[20,13],[20,12],
    [20,20],[21,20],[21,12],[21,11],[21,19],[22,20],[22,12],[22,11],[22,10],[22,19],[23,11],[23,10],[23,9],[23,19],[24,10],
    [24,9],[24,19]]
right_hem_channel_names = ["S" + str(value[0]) + "_D" + str(value[1]) + " hbo" for idx, value in enumerate(right_hem_channels)]
group_df = pd.DataFrame()  # To store channel level results


# load in an example raw_haemo for this montage
data = mne.io.read_raw_nirx(f"{curr_fnirs_data_folders[0]}/{curr_fnirs_data_folders[0][-14:]}_config.hdr",
                                    verbose=False, preload=True)
data_snirf = mne.io.read_raw_snirf(f"{curr_fnirs_data_folders[0]}/{curr_fnirs_data_folders[0][-14:]}.snirf",
                                    optode_frame="mri", preload=True)
events, event_dict = mne.events_from_annotations(data, verbose=False)
raw_haemo_for_plotting, null = preprocess_NIRX(data, data_snirf, event_dict,
                                           save=False,
                                           savename='non.fif',
                                           plot_steps=False,
                                           crop=False, crop_low=0, crop_high=0,
                                           events_modification=False, reject=True,
                                           short_regression=False, events_from_snirf=False,
                                           drop_short=False, negative_enhancement=False,
                                           snr_thres=3, sci_thres=0.8, filter_type='iir', filter_limits=[0.01,0.5])
raw_haemo_for_plotting.info['bads'] = []






for ii, subject_num in enumerate(range(n_subjects)):
    subject = curr_subject_ID[ii]
    individual_results = pd.read_csv(mild_master_root + "/RESULTS DATA/" + subject + "_results.csv")
    group_df = pd.concat([group_df, individual_results], ignore_index=True)

group_results = group_df.query("Condition in ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']")
group_results.to_csv(mild_master_root + "/RESULTS DATA/group_results.csv")


caxis_min = -0.1
caxis_max = 0.1
groups_single_chroma = dict(
    Left_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), left_hem_channels,
                                      on_missing='warning'),
    Right_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), right_hem_channels,
                                       on_missing='warning'))



group_theta_for_topoplot = group_results.query("Chroma in ['hbo']").groupby(by=['ch_name','Condition'],as_index=False)['theta'].mean()
fig, topo_axes = plt.subplots(nrows=1, ncols=4,figsize=(18,10))

this_info_left = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_left.drop_channels([val for idx, val in enumerate(this_info_left.ch_names) if val not in left_hem_channel_names])
this_info_left.drop_channels([i for i in this_info_left.ch_names if i not in np.unique(group_theta_for_topoplot['ch_name'])])
this_info_left = this_info_left.info

this_info_right = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_right.drop_channels([val for idx, val in enumerate(this_info_right.ch_names) if val not in right_hem_channel_names])
this_info_right.drop_channels([i for i in this_info_right.ch_names if i not in np.unique(group_theta_for_topoplot['ch_name'])])
this_info_right = this_info_right.info

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[0],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[0],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[1],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[1],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[2],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[2],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[3],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[3],
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_beta.png")
plt.close(fig)







# ---------------------------------------------------------------
# -----------------     Topomap of Mean HbO             ---------
#----------------------------------------------------------------
caxis_min = -0.1
caxis_max = 0.1
group_mean_hbo_for_topoplot = group_results.query("Chroma in ['hbo']").groupby(by=['ch_name','Condition'],as_index=False)['mean_hbo'].mean()
group_mean_hbo_for_topoplot.loc[np.isnan(group_mean_hbo_for_topoplot['mean_hbo']),"mean_hbo"] = 0

fig, topo_axes = plt.subplots(nrows=1, ncols=4,figsize=(18,10))

this_info_left = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_left.drop_channels([val for idx, val in enumerate(this_info_left.ch_names) if val not in left_hem_channel_names])
this_info_left.drop_channels([i for i in this_info_left.ch_names if i not in np.unique(group_mean_hbo_for_topoplot['ch_name'])])
this_info_left = this_info_left.info

this_info_right = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_right.drop_channels([val for idx, val in enumerate(this_info_right.ch_names) if val not in right_hem_channel_names])
this_info_right.drop_channels([i for i in this_info_right.ch_names if i not in np.unique(group_mean_hbo_for_topoplot['ch_name'])])
this_info_right = this_info_right.info

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[0],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[0],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[1],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[1],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[2],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[2],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[3],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[3],
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_mean_hbo.png")
plt.close(fig)