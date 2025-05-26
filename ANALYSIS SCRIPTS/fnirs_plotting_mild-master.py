
import numpy as np
import mne
import math
import matplotlib
from matplotlib import pyplot as plt
import os
import pandas as pd
from collections import defaultdict


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

import statsmodels.formula.api as smf

from mne.preprocessing.nirs import source_detector_distances, _channel_frequencies, _check_channels_ordered
from mne.channels.layout import find_layout
from copy import deepcopy
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
    mild_master_root = 'C:/Users/benri/Documents/GitHub/MILD-Master'
    mne.set_config('MNE_BROWSER_BACKEND', 'qt')
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






for ii, subject_num in enumerate(range(n_subjects)):
    subject = curr_subject_ID[ii]
    individual_results = pd.read_csv(mild_master_root + "/RESULTS DATA/" + subject + "_results.csv")
    group_df = pd.concat([group_df, individual_results], ignore_index=True)

    ## TOPOMAP OF INDIVIDUAL RESULTS
    caxis_min = -0.3
    caxis_max = 0.3
    groups_single_chroma = dict(
        Left_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), left_hem_channels,
                                          on_missing='warning'),
        Right_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), right_hem_channels,
                                           on_missing='warning'))

    individual_mean_hbo_for_topoplot = individual_results.query("Chroma in ['hbo']").groupby(by=['ch_name', 'Condition'], as_index=False)[
            'mean_hbo'].mean()
    fig, topo_axes = plt.subplots(nrows=1, ncols=4, figsize=(18, 10))

    this_info_left = raw_haemo_for_plotting.copy().pick(picks="hbo")
    this_info_left.drop_channels(
        [val for idx, val in enumerate(this_info_left.ch_names) if val not in left_hem_channel_names])
    this_info_left.drop_channels(
        [i for i in this_info_left.ch_names if i not in np.unique(individual_mean_hbo_for_topoplot['ch_name'])])
    this_info_left = this_info_left.info

    this_info_right = raw_haemo_for_plotting.copy().pick(picks="hbo")
    this_info_right.drop_channels(
        [val for idx, val in enumerate(this_info_right.ch_names) if val not in right_hem_channel_names])
    this_info_right.drop_channels(
        [i for i in this_info_right.ch_names if i not in np.unique(individual_mean_hbo_for_topoplot['ch_name'])])
    this_info_right = this_info_right.info

    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=5_az=0']").query(
        "ch_name in @this_info_left['ch_names']")['mean_hbo'],
                         this_info_left, sensors=True, axes=topo_axes[0], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=5_az=0']").query(
        "ch_name in @this_info_right['ch_names']")['mean_hbo'],
                         this_info_right, sensors=True, axes=topo_axes[0], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))

    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query(
        "ch_name in @this_info_left['ch_names']")['mean_hbo'],
                         this_info_left, sensors=True, axes=topo_axes[1], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query(
        "ch_name in @this_info_right['ch_names']")['mean_hbo'],
                         this_info_right, sensors=True, axes=topo_axes[1], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))

    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query(
        "ch_name in @this_info_left['ch_names']")['mean_hbo'],
                         this_info_left, sensors=True, axes=topo_axes[2], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query(
        "ch_name in @this_info_right['ch_names']")['mean_hbo'],
                         this_info_right, sensors=True, axes=topo_axes[2], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query(
        "ch_name in @this_info_left['ch_names']")['mean_hbo'],
                         this_info_left, sensors=True, axes=topo_axes[3], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query(
        "ch_name in @this_info_right['ch_names']")['mean_hbo'],
                         this_info_right, sensors=True, axes=topo_axes[3], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    plt.savefig(mild_master_root + "/CASUAL FIGURES/" + f"{subject}_mean_hbo_topoplot.png")
    plt.close(fig)

    individual_theta_for_topoplot = \
    individual_results.query("Chroma in ['hbo']").groupby(by=['ch_name', 'Condition'], as_index=False)[
        'theta'].mean()
    fig, topo_axes = plt.subplots(nrows=1, ncols=4, figsize=(18, 10))

    this_info_left = raw_haemo_for_plotting.copy().pick(picks="hbo")
    this_info_left.drop_channels(
        [val for idx, val in enumerate(this_info_left.ch_names) if val not in left_hem_channel_names])
    this_info_left.drop_channels(
        [i for i in this_info_left.ch_names if i not in np.unique(individual_theta_for_topoplot['ch_name'])])
    this_info_left = this_info_left.info

    this_info_right = raw_haemo_for_plotting.copy().pick(picks="hbo")
    this_info_right.drop_channels(
        [val for idx, val in enumerate(this_info_right.ch_names) if val not in right_hem_channel_names])
    this_info_right.drop_channels(
        [i for i in this_info_right.ch_names if i not in np.unique(individual_theta_for_topoplot['ch_name'])])
    this_info_right = this_info_right.info

    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=5_az=0']").query(
        "ch_name in @this_info_left['ch_names']")['theta'],
                         this_info_left, sensors=True, axes=topo_axes[0], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=5_az=0']").query(
        "ch_name in @this_info_right['ch_names']")['theta'],
                         this_info_right, sensors=True, axes=topo_axes[0], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))

    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query(
        "ch_name in @this_info_left['ch_names']")['theta'],
                         this_info_left, sensors=True, axes=topo_axes[1], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query(
        "ch_name in @this_info_right['ch_names']")['theta'],
                         this_info_right, sensors=True, axes=topo_axes[1], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))

    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query(
        "ch_name in @this_info_left['ch_names']")['theta'],
                         this_info_left, sensors=True, axes=topo_axes[2], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query(
        "ch_name in @this_info_right['ch_names']")['theta'],
                         this_info_right, sensors=True, axes=topo_axes[2], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query(
        "ch_name in @this_info_left['ch_names']")['theta'],
                         this_info_left, sensors=True, axes=topo_axes[3], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    mne.viz.plot_topomap(individual_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query(
        "ch_name in @this_info_right['ch_names']")['theta'],
                         this_info_right, sensors=True, axes=topo_axes[3], contours=0,
                         extrapolate='local', image_interp='linear', vlim=(caxis_min, caxis_max))
    plt.savefig(mild_master_root + "/CASUAL FIGURES/" + f"{subject}_theta_topoplot.png")
    plt.close(fig)

group_results = group_df.query("Condition in ['az_itd=5_az=0','az_itd=15_az=0','az_itd=0_az=5','az_itd=0_az=15']")
# group_results.to_csv(mild_master_root + "/RESULTS DATA/group_results.csv")

ch_model = smf.mixedlm("theta ~ ch_name:Chroma:Condition",group_results,groups=group_results["ID"]).fit(method="nm")
ch_model_df = statsmodels_to_results(ch_model)

caxis_min = -0.05
caxis_max = 0.05
groups_single_chroma = dict(
    Left_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), left_hem_channels,
                                      on_missing='warning'),
    Right_Hemisphere=picks_pair_to_idx(raw_haemo_for_plotting.copy().pick(picks='hbo'), right_hem_channels,
                                       on_missing='warning'))



# ---------------------------------------------------------------
# -----------------     Topomap of Mean Beta            ---------
#----------------------------------------------------------------

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
                     this_info_left,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['theta'],
                     this_info_left,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_theta_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['theta'],
                     this_info_right,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_beta_glm_dur_11.png")
plt.close(fig)

# ---------------------------------------------------------------
# -----------------     Topomap of beta Value Histograms   ---------
#----------------------------------------------------------------

layout = find_layout(get_long_channels(raw_haemo_for_plotting).info)
layout = deepcopy(layout)
layout.pos[:, :2] -= layout.pos[:, :2].min(0)
layout.pos[:, :2] /= layout.pos[:, :2].max(0)
positions = layout.pos[:, :2] * 0.9
# set up subplots
fig = plt.figure(figsize=(5, 4), dpi=200)

width, height = 0.05, 0.05
lims = dict(hbo=[-0.2, 0.2], hbr=[-0.2, -0.2])
# for each channel

unique_positions = np.unique(positions)

unique_markers = np.zeros(np.shape(unique_positions))
for ichannel in range(len(layout.pos)):

    this_channel_name = layout.names[ichannel]
    print(this_channel_name)
    pos = positions[ichannel, :]

    # plot --- [lowerCorner_x, lowerCorner_y, width, height]
    ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])

    this_color = "w"
    if "hbo" in this_channel_name:
        this_pick = "hbo"
        this_color = "r"
    elif "hbr" in this_channel_name:
        this_pick = "hbr"
        this_color = "b"
    plt.hist(group_results.query(f"ch_name in ['{this_channel_name}']").query(f"Chroma in ['{this_pick}']")['theta'], axes=ax, color=this_color, range = [lims['hbo'][0],lims['hbo'][1]])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=0.1, length=2, width=0.5, labelcolor='w')
    ax.set_xlim([lims['hbo'][0],lims['hbo'][1]])
    ax.set_ylim([0,100])
    ax.set_title(f'{layout.names[ichannel][:-4]}', fontsize=3, pad=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_facecolor("none")

# add an empty plot with labels
ax = fig.add_axes([0.5, 0.075, 1.5 * width, 1.5 * height])
ax.set_xlabel('beta value', fontsize=4)
ax.set_ylabel('Freq. of occurrence (count)', fontsize=4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=4)
ax.set_xlim([lims['hbo'][0], lims['hbo'][1]])
ax.set_ylim([0, 100])


plt.savefig(mild_master_root + f"/CASUAL FIGURES/beta_value_histograms_topoplot.png")
plt.close(fig)







# ---------------------------------------------------------------
# -----------------     Topomap of t values             ---------
#----------------------------------------------------------------
caxis_min = -8
caxis_max = 8
group_t_for_topoplot = group_results.query("Chroma in ['hbo']").groupby(by=['ch_name','Condition'],as_index=False)['t'].mean()
fig, topo_axes = plt.subplots(nrows=1, ncols=4,figsize=(18,10))

this_info_left = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_left.drop_channels([val for idx, val in enumerate(this_info_left.ch_names) if val not in left_hem_channel_names])
this_info_left.drop_channels([i for i in this_info_left.ch_names if i not in np.unique(group_t_for_topoplot['ch_name'])])
this_info_left = this_info_left.info

this_info_right = raw_haemo_for_plotting.copy().pick(picks="hbo")
this_info_right.drop_channels([val for idx, val in enumerate(this_info_right.ch_names) if val not in right_hem_channel_names])
this_info_right.drop_channels([i for i in this_info_right.ch_names if i not in np.unique(group_t_for_topoplot['ch_name'])])
this_info_right = this_info_right.info

mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_left['ch_names']")['t'],
                     this_info_left,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['t'],
                     this_info_right,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['t'],
                     this_info_left,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['t'],
                     this_info_right,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['t'],
                     this_info_left,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['t'],
                     this_info_right,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['t'],
                     this_info_left,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_t_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['t'],
                     this_info_right,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_t_values.png")
plt.close(fig)




# ---------------------------------------------------------------
# -----------------     Topomap of Mean t-value Contrasts           ---------
#----------------------------------------------------------------
small_itd_data_left = group_t_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_left['ch_names']")['t'].reset_index(drop=True)
large_itd_data_left = group_t_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['t'].reset_index(drop=True)
small_itd_data_right = group_t_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['t'].reset_index(drop=True)
large_itd_data_right = group_t_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['t'].reset_index(drop=True)

small_ild_data_left = group_t_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['t'].reset_index(drop=True)
large_ild_data_left = group_t_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['t'].reset_index(drop=True)
small_ild_data_right = group_t_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['t'].reset_index(drop=True)
large_ild_data_right = group_t_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['t'].reset_index(drop=True)


caxis_min = -3
caxis_max = 3
fig, topo_contrast_axes = plt.subplots(nrows=3, ncols=2,figsize=(18,10))

# Large ITD - Small ITD
mne.viz.plot_topomap(large_itd_data_left - small_itd_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[0, 0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(large_itd_data_right - small_itd_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[0, 0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[0,0].set_title("Large ITD - Small ITD")

# Large ILD - Small ILD
mne.viz.plot_topomap(large_ild_data_left - small_ild_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[0, 1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(large_ild_data_right - small_ild_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[0, 1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[0, 1].set_title("Large ILD - Small ILD")


# Small ITD - Small ILD
mne.viz.plot_topomap(small_itd_data_left - small_ild_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[1,0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(small_itd_data_right - small_ild_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[1,0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[1,0].set_title("Small ITD - Small ILD")

# Large ITD - Large ILD
mne.viz.plot_topomap(large_itd_data_left - large_ild_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[1,1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(large_itd_data_right - large_ild_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[1,1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[1,1].set_title("Large ITD - Large ILD")


# Small ITD - Large ILD
mne.viz.plot_topomap(small_itd_data_left - large_ild_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[2,0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(small_itd_data_right - large_ild_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[2,0],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[2,0].set_title("Small ITD - Large ILD")

# Large ITD - Small ILD
mne.viz.plot_topomap(large_itd_data_left - small_ild_data_left,
                     this_info_left,sensors=True, axes = topo_contrast_axes[2,1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(large_itd_data_right - small_ild_data_right,
                     this_info_right,sensors=True, axes = topo_contrast_axes[2,1],contours = 0,
                     extrapolate='local',image_interp='linear',vlim=(caxis_min,caxis_max))
topo_contrast_axes[2,1].set_title("Large ITD - Small ILD")

plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_t_contrasts.png")
plt.close(fig)

# ---------------------------------------------------------------
# -----------------     Topomap of t Value Histograms   ---------
#----------------------------------------------------------------

layout = find_layout(get_long_channels(raw_haemo_for_plotting).info)
layout = deepcopy(layout)
layout.pos[:, :2] -= layout.pos[:, :2].min(0)
layout.pos[:, :2] /= layout.pos[:, :2].max(0)
positions = layout.pos[:, :2] * 0.9
# set up subplots
fig = plt.figure(figsize=(5, 4), dpi=200)

width, height = 0.05, 0.05
lims = dict(hbo=[-45, 45], hbr=[-45, 45])
# for each channel

unique_positions = np.unique(positions)

unique_markers = np.zeros(np.shape(unique_positions))
for ichannel in range(len(layout.pos)):

    this_channel_name = layout.names[ichannel]
    print(this_channel_name)
    pos = positions[ichannel, :]

    # plot --- [lowerCorner_x, lowerCorner_y, width, height]
    ax = fig.add_axes([pos[0] + width / 2, pos[1], width, height])

    this_color = "w"
    if "hbo" in this_channel_name:
        this_pick = "hbo"
        this_color = "r"
    elif "hbr" in this_channel_name:
        this_pick = "hbr"
        this_color = "b"
    plt.hist(group_results.query(f"ch_name in ['{this_channel_name}']").query(f"Chroma in ['{this_pick}']")['t'], axes=ax, color=this_color, range = [lims['hbo'][0],lims['hbo'][1]])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=0.1, length=2, width=0.5, labelcolor='w')
    ax.set_title(f'{layout.names[ichannel][:-4]}', fontsize=3, pad=0)
    ax.set_xlim([lims['hbo'][0], lims['hbo'][1]])
    ax.set_ylim([0, 70])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_facecolor("none")

# add an empty plot with labels
ax = fig.add_axes([0.5, 0.075, 1.5 * width, 1.5 * height])
ax.set_xlabel('t value', fontsize=4)
ax.set_ylabel('Freq. of occurrence (count)', fontsize=4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=4)
ax.set_xlim([lims['hbo'][0],lims['hbo'][1]])
ax.set_ylim([0,70])


plt.savefig(mild_master_root + f"/CASUAL FIGURES/t_value_histograms_topoplot.png")
plt.close(fig)





# ---------------------------------------------------------------
# -----------------     Topomap of Mean HbO             ---------
#----------------------------------------------------------------
caxis_min = -0.08
caxis_max = 0.08
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
                     this_info_left,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=5_az=0']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[0],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=15_az=0']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[1],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=5']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[2],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))

mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_left['ch_names']")['mean_hbo'],
                     this_info_left,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
mne.viz.plot_topomap(group_mean_hbo_for_topoplot.query("Condition in ['az_itd=0_az=15']").query("ch_name in @this_info_right['ch_names']")['mean_hbo'],
                     this_info_right,sensors=True, axes = topo_axes[3],contours = 0,
                     extrapolate='local',image_interp = 'linear',vlim=(caxis_min,caxis_max))
plt.savefig(mild_master_root + "/CASUAL FIGURES/group_topoplot_mean_hbo.png")
plt.close(fig)


import mne_nirs
fig, surface_axes = plt.subplots(nrows=1, ncols=1,figsize=(18,10))

group_theta_for_surface=  group_results.groupby(by=['ch_name','Condition'],as_index=False)['theta'].mean()
mne_nirs.visualisation.plot_glm_surface_projection(raw_haemo_for_plotting, group_theta_for_surface.query("Condition in ['az_itd=0_az=15']"), picks='hbo', value="theta")
