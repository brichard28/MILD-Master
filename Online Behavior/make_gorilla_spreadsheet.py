# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:56:07 2024

@author: benri
"""
import os
import pandas as pd
import soundfile as sf
import random
import numpy as np

data = {'display':["instruction"],'cue':[" "],'P1':[" "],'P2':[" "],
                                            'P3':[" "],'P4':[" "],'P5':[" "],'P6':[" "],'P7':[" "],'P8':[" "],'P9':[" "],
                                            'P10':[" "],'P11':[" "],'P12':[" "],'P13':[" "],'P14':[" "],'P15':[" "],'P16':[" "]}

spreadsheet_df = pd.DataFrame(data)

# Set instructions row
# Specify variables, generate conditions
trial_type = "training"
audio_folder = "s_tokentest7"
base_filepath = f'C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\Online Behavior\\stim\\mild-master-online\\{audio_folder}\\'

possible_conditions = os.listdir(base_filepath)

n_trials_per_condition = 8
n_conditions = len(possible_conditions)
n_trials = n_trials_per_condition * n_conditions

conditions = []
for i in range(n_trials_per_condition):
    for cond in possible_conditions:
        conditions.append(cond) 
    
random.shuffle(conditions)

for itrial in range(n_trials):
    spreadsheet_df['display'] = pd.concat([spreadsheet_df['display'], pd.Series(trial_type)], ignore_index=True)

    # Find condition
    this_cond = conditions[itrial]
    
    ## Load in audio files

    # Load in cue audio
    cue_filename = base_filepath + this_cond + "\\cue_" +  this_cond + ".wav"
    spreadsheet_df['cue'].append(pd.Series(cue_filename))
    
    n_tokens = 16
    token_pair_keys = ['P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15','P16']
    
    curr_possible_token_filenames = os.listdir(base_filepath + this_cond)
    curr_possible_token_filenames_bash = [name for name in curr_possible_token_filenames if not 'cue' in name and 'bash' in name]
    curr_possible_token_filenames_nonbash = [name for name in curr_possible_token_filenames if not 'cue' in name and not 'bash' in name]
    
    curr_possible_token_filenames_target_bash = [name for name in curr_possible_token_filenames_bash if 't_bash' in name]
    curr_possible_token_filenames_masker_bash = [name for name in curr_possible_token_filenames_bash if 'm_bash' in name]
    
    
    # Decide where bash goes on this trial
    total_number_bashs = np.random.randint(5,7)
    min_number_bashs = 1
    max_number_bashs = total_number_bashs - min_number_bashs
    n_close_idx = 2
    
    # Next, find indices of bashs for target
    possible_target_bash_indices = np.arange(1,n_tokens)#np.arange(1,len(t_inds)-1) # cannot have bash on first or last one
    valid_bash_indices_target = False
    while not valid_bash_indices_target:
        bash_indices_target = np.sort(np.random.choice(possible_target_bash_indices, size=np.random.randint(min_number_bashs, high=max_number_bashs + 1), replace=False))
        if (np.diff(bash_indices_target) >= n_close_idx).all():
            valid_bash_indices_target = True
    # Generate bash indices for masker that are never within n_close_idx of a target bash
    possible_masker_bash_indices = possible_target_bash_indices.copy()
    possible_masker_bash_indices = np.delete(possible_masker_bash_indices, bash_indices_target - 1)
    valid_bash_indices_masker = False
    while not valid_bash_indices_masker:
        bash_indices_masker = np.sort(np.random.choice(possible_masker_bash_indices, size= total_number_bashs - len(bash_indices_target), replace=False))
        if (np.diff(np.sort(np.append(bash_indices_target,bash_indices_masker))) >= n_close_idx).all() and not np.isin(bash_indices_target,bash_indices_masker).any():
            valid_bash_indices_masker = True
    
    
    
    for itoken in range(n_tokens):
        if itoken in bash_indices_target:
            this_filename = random.choice(curr_possible_token_filenames_target_bash)
        elif itoken in bash_indices_masker:
            this_filename = random.choice(curr_possible_token_filenames_masker_bash)
        else:
            this_filename = random.choice(curr_possible_token_filenames_nonbash)
        
        spreadsheet_df[token_pair_keys[itoken]].append(pd.Series(this_filename))

        
        
        
        
    
    
