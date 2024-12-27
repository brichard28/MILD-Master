# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:21:24 2024

@author: benri

Making target-masker token pair audios for use in mild-master-online
"""
import os, sys
import getopt
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
#from pysndfile.sndio import read as sfread
#from pysndfile.sndio import write as sfwrite
import psylab
import ir
#import spyral
import ild
from soundfile import read
from soundfile import write
import random

def main(subjectID):
    exp = 'mild-master-online'
    
    variables = {
        "side":      ['r',   'r',      'l',      'r',      'l',     'r',      'l',      'r',     'l',      'r',     'l',      'r',     'l',   ],
        "itd":       ['500',   '50',     '50',     '500',    '500',   '0',      '0',      '0',     '0',      '0',     '0',      '0',     '0', ],  
        "az":        ['0',   '0',      '0',      '0',      '0',     '5',      '5',      '10',    '10',     '15',    '15',     '20',     '20', ],
        "mag":       ['0',   '0',      '0',      '0',      '0' ,    '0',      '0',      '0',     '0',      '0',     '0',      '0',     '0',   ],
        		}
    
    
    t_path = 'C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\Online Behavior\\stim\\bashdashgash'
    m_path = 'C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\Online Behavior\\stim\\bashdashgash'
    out_path_stim = f'.\\stim\\{exp}\\s_{subjectID}'
    out_path_data = f'.\\data\\{exp}'
    
    t_name = 'ben'
    m_name = 'ben'
    
    t_s = ['b', 'd', 'g']
    m_s = ['b', 'd', 'g']
     
    loc_cue_word = 'b'
    
    tm_delay = 0.25     # Delay between lead token start and lag token start.
     
     
    # Load hrtf data
    hrtf_path = 'stim/hrtf_kemar_0el.npy'
    dat = np.load(hrtf_path)
    
    # Ensure output directory exists
    if not os.path.exists(out_path_stim):
        os.makedirs(out_path_stim)
        
    # ILD stuff
    ild_bands = psylab.signal.logspace(62.5,8000,8) # psylab.signal.frequency.logspace(70, 2240, 5+1) # 5 1-oct bands
    ild_max_itd_us = 750
    ild_max_ild_db = 32
    ild_wsize = 882
    ild_stepsize = 882
    ild_attack = 882
    ild_exponent = 1.
    
    nofun = np.zeros(ild_max_itd_us)

    itd_ratio = np.linspace(0,ild_max_itd_us,ild_max_itd_us+1)/ild_max_itd_us
    fun = itd_ratio**ild_exponent

    # Here, no positive gain, only attenuation
    ild_table_l = np.concatenate((nofun, -fun)) #np.concatenate((nofun, nofun))
    ild_table_r = np.concatenate(((1-fun)*-1, nofun)) #np.concatenate((nofun, nofun))
    
     
    #t_toks_ae = []
    t_toks_ash = []
    for token in t_s:
        #this_t_data_ae,fs = read(os.path.join(t_path, t_name, token+"ae_"+m_name+".wav"))
        this_t_data_ash,fs = read(os.path.join(t_path, t_name, token+"ash_"+m_name+".wav"))
        #t_toks.append(this_t_data[:,0])
        #t_toks_ae.append(this_t_data_ae)
        t_toks_ash.append(this_t_data_ash)
        
    #m_toks_ae = []
    m_toks_ash = []
    for token in m_s:
        #this_m_data_ae,fs  = read(os.path.join(m_path, m_name, token+"ae_"+m_name+".wav"))
        this_m_data_ash,fs = read(os.path.join(m_path, m_name, token+"ash_"+m_name+".wav"))
        #m_toks.append(this_m_data[:,0])
        #m_toks_ae.append(this_m_data_ae)
        m_toks_ash.append(this_m_data_ash)
        #this_m_data, fs = read(os.path.join(m_path, m_name, token+".wav"))
        #m_toks.append(this_m_data[:,0])
    
    # Get location cue
    #cue_data_ae = t_toks_ae[t_s.index(loc_cue_word)]
    cue_data = t_toks_ash[t_s.index(loc_cue_word)]
    cue_dur = psylab.signal.samp2ms(len(cue_data), fs)
    
    levels = {}
    conditions = list(zip(*list(variables.values())))
    cn = 1
    #for side,masker,itd,ild,az_mag,magnified,control in zip(sides,maskers,itds,ilds,az_mags,magnifieds,controls):
    for condition in conditions:
        for k,v in zip(list(variables.keys()), condition):
            levels[k] = v
        
        this_sf_path_out_sub = "_".join([f"{k}={v}" for k,v in levels.items()])
        this_sf_path_out = os.path.join(out_path_stim, this_sf_path_out_sub)
        if not os.path.exists(this_sf_path_out):
            os.makedirs(this_sf_path_out)
        if not os.path.exists(out_path_data):
            os.makedirs(out_path_data)
            
        for i_t, t_token in enumerate(t_toks_ash): # For each token in possible target tokens....
            for i_m, m_token in enumerate(m_toks_ash): # For each masker in possible masker tokens...
                if i_t == i_m == 0: # both tokens are bash, we skip
                    continue
                for t_lead in range(2): # Whether target will lead or lag for this trial
                    # Make L/R copies of target and makser
                    t_data_l = t_token.copy()
                    t_data_r = t_token.copy()
                    m_data_l = m_token.copy()
                    m_data_r = m_token.copy()
                    c_data_l = cue_data.copy()
                    c_data_r = cue_data.copy()
                    
                    # Apply ITD
                    if 'itd' in levels and levels['itd'] != '0':
                        itd_data = np.zeros(int(int(levels['itd']) / 1000000 * fs))
                    else:
                        itd_data = np.array(())
                        
                    if levels['side'] == 'l':
                        t_data_l = np.concatenate((t_data_l,itd_data))
                        t_data_r = np.concatenate((itd_data,t_data_r))
                        m_data_l = np.concatenate((itd_data,m_data_l))
                        m_data_r = np.concatenate((m_data_r,itd_data))
                        c_data_l = np.concatenate((c_data_l,itd_data))
                        c_data_r = np.concatenate((itd_data,c_data_r))
                        
                    elif levels['side'] == 'r':
                        t_data_l = np.concatenate((itd_data,t_data_l))
                        t_data_r = np.concatenate((t_data_r,itd_data))
                        m_data_l = np.concatenate((m_data_l,itd_data))
                        m_data_r = np.concatenate((itd_data,m_data_r))
                        c_data_l = np.concatenate((itd_data,c_data_l))
                        c_data_r = np.concatenate((c_data_r,itd_data))
                        
                        
                    # Apply Natural ILD
                    
                    if 'az' in levels:
                        # Natural ILDs; apply the magnitude spectrum of the specified ir
                        # MITs KEMAR HRTF set is mono and thus symmetry is assumed
                        # Positive az's are to the right
                        # Use the same ear for both sides but flip the az. 
                        # So for a source az = +60, left = 60, right = 360-60 = 300
                        hrtf_r_ind = int((int(levels['az']) % 360)/5 % (360/5))
                        hrtf_l_ind = int((360-(int(levels['az']) % 360))/5 % (360/5))
                        
                        hrtf_l_data = dat[:, hrtf_l_ind]
                        hrtf_r_data = dat[:, hrtf_r_ind]
        
                        # Apply full hrtf to stim
                        # We do this because we need itds in order
                        # to estimate gain track for magnification
                        c_f_l = np.convolve(c_data_l,hrtf_l_data)
                        c_f_r = np.convolve(c_data_r,hrtf_r_data)
                        t_f_l = np.convolve(t_data_l,hrtf_l_data)
                        t_f_r = np.convolve(t_data_r,hrtf_r_data)
                        m_f_l = np.convolve(m_data_l,hrtf_r_data)
                        m_f_r = np.convolve(m_data_r,hrtf_l_data)
        
                        # Apply magspec from hrtf separately
                        # This will receive magnification
                        c_m_l = ir.apply_mag_spec(c_data_l,hrtf_l_data)
                        c_m_r = ir.apply_mag_spec(c_data_r,hrtf_r_data)
                        t_m_l = ir.apply_mag_spec(t_data_l,hrtf_l_data)
                        t_m_r = ir.apply_mag_spec(t_data_r,hrtf_r_data)
                        m_m_l = ir.apply_mag_spec(m_data_l,hrtf_r_data)
                        m_m_r = ir.apply_mag_spec(m_data_r,hrtf_l_data)
                    
                        # Mix Time Delayed target and masker
                        time_delay_data = np.zeros(int(tm_delay * fs))
                        if t_lead == 0:
                            t_m_l = np.concatenate((time_delay_data,t_m_l))
                            t_m_r = np.concatenate((time_delay_data,t_m_r))
                            m_m_l = np.concatenate((m_m_l,time_delay_data))
                            m_m_r = np.concatenate((m_m_r,time_delay_data))
                            
                            t_f_l = np.concatenate((time_delay_data,t_f_l))
                            t_f_r = np.concatenate((time_delay_data,t_f_r))
                            m_f_l = np.concatenate((m_f_l,time_delay_data))
                            m_f_r = np.concatenate((m_f_r,time_delay_data))
                        elif t_lead == 1:
                            t_m_l = np.concatenate((t_m_l,time_delay_data))
                            t_m_r = np.concatenate((t_m_r,time_delay_data))
                            m_m_l = np.concatenate((time_delay_data,m_m_l))
                            m_m_r = np.concatenate((time_delay_data,m_m_r))
                            
                            t_f_l = np.concatenate((t_f_l,time_delay_data))
                            t_f_r = np.concatenate((t_f_r,time_delay_data))
                            m_f_l = np.concatenate((time_delay_data,m_f_l))
                            m_f_r = np.concatenate((time_delay_data,m_f_r))
                            
                            
                    
                        # Apply Magnification if necessary
                        if levels['mag'] == '0':
                            # No magnification, so t_m and m_m are the final veerion
                            out_l = t_m_l + m_m_l
                            out_r = t_m_r + m_m_r 
                            
                            c_out_l = c_m_l
                            c_out_r = c_m_r
                            
                        elif levels['mag'] == '1':
                            tm_f_l = t_f_l + m_f_l
                            tm_f_r = t_f_r + m_f_r
                            tm_m_l = t_m_l + m_m_l
                            tm_m_r = t_m_r + m_m_r
                            
                            c_f_bp_l = psylab.signal.filter_bank(c_f_l, fs, 4, ild_bands)
                            c_f_bp_r = psylab.signal.filter_bank(c_f_r, fs, 4, ild_bands)
                            tm_f_bp_l = psylab.signal.filter_bank(tm_f_l, fs, 4, ild_bands)
                            tm_f_bp_r = psylab.signal.filter_bank(tm_f_r, fs, 4, ild_bands)
                    
                            c_m_bp_l = psylab.signal.filter_bank(c_m_l, fs, 4, ild_bands)
                            c_m_bp_r = psylab.signal.filter_bank(c_m_r, fs, 4, ild_bands)
                            tm_m_bp_l = psylab.signal.filter_bank(tm_m_l, fs, 4, ild_bands)
                            tm_m_bp_r = psylab.signal.filter_bank(tm_m_r, fs, 4, ild_bands)
                            
                            c_bp_g_l = np.zeros(shape=c_f_bp_l.shape)
                            c_bp_g_r = np.zeros(shape=c_f_bp_r.shape)
                            tm_m_bp_g_l = np.zeros(shape=tm_m_bp_l.shape)
                            tm_m_bp_g_r = np.zeros(shape=tm_m_bp_r.shape)
                            
                            for i in np.arange(len(ild_bands)-1):
        
                                c_itd_track, c_cc = ild.gen_itd_track(c_f_bp_l[:,i], c_f_bp_r[:,i], fs=fs, max_itd_us=ild_max_itd_us,
                                                    wsize=ild_wsize, stepsize=ild_stepsize, thresh=0, return_cc=True)
                                c_gain_l, c_gain_r = ild.gen_gain_tracks(c_f_bp_l[:,i].shape[0], c_itd_track, ild_table_l, ild_table_r,
                                                    max_ild_db=ild_max_ild_db, max_itd_us=ild_max_itd_us, wsize=ild_wsize, stepsize=ild_stepsize,attack=ild_attack,
                                                    nr_cc=None)
                
                                itd_track, cc = ild.gen_itd_track(tm_f_bp_l[:,i], tm_f_bp_r[:,i], fs=fs, max_itd_us=ild_max_itd_us,
                                                    wsize=ild_wsize, stepsize=ild_stepsize, thresh=0, return_cc=True)
                                gain_l, gain_r = ild.gen_gain_tracks(tm_f_bp_l[:,i].shape[0], itd_track, ild_table_l, ild_table_r,
                                                    max_ild_db=ild_max_ild_db, max_itd_us=ild_max_itd_us, wsize=ild_wsize, stepsize=ild_stepsize,attack=ild_attack,
                                                    nr_cc=None)
                                
                                c_bp_g_l[:,i] = ild.apply_gain_track(c_m_bp_l[:,i], c_gain_l)
                                c_bp_g_r[:,i] = ild.apply_gain_track(c_m_bp_r[:,i], c_gain_r)
                                tm_m_bp_g_l[:,i] = ild.apply_gain_track(tm_m_bp_l[:,i], gain_l)
                                tm_m_bp_g_r[:,i] = ild.apply_gain_track(tm_m_bp_r[:,i], gain_r)
                            
                            out_l = tm_m_bp_g_l.sum(axis=1)
                            out_r = tm_m_bp_g_r.sum(axis=1)
                            
                            c_out_l = c_bp_g_l.sum(axis=1)
                            c_out_r = c_bp_g_r.sum(axis=1)
                            # Apply gain track not to hrtf'd data, but to mag-spec'd data
                    else: # no ILDs applied, use ITD only
                        out_l = t_data_l + m_data_l
                        out_r = t_data_r + m_data_r
                        
                        c_out_l = c_data_l
                        c_out_r = c_data_r
                        
                    filename = f"t_{t_s[i_t]}ash_m_{m_s[i_m]}ash_target_lead_{t_lead}_{this_sf_path_out_sub}.wav"
                    this_sf_filepath_out = os.path.join(this_sf_path_out, filename)
                    write(this_sf_filepath_out,np.transpose(np.stack((out_l,out_r))),fs)
                    #elif t_lead == 1: # target lead
        c_filename = f"cue_{this_sf_path_out_sub}.wav"
        this_sf_filepath_out = os.path.join(this_sf_path_out, c_filename)
        write(this_sf_filepath_out,np.transpose(np.stack((c_out_l,c_out_r))),fs)
        
        
        
    
    
if __name__ == '__main__':
    subjectID = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:", ["help", "subjectID="])
    except (getopt.error, msg):
        print("for help use --help")
        sys.exit(2)
    for var, val in opts:
        if var in ("--help", "-h"):
            print(__doc__)
            sys.exit(0)
        elif var in ("--subjectID", "-s"):
            subjectID = val

    if not subjectID:
        subjectID = input('Enter subject ID: ')
    if subjectID == '':
        print("No subject ID found, exiting...")
        sys.exit(0)

    main(subjectID)