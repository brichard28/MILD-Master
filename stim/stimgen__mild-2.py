


import os, sys
import time
import getopt
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal
from pysndfile.sndio import read as sfread
from pysndfile.sndio import write as sfwrite
import psylab
import ir
#import spyral
import ild

__doc__ = """
Generates stimuli for nirs/eeg study looking at ild-based informational masking release.

Word array should have 18 words, 2-4 of which should be colors. 

Each word is 300 ms long, and there should be 300 ms silence between, for 10.8 sec total.

4 stim per condition

Color words were added to the masker stream, as in zhang et al.

"""

def main(subjectID):

    exp = 'mild-2'

    # IMPORTANT: ind=0 of each variable should always be the practice condition
    # az_mags: Which azimuth to take magnitude spectrum to apply
    # magnifieds: Magnify ILDs?
    variables = {
        "side":      ['l',    'l',      'r',      'l',      'r',      'l',      'r',      'l',      'r',      ],
        "masker":    ['none', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', ],
        "itd":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
        "ild":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
        "az":        ['0',    '10',     '10',     '10',     '10',     '30',     '30',     '30',     '30',     ],
        "mag":       ['0',    '0',      '0',      '1',      '1',      '0',      '0',      '1',      '1',      ],
        "control":   ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
                }

    #variables = {
    #    "side":      ['l',    'l',      'r',      'l',      'r',      'l',      'r',      'l',      'r',      ],
    #    "masker":    ['none', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', 'speech', ],
    #    "itd":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
    #    "ild":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
    #    "az":        ['0',    '10',     '10',     '10',     '10',     '30',     '30',     '30',     '30',     ],
    #    "mag":       ['0',    '0',      '0',      '1',      '1',      '0',      '0',      '1',      '1',      ],
    #    "control":   ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
    #            }

    n = 15             # # of runs, Zhang used 6. Here we use 4 since side is an IV (8 total, 4 per side)
    nwords = 18        # # of words in a list, Zhang used 25 (15 sec) (18 * (.3word + .3 silence) = 10.8s)
    ncolors = [3, 4]   # list of the range of possible numbers of colors, Zhang used 3-5
    word_dur = .3
    silence = .3       # Amount of silence between words
    masker_delay = .15 # Delay between target start and masker start. Has been effectively zero <= nirs-im-8
                       # Won't change overall length. IE, will eat into silence var above; so if word dur == .3
                       # silence == .3, and masker_delay == .15, there will be .15 s of silence between offset
                       # of a masker and onset of next target
    t_talker = 'Bob' # Target Talker (Bob or Mike)
    mask_talker = 'Bob' # Masker Talker

    target_first = -1  # if masker_delay > 0; 1 = target first; 0 = masker first; -1 = random

    # location cue; a tone to indicate which location to listen to.
    loc_cue_f = 500
    loc_cue_dur = 500
    loc_cue_isi = 2000 #[500,2500] # ms; use a number for fixed, a list of size 2 for random

    # Soundfile paths


    #sf_path_out = os.path.expanduser(f'~/work/Python/stim/bobmike/{exp}/s_{subjectID}')
    #data_path_out = os.path.expanduser(f'~/work/Python/data/{exp}')
    sf_path_out = f'./stim/bobmike/{exp}/s_{subjectID}'
    data_path_out = f'./data/{exp}'
    sf_path_in = './stim/bobmike/unprocessed'
    sf_fs = 44100
    cue_word = 'shoe'

    # Token info
    colors = ['blue', 'green', 'red',    'white']
    objects = ['bag', 'card',  'chairs', 'desks', 'glove', 'hat', 'pen', 'shoe', 'sock', 'spoons', 'tables', 'toy']

    snr = 0

    # ILD stuff
    ild_bands= psylab.signal.logspace(62.5,8000,8) # psylab.signal.frequency.logspace(70, 2240, 5+1) # 5 1-oct bands
    ild_max_itd_us = 1000
    ild_max_ild_db = 32

    ild_wsize = 882
    ild_stepsize = 441
    ild_attack = 0
    ild_exponent = 1.

    nofun = np.zeros(ild_max_itd_us)

    itd_ratio = np.linspace(0,ild_max_itd_us,ild_max_itd_us+1)/ild_max_itd_us
    fun = itd_ratio**ild_exponent

    # Here, no positive gain, only attenuation
    ild_table_l = np.concatenate((nofun, -fun)) #np.concatenate((nofun, nofun))
    ild_table_r = np.concatenate(((1-fun)*-1, nofun)) #np.concatenate((nofun, nofun))

    # Vocoder parameters From Zhang, Alamatsaz & Ihlefeld 2021
    # Use filtfilt for 72 dB/oct and 0 phase shift
    v_channel_fcs = psylab.signal.logspace(300,10000,16)
    v_channel_cos = np.zeros((2, 16))

    # format: (b,a)/(:)/(hi,lo)/(channel)
    v_channel_coefs = np.zeros((2,7,2,16))
    for i in np.arange(v_channel_fcs.size):
        f1,f2 = psylab.signal.oct2f(v_channel_fcs[i], 1/20.)
        v_channel_cos[:,i] = f1,f2
        v_channel_coefs[:,:,0,i] = scipy.signal.butter(6, f1/sf_fs/2., btype='highpass')
        v_channel_coefs[:,:,1,i] = scipy.signal.butter(6, f2/sf_fs/2., btype='lowpass')

    hrtf_path = 'stim/hrtf_kemar_0el.npy'
    dat = np.load(hrtf_path)

    if not os.path.exists(sf_path_out):
        os.makedirs(sf_path_out)
    if not os.path.exists(data_path_out):
        os.makedirs(data_path_out)
#    fid_color_times_t = open(os.path.join(data_path_out, f'target_color_times__s_{subjectID}.txt'), 'w')
#    fid_color_times_m = open(os.path.join(data_path_out, f'masker_color_times__s_{subjectID}.txt'), 'w')
#    fid_target_words = open(os.path.join(data_path_out, f'target_words__s_{subjectID}.txt'), 'w')
    fid_all_times = open(os.path.join(data_path_out, f'{exp}__s_{subjectID}__{time.strftime("%Y-%m-%d")}__word_times.csv'), 'w')

    #cue_data = psylab.signal.tone(loc_cue_f, sf_fs, loc_cue_dur, amp=.3)
    #cue_data = psylab.signal.ramps(cue_data,sf_fs)
    cue_data,fs,_enc = sfread(os.path.join(sf_path_in,  f"{t_talker.lower()}_object",   f"{cue_word}_short.wav"))

    rng = np.random.default_rng()

    header = "Condition,Run,Type"
    for k in np.arange(nwords):
        header += f",Word{k},Time{k}"
    header += "\n"
    fid_all_times.write(header)

    #for masker,itd,ild in zip(maskers, itds, ilds):
    pract = True
    levels = {}
    conditions = list(zip(*list(variables.values())))
    cn = 1
    #for side,masker,itd,ild,az_mag,magnified,control in zip(sides,maskers,itds,ilds,az_mags,magnifieds,controls):
    for condition in conditions:
        for k,v in zip(list(variables.keys()), condition):
            levels[k] = v
        if pract:
            this_sf_path_out_sub = 'practice'
            pract = False
        else:
            this_sf_path_out_sub = "_".join([f"{k}={v}" for k,v in levels.items()])
            #this_sf_path_out_sub = f"m_{masker}__ild_{ild}__itd_{itd}__targ_{side}__control_{control}" # For textfiles
        this_sf_path_out = os.path.join(sf_path_out, this_sf_path_out_sub)
        if not os.path.exists(this_sf_path_out):
            os.makedirs(this_sf_path_out)
        print(f"{cn} of {len(conditions)}")
        cn += 1

        if isinstance(loc_cue_isi,list):
            # list means random in specified range
            cue_isi_dur = rng.uniform(loc_cue_isi[0], loc_cue_isi[1], size=1)[0]
        else:
            cue_isi_dur = loc_cue_isi
        cue_isi_samp = psylab.signal.ms2samp(cue_isi_dur, sf_fs)
        cue_isi_data = np.zeros(cue_isi_samp)

        if levels['itd'] != '0':
            itd_data = np.zeros(int(int(levels['itd']) / 1000000 * sf_fs))
        else:
            itd_data = np.array(())

        lead_silence_data_pre = np.zeros(0)
        lead_silence_data_post = np.zeros(np.int32(sf_fs*silence)-len(itd_data))
        lag_silence_data_pre = np.zeros(np.int32(sf_fs*masker_delay))
        lag_silence_data_post = np.zeros(np.int32(sf_fs*(silence-masker_delay)))

        # Generate lists of inds, ensure: 
        #   - each ind is unique across lists
        #   - no targ ind is same as previous
        for k in np.arange(n):
            filename = f"{k}_{this_sf_path_out_sub}.wav"
            this_sf_filepath_out = os.path.join(this_sf_path_out, filename)
            this_sf_filepath_out_sub = os.path.join(this_sf_path_out_sub, filename)
            print(f"  {filename}")

            # used for writing words and times to fid_all_times. In S
            curr_time = ((word_dur*1000) + cue_isi_dur) / 1000

            t_str = f"{this_sf_path_out_sub},{k},Target,"
            mask_str = f"{this_sf_path_out_sub},{k},Masker,"
            nobjectlists = nwords // len(objects)
            robjectlists = nwords % len(objects)
            objects_25 = []
            for j in range(nobjectlists):
                objects_25.extend(objects)
            for j in range(robjectlists):
                objects_25.append(objects[np.random.choice(len(objects))])
            #objects_25 = objects + objects + [objects[np.random.choice(len(objects))]]
            prev_targ = None
            t_arr =  []
            #m1_arr = []
            m2_arr = []
            sf_t = []
            #sf_m1 = []
            sf_m2 = []
            for j in np.arange(nwords):
                this_i = np.random.choice(nwords, 3, replace=False)
                if prev_targ:
                    while this_i[0] == prev_targ:
                        this_i = np.random.choice(nwords, 3, replace=False)
                prev_targ = this_i[0]
                t_arr.append(this_i[0])
                #m1_arr.append(this_i[1])
                m2_arr.append(this_i[2])

            t_words =  []
            #m1_words = []
            m2_words = []

            t_times = []
            #m1_times = []
            m2_times = []

            t_data =  np.zeros(sf_fs)
            #m1_data = np.zeros(sf_fs)
            m2_data = np.zeros(sf_fs)

            colors_n_t = np.random.randint(ncolors[0],ncolors[1]+1) # between 3 - 5 colors
            colors_n_m = np.random.randint(ncolors[0],ncolors[1]+1)
            colors_ind = np.random.choice(nwords, colors_n_t+colors_n_m, replace=False)
            colors_i_t = colors_ind[:colors_n_t]
            colors_i_m = colors_ind[colors_n_t:]
            color_times_t = []
            color_times_m = []

            for i in np.arange(len(t_arr)):

                # TODO: this_t_i and t_arr do not seem to be used, as 
                # we are choosing from objects or colors at random below
                this_t_i =  t_arr[i]
                #this_m1_i = m1_arr[i]
                this_m2_i = m2_arr[i]
                if i in colors_i_t:
                    this_t_word = colors[np.random.randint(len(colors))]
                    this_t_kind = f'{t_talker.lower()}_color'
                    color_times_t.append(str(np.round(i*(word_dur+silence), 3))) # 600 ms between words
                else:
                    #this_t_word = objects[np.random.randint(len(objects))]
                    this_t_word = objects_25[this_t_i]
                    this_t_kind = f'{t_talker.lower()}_object'
                t_words.append(this_t_word)

                if levels['masker'] == 'speech':
                    if i in colors_i_m:
                        #this_m1_word = colors[np.random.randint(len(colors))]
                        #this_m1_kind = f'{mask_talker.lower()}_color'
                        this_m2_word = colors[np.random.randint(len(colors))]
                        this_m2_kind = f'{mask_talker.lower()}_color'
                        color_times_m.append(str(np.round(i*(word_dur+silence+masker_delay), 3))) # 600 ms between words
                    else:
                        #this_m1_word = objects_25[this_m1_i]
                        #this_m1_kind = f'{mask_talker.lower()}_object'
                        this_m2_word = objects_25[this_m2_i]
                        this_m2_kind = f'{mask_talker.lower()}_object'
                elif levels['masker'] == 'noise':
                    #this_m1_word = "noise"
                    #this_m1_kind = ""
                    this_m2_word = "noise"
                    this_m2_kind = ""
                else:
                    #this_m1_word = "none"
                    #this_m1_kind = ""
                    this_m2_word = "none"
                    this_m2_kind = ""

                #m1_words.append(this_m1_word)
                m2_words.append(this_m2_word)
                this_t_sf = os.path.join(sf_path_in,  this_t_kind,   f"{this_t_word}_short.wav")
                sf_t.append(this_t_sf)
                this_t_data,fs,_enc = sfread(this_t_sf)
            #    target_first = -1  # if masker_delay > 0; 1 = target first; 0 = masker first; -1 = random
                if target_first in [0, 1]:
                    t_f = target_first
                else:
                    t_f = np.random.choice([0,1])
                if t_f == 1:
                    t_data = np.concatenate((t_data, lead_silence_data_pre, this_t_data, lead_silence_data_post))
                    this_t_time = np.round(curr_time, 4)
                else:
                    t_data = np.concatenate((t_data, lag_silence_data_pre, this_t_data, lag_silence_data_post))
                    this_t_time = np.round(curr_time + masker_delay, 4)

                if this_m2_word == 'noise':
                    #this_m1_data = psylab.signal.white(this_t_data.shape[0])
                    this_m2_data = psylab.signal.white(this_t_data.shape[0])
                    #this_m1_data *= psylab.signal.rms(this_t_data) / psylab.signal.rms(this_m1_data)
                    this_m2_data *= psylab.signal.rms(this_t_data) / psylab.signal.rms(this_m2_data)
                elif this_m2_word == 'none':
                    #this_m1_data = np.zeros(this_t_data.shape[0])
                    this_m2_data = np.zeros(this_t_data.shape[0])
                else:
                    #print(this_m1_word)
                    #this_m1_sf = os.path.join(sf_path_in, this_m1_kind, f"{this_m1_word}_short.wav")
                    #sf_m1.append(this_m1_sf)
                    #this_m1_data,fs,_enc = sfread(this_m1_sf)

                    this_m2_sf = os.path.join(sf_path_in, this_m2_kind, f"{this_m2_word}_short.wav")
                    sf_m2.append(this_m2_sf)
                    this_m2_data,fs,_enc = sfread(this_m2_sf)

                if t_f == 1:
                    #m1_data = np.concatenate((m1_data, lag_silence_data_pre, this_m1_data, lag_silence_data_post))
                    m2_data = np.concatenate((m2_data, lag_silence_data_pre, this_m2_data, lag_silence_data_post))
                    this_mask_time = np.round(curr_time + masker_delay, 4)
                else:
                    #m1_data = np.concatenate((m1_data, lead_silence_data_pre, this_m1_data, lead_silence_data_post))
                    m2_data = np.concatenate((m2_data, lead_silence_data_pre, this_m2_data, lead_silence_data_post))
                    this_mask_time = np.round(curr_time, 4)

                t_times.append(this_t_time)
                m2_times.append(this_mask_time)
                #t_times += f",{this_t_word},{this_t_time}"
                #mask_times += f",{this_m2_word},{this_mask_time}"

                curr_time += word_dur + silence

            # Done with stim creation; begin processing

            c_data_l = cue_data.copy()
            c_data_r = cue_data.copy()
            t_data_l = t_data.copy()
            t_data_r = t_data.copy()
            #m1_data_l = m1_data.copy()
            #m1_data_r = m1_data.copy()
            m2_data_l = m2_data.copy()
            m2_data_r = m2_data.copy()
            #if ild == 'inf':
            #    t_data_r = np.zeros(t_data_r.shape)
            #elif ild != '0':
            #    t_data_r = psylab.signal.atten(t_data_r, int(ild))

            if levels['masker'] == 'none':
                t_data_l = psylab.signal.atten(t_data_l, -6) # W no maskers, the overall level is low
                t_data_r = psylab.signal.atten(t_data_r, -6) # Probably because there are only 7 bands total
                #m1_data_l = np.zeros(t_data.shape)
                #m1_data_r = np.zeros(t_data.shape)
                m2_data_l = np.zeros(t_data.shape)
                m2_data_r = np.zeros(t_data.shape)
            else:
                #m1_data_l = psylab.signal.atten(m1_data_l,snr) #+3) # +3 to account for 2 maskers
                #m1_data_r = psylab.signal.atten(m1_data_r,snr) #+3)
                m2_data_l = psylab.signal.atten(m2_data_l,snr) #+3)
                m2_data_r = psylab.signal.atten(m2_data_r,snr) #+3)

                ## Use 500-us ITD for cue, because it is a tone
                ## and other cues don't always lateralize it sufficiently
                #cue_data_l = np.concatenate((cue_data_l, np.zeros(int(500 / 1000000 * sf_fs)) ))
                #cue_data_r = np.concatenate((np.zeros(int(500 / 1000000 * sf_fs)), cue_data_r))

                # Apply any ITD
                # 0 is handled correctly
                # We only use 1 masker for nirs-im-8, so ignore m1 and apply spatial cues to targ
                c_data_l = np.concatenate((c_data_l, itd_data))
                c_data_r = np.concatenate((itd_data, c_data_r))
                t_data_l = np.concatenate((t_data_l, itd_data))
                t_data_r = np.concatenate((itd_data, t_data_r))
                #m1_data_l = np.concatenate((itd_data, m1_data_l))
                #m1_data_r = np.concatenate((m1_data_r, itd_data))
                m2_data_l = np.concatenate((itd_data, m2_data_l))
                m2_data_r = np.concatenate((m2_data_r, itd_data))

                if 'ild' in levels:
                    if levels['ild'] == 'inf':
                        c_data_r = np.zeros(len(c_data_r))
                        t_data_r = np.zeros(len(t_data_r))
                        #m1_data_r = np.zeros(len(m1_data_r))
                        m2_data_l = np.zeros(len(m2_data_l))
                    
                    else:
                        this_ild = float(levels['ild']) / 2
                        c_data_l = psylab.signal.atten(c_data_l,-this_ild)
                        c_data_r = psylab.signal.atten(c_data_r,this_ild)
                        t_data_l = psylab.signal.atten(t_data_l,-this_ild)
                        t_data_r = psylab.signal.atten(t_data_r,this_ild)
                        #m1_data_l = psylab.signal.atten(m1_data_l,this_ild)
                        #m1_data_r = psylab.signal.atten(m1_data_r,-this_ild)
                        m2_data_l = psylab.signal.atten(m2_data_l,-this_ild)
                        m2_data_r = psylab.signal.atten(m2_data_r,this_ild)

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
                    c_f_l = np.convolve(c_data_l,hrtf_l_data)
                    c_f_r = np.convolve(c_data_r,hrtf_r_data)
                    t_f_l = np.convolve(t_data_l,hrtf_l_data)
                    t_f_r = np.convolve(t_data_r,hrtf_r_data)
                    #m1_f_l = np.convolve(m1_data_l,hrtf_l_data)
                    #m1_f_r = np.convolve(m1_data_r,hrtf_r_data)
                    m2_f_l = np.convolve(m2_data_l,hrtf_r_data)
                    m2_f_r = np.convolve(m2_data_r,hrtf_l_data)

                    c_m_l = ir.apply_mag_spec(c_data_l,hrtf_l_data)
                    c_m_r = ir.apply_mag_spec(c_data_r,hrtf_r_data)
                    t_m_l = ir.apply_mag_spec(t_data_l,hrtf_l_data)
                    t_m_r = ir.apply_mag_spec(t_data_r,hrtf_r_data)
                    #m1_m_l = ir.apply_mag_spec(m1_data,hrtf_r_data)
                    #m1_m_r = ir.apply_mag_spec(m1_data,hrtf_l_data)
                    m2_m_l = ir.apply_mag_spec(m2_data_l,hrtf_r_data)
                    m2_m_r = ir.apply_mag_spec(m2_data_r,hrtf_l_data)

                    # Mix
                    t_m_l, t_m_r, m2_m_l, m2_m_r, t_f_l, t_f_r, m2_f_l, m2_f_r = psylab.signal.zeropad(
                    t_m_l, t_m_r, m2_m_l, m2_m_r, t_f_l, t_f_r, m2_f_l, m2_f_r)

                    # Filterbank
                    c_m_bp_l = psylab.signal.filter_bank(c_m_l, fs, 4, ild_bands)
                    c_m_bp_r = psylab.signal.filter_bank(c_m_r, fs, 4, ild_bands)
                    t_m_bp_l = psylab.signal.filter_bank(t_m_l, fs, 4, ild_bands)
                    t_m_bp_r = psylab.signal.filter_bank(t_m_r, fs, 4, ild_bands)
                    #m1_m_bp_l = psylab.signal.filter_bank(m1_m_l, fs, 4, ild_bands)
                    #m1_m_bp_r = psylab.signal.filter_bank(m1_m_r, fs, 4, ild_bands)
                    m2_m_bp_l = psylab.signal.filter_bank(m2_m_l, fs, 4, ild_bands)
                    m2_m_bp_r = psylab.signal.filter_bank(m2_m_r, fs, 4, ild_bands)

                    if levels['mag'] == '0':
                        c_data_l = c_m_bp_l.sum(axis=1)
                        c_data_r = c_m_bp_r.sum(axis=1)
                        t_data_l = t_m_bp_l.sum(axis=1)
                        t_data_r = t_m_bp_r.sum(axis=1)
                        #m1_data_l = m1_m_bp_l.sum(axis=1)
                        #m1_data_r = m1_m_bp_r.sum(axis=1)
                        m2_data_l = m2_m_bp_l.sum(axis=1)
                        m2_data_r = m2_m_bp_r.sum(axis=1)
                    else:
                        tm_f_l = t_f_l + m2_f_l #+ m2_f_l
                        tm_f_r = t_f_r + m2_f_r #+ m2_f_r
                        tm_m_l = t_m_l + m2_m_l #+ m2_m_l
                        tm_m_r = t_m_r + m2_m_r #+ m2_m_r

                        c_f_bp_l = psylab.signal.filter_bank(c_f_l, fs, 4, ild_bands)
                        c_f_bp_r = psylab.signal.filter_bank(c_f_r, fs, 4, ild_bands)
                        tm_f_bp_l = psylab.signal.filter_bank(tm_f_l, fs, 4, ild_bands)
                        tm_f_bp_r = psylab.signal.filter_bank(tm_f_r, fs, 4, ild_bands)

                        c_m_bp_g_l = np.zeros(shape=c_m_bp_l.shape)
                        c_m_bp_g_r = np.zeros(shape=c_m_bp_r.shape)
                        t_m_bp_g_l = np.zeros(shape=t_m_bp_l.shape)
                        t_m_bp_g_r = np.zeros(shape=t_m_bp_r.shape)
                        #m1_m_bp_g_l = np.zeros(shape=m1_m_bp_l.shape)
                        #m1_m_bp_g_r = np.zeros(shape=m1_m_bp_r.shape)
                        m2_m_bp_g_l = np.zeros(shape=m2_m_bp_l.shape)
                        m2_m_bp_g_r = np.zeros(shape=m2_m_bp_r.shape)

                        # Loop through bands
                        for i in np.arange(len(ild_bands)-1):

                            c_itd_track, c_cc = ild.gen_itd_track(c_f_bp_l[:,i], c_f_bp_r[:,i], fs=fs, max_itd_us=ild_max_itd_us,
                                                wsize=ild_wsize, stepsize=ild_stepsize, thresh=0, return_cc=True)
                            c_gain_l, c_gain_r = ild.gen_gain_tracks(c_f_bp_l[:,i].shape[0], c_itd_track, ild_table_l, ild_table_r,
                                                max_ild_db=ild_max_ild_db, max_itd_us=ild_max_itd_us, wsize=ild_wsize, stepsize=ild_stepsize,
                                                nr_cc=None)
            
                            itd_track, cc = ild.gen_itd_track(tm_f_bp_l[:,i], tm_f_bp_r[:,i], fs=fs, max_itd_us=ild_max_itd_us,
                                                wsize=ild_wsize, stepsize=ild_stepsize, thresh=0, return_cc=True)
                            gain_l, gain_r = ild.gen_gain_tracks(tm_f_bp_l[:,i].shape[0], itd_track, ild_table_l, ild_table_r,
                                                max_ild_db=ild_max_ild_db, max_itd_us=ild_max_itd_us, wsize=ild_wsize, stepsize=ild_stepsize,
                                                nr_cc=None)
                            # Apply gain track not to hrtf'd data, but to mag-spec'd data
                            # This allows magnified ILD cues to be applied without itd cues present
                            # Which is needed for unvocoded stimuli
                            c_m_bp_g_l[:,i] = ild.apply_gain_track(c_m_bp_l[:,i], c_gain_l)
                            c_m_bp_g_r[:,i] = ild.apply_gain_track(c_m_bp_r[:,i], c_gain_r)
                            t_m_bp_g_l[:,i] = ild.apply_gain_track(t_m_bp_l[:,i], gain_l)
                            t_m_bp_g_r[:,i] = ild.apply_gain_track(t_m_bp_r[:,i], gain_r)
                            #m1_m_bp_g_l[:,i] = ild.apply_gain_track(m1_m_bp_l[:,i], gain_l)
                            #m1_m_bp_g_r[:,i] = ild.apply_gain_track(m1_m_bp_r[:,i], gain_l)
                            m2_m_bp_g_l[:,i] = ild.apply_gain_track(m2_m_bp_l[:,i], gain_l)
                            m2_m_bp_g_r[:,i] = ild.apply_gain_track(m2_m_bp_r[:,i], gain_r)

                        c_data_l = c_m_bp_g_l.sum(axis=1)
                        c_data_r = c_m_bp_g_r.sum(axis=1)
                        t_data_l = t_m_bp_g_l.sum(axis=1)
                        t_data_r = t_m_bp_g_r.sum(axis=1)
                        #m1_data_l = m1_m_bp_g_l.sum(axis=1)
                        #m1_data_r = m1_m_bp_g_r.sum(axis=1)
                        m2_data_l = m2_m_bp_g_l.sum(axis=1)
                        m2_data_r = m2_m_bp_g_r.sum(axis=1)

            # Add location cue; combine target and masker
            sig_l = np.concatenate((c_data_l, cue_isi_data, t_data_l + m2_data_l))
            sig_r = np.concatenate((c_data_r, cue_isi_data, t_data_r + m2_data_r))
            if levels['side'] == 'r':
                sig_out = np.column_stack((sig_r, sig_l))
            else:
                sig_out = np.column_stack((sig_l, sig_r))

            sfwrite(this_sf_filepath_out, sig_out, fs, format='wav')
            #t_str += f"{','.join(t_words)},{','.join(map(str, t_times))}\n"
            #mask_str += f"{','.join(m2_words)},{','.join(map(str, m2_times))}\n"
            t_str += ','.join([','.join(map(str, i)) for i in zip(t_words, t_times)]) + "\n"
            mask_str += ','.join([','.join(map(str, i)) for i in zip(m2_words, m2_times)]) + "\n"
            fid_all_times.write(t_str)
            fid_all_times.write(mask_str)


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
