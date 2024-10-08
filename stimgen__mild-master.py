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

__doc__ = """
Generates stimuli for eeg study looking at ild-based informational masking release
and magnified ILD cues. Conditions include ITD50, ITD500, ILDAZ10, ILDAZ10MAG.

"""

def main(subjectID):

    exp = 'mild-master'

    # IMPORTANT: ind=0 of each variable should always be the practice condition
    # variables = {
    #     "side":      ['l',    'l',      'r',      'l',      'r',      'l',      'r',      'l',      'r',      ],
    #     "itd":       ['0',    '50',     '50',     '500',    '500',    '0',      '0',      '0',      '0',      ],
    #     "az":        ['0',    '0',      '0',      '0',      '0',      '15',      '15',      '15',      '15',      ],
    #     "mag":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '1',      '1',      ],
    #             }
    
    # variables for testing various ild azimuths
    # variables = {
    #     "side":      ['l',    'l', 'l',      'r',      'l',      'r',      'l',      'r',      'l',      'r',     'r' ],
    #     "itd":       ['0',    '0', '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0'],
    #     "az":        ['0',    '0', '0',      '5',      '5',      '10',      '10',      '30',     '30', '60', '60'  ],
    #     "mag":       ['0',    '0', '1',      '0',      '1',      '0',      '1',      '0',      '1',      '0',      '1'],
    #             }
    
    variables = {
        "side":      ['l', 'l',    'l' ],
        "itd":       ['0', '0',    '0'],
        "az":        ['0', '5',    '60'],
        "mag":       ['0', '1',    '1'],
                }
    
    # variables for ITD pilot
    # variables = {
    #     "side":      ['l',    'l',      'r',      'l',      'r',      'l',      'r',      'l',      'r',      ],
    #     "itd":       ['0',    '50',     '50',     '100',    '100',    '200',    '200',    '400',    '400',      ],
    #     "az":        ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
    #     "mag":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],
    #             }
    
        #"masker":    ['none', 'sp',     'sp',     'sp',     'sp',     'sp',     'sp',     'sp',     'sp',     'sp',     'sp', ],
        #"ild":       ['0',    '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      '0',      ],

    fs = 44100          # fs
    n = 1             # # of trials
    ntokens = 7         # # of tokens in a trial
    snr = 0
    t_path = './stim/bashdashgash'
    m_path = './stim/bashdashgash'
    out_path_stim = f'./stim/{exp}/s_{subjectID}'
    out_path_data = f'./data/{exp}'

    t_name = 'ben'
    m_name = 'ben'

    t_s = ['bash', 'dash', 'gash']
    m_s = ['bash', 'dash', 'gash']
    
    tm_delay = 0.25     # Delay between lead token start and lag token start.
    tt_delay = 0.05     # Delay between lag token end and next lead target start.
    t_first = -1        # if tm_delay > 0; 1 = target first; 0 = masker first; -1 = random

    """
        L  Ba                Ga
           |__________|      |__________|
        R           Da                Ba          
                    |__________|      |__________|

           |        | |      | |      | |
     Time  0      300 350  600 650  900 950
 tm_delay  |________|
 tt_delay             |______|
    """
    loc_cue_isi = 1.5
    loc_cue_word = 'bash'

    pract_atten = 5

    # Read in tokens
    t_toks = []
    for token in t_s:
        #this_t_data,fs,_enc = sfread(os.path.join(t_path, t_name, token+".wav"))
        #t_toks.append(this_t_data[:,0])
        
        this_t_data, fs = read(os.path.join(t_path, t_name, token+".wav"))
        t_toks.append(this_t_data[:,0])

    m_toks = []
    for token in m_s:
        #this_m_data,fs,_enc = sfread(os.path.join(m_path, m_name, token+".wav"))
        #m_toks.append(this_m_data[:,0])
        
        this_m_data, fs = read(os.path.join(m_path, m_name, token+".wav"))
        m_toks.append(this_m_data[:,0])

    # Get location cue
    cue_data = t_toks[t_s.index(loc_cue_word)]
    cue_dur = psylab.signal.samp2ms(len(cue_data), fs)

    tm_delay_data = np.zeros(psylab.signal.ms2samp(tm_delay*1000, fs))
    if tt_delay < 0:
        tt_delay_data = np.array(())
    else:
        tt_delay_data = np.zeros(psylab.signal.ms2samp(tt_delay*1000, fs))

    # ILD stuff
    ild_bands= psylab.signal.logspace(62.5,8000,8) # psylab.signal.frequency.logspace(70, 2240, 5+1) # 5 1-oct bands
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


    ## Vocoder stuff
    #voc_ch = 12
    #voc_f_lo = 125
    #voc_f_hi = 8000
    #voc_bands = psylab.signal.logspace(voc_f_lo,voc_f_hi,voc_ch+1)
    #
    ## format: (b,a)/(:)/(hi,lo)/(channel)
    #v_channel_coefs = np.zeros((2,7,2,voc_ch))
    #for i in np.arange(voc_bands.size-1):
    #    f1 = voc_bands[i]
    #    f2 = voc_bands[i+1]
    #    v_channel_coefs[:,:,0,i] = scipy.signal.butter(6, f1/fs/2., btype='highpass')
    #    v_channel_coefs[:,:,1,i] = scipy.signal.butter(6, f2/fs/2., btype='lowpass')

    # Load hrtf data
    hrtf_path = 'stim/hrtf_kemar_0el.npy'
    dat = np.load(hrtf_path)

    # Ensure output directory exists
    if not os.path.exists(out_path_stim):
        os.makedirs(out_path_stim)

    # Use a bool for practice condition. 
    # Enter main loop and check, then immediately set to false
    pract = True
    levels = {}
    conditions = list(zip(*list(variables.values())))
    cn = 1
    #for side,masker,itd,ild,az_mag,magnified,control in zip(sides,maskers,itds,ilds,az_mags,magnifieds,controls):
    
    for condition in conditions:
        all_itd_tracks = []
        for k,v in zip(list(variables.keys()), condition):
            levels[k] = v
        if pract:
            this_sf_path_out_sub = 'practice'
            #pract = False
        else:
            this_sf_path_out_sub = "_".join([f"{k}={v}" for k,v in levels.items()])
            #this_sf_path_out_sub = f"m_{masker}__ild_{ild}__itd_{itd}__targ_{side}__control_{control}" # For textfiles
        this_sf_path_out = os.path.join(out_path_stim, this_sf_path_out_sub)
        if not os.path.exists(this_sf_path_out):
            os.makedirs(this_sf_path_out)
        if not os.path.exists(out_path_data):
            os.makedirs(out_path_data)

        print(f"{cn} of {len(conditions)}")
        cn += 1

        if isinstance(loc_cue_isi,list):
            # list means random in specified range
            cue_isi_dur = np.random.uniform(loc_cue_isi[0], loc_cue_isi[1], size=1)[0]
        else:
            cue_isi_dur = loc_cue_isi
        c_isi_samp = psylab.signal.ms2samp(cue_isi_dur*1000, fs)
        c_isi_data = np.zeros(c_isi_samp)

        if 'itd' in levels and levels['itd'] != '0':
            itd_data = np.zeros(int(int(levels['itd']) / 1000000 * fs))
        else:
            itd_data = np.array(())

        for k in np.arange(n):
            t_words =  []
            #m1_words = []
            m2_words = []

            t_times = []
            #m1_times = []
            m2_times = []

            # if t_first == 1:
            #     t_data = np.array(())
            #     m2_data = tm_delay_data
            # elif t_first == 0:
            #     t_data = tm_delay_data
            #     m2_data = np.array(())
            # else:
            #     # TODO: add jitter condition
            #     # choose a random leader
            #     t_lead = random.randint(0, 1) # 1 target leads, 0 masker leads
            #     if t_lead == 1:
            #         t_data = np.array(())
            #         m2_data = tm_delay_data
            #     else:
            #         t_data = tm_delay_data
            #         m2_data = np.array(())
                
            t_data = np.array(())
            m2_data = np.array(())
            # Create order of bash-dash-gashs
            min_number_bashs = 1
            max_number_bashs = 2
            total_number_bashs = min_number_bashs + max_number_bashs
            n_close_idx = 2
            bash_index = 0
            # First, create just dashs and gashs 
            t_inds = np.random.randint(1, high=3, size=ntokens)
            m2_inds = np.random.randint(1, high= 3, size=ntokens)
            
            # Next, find indices of bashs for target
            possible_target_bash_indices = np.arange(0,len(t_inds))#np.arange(1,len(t_inds)-1) # cannot have bash on first or last one
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
            # Assign bash indices to existing inds
            t_inds[bash_indices_target] = bash_index
            m2_inds[bash_indices_masker] = bash_index


            #while True:
            #    m2_inds = np.random.randint(len(m_s), size=ntokens)
            #    if not (t_inds==m2_inds).any():
            #        break
            
            for j in np.arange(ntokens):
                if t_first == 0: # DONT USE THIS
                    t_data = np.append(t_data, tm_delay_data)
                elif t_first == 1: # DONT USE THIS
                    m2_data = np.append(m2_data, tm_delay_data)
                else: # USE THIS ONE
                    t_lead = random.randint(0, 1) # 1 target leads, 0 masker leads
                    if t_lead == 1: # target leads
                        m2_data = np.append(m2_data, tm_delay_data)
                        t_times.append(len(t_data))
                        m2_times.append(len(m2_data))
                        
                        t_data = np.append(t_data, t_toks[t_inds[j]])
                        t_words.append(t_s[t_inds[j]])
                        m2_data = np.append(m2_data, m_toks[m2_inds[j]])
                        m2_words.append(m_s[m2_inds[j]])
                        t_data = np.append(t_data, tm_delay_data)
                    elif t_lead == 0:
                        t_data = np.append(t_data, tm_delay_data)
                        t_times.append(len(t_data))
                        m2_times.append(len(m2_data))
                        
                        
                        t_data = np.append(t_data, t_toks[t_inds[j]])
                        t_words.append(t_s[t_inds[j]])
                        m2_data = np.append(m2_data, m_toks[m2_inds[j]])
                        m2_words.append(m_s[m2_inds[j]])
                        m2_data = np.append(m2_data, tm_delay_data)
                
                t_data = np.append(t_data, tt_delay_data)
                m2_data = np.append(m2_data, tt_delay_data)

            filename = f"{k}__{this_sf_path_out_sub}__l={'-'.join(t_words)}__r={'-'.join(m2_words)}.wav"
            this_sf_filepath_out = os.path.join(this_sf_path_out, filename)
            this_sf_filepath_out_sub = os.path.join(this_sf_path_out_sub, filename)
            print(f"  {filename}")

            # used for writing words and times to fid_all_times. In S
            curr_time = (cue_dur + cue_isi_dur) / 1000

            # When side is r, the masker is now the target and visa versa
            if levels['side'] == 'r':
                t_str = f"{k}_{this_sf_path_out_sub},{filename},Masker,"
                m_str = f"{k}_{this_sf_path_out_sub},{filename},Target,"
            else:
                t_str = f"{k}_{this_sf_path_out_sub},{filename},Target,"
                m_str = f"{k}_{this_sf_path_out_sub},{filename},Masker,"

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

            if pract: #if levels['masker'] == 'none':
                c_data_l = psylab.signal.atten(c_data_l, -6)
                c_data_r = psylab.signal.atten(c_data_r, -6)
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

                # Use 500-us ITD for cue, because it is a tone
                # and other cues don't always lateralize it sufficiently
                #cue_data_l = np.concatenate((cue_data_l, np.zeros(int(500 / 1000000 * fs)) ))
                #cue_data_r = np.concatenate((np.zeros(int(500 / 1000000 * fs)), cue_data_r))

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
                    # We do this because we need itds in order
                    # to estimate gain track for magnification
                    c_f_l = np.convolve(c_data_l,hrtf_l_data)
                    c_f_r = np.convolve(c_data_r,hrtf_r_data)
                    t_f_l = np.convolve(t_data_l,hrtf_l_data)
                    t_f_r = np.convolve(t_data_r,hrtf_r_data)
                    #m1_f_l = np.convolve(m1_data_l,hrtf_l_data)
                    #m1_f_r = np.convolve(m1_data_r,hrtf_r_data)
                    m2_f_l = np.convolve(m2_data_l,hrtf_r_data)
                    m2_f_r = np.convolve(m2_data_r,hrtf_l_data)

                    # Apply magspec from hrtf separately
                    # This will receive magnification
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
                        c_m_bp_l = psylab.signal.filter_bank(c_f_l, fs, 4, ild_bands)
                        c_m_bp_r = psylab.signal.filter_bank(c_f_r, fs, 4, ild_bands)

                        tm_f_bp_l = psylab.signal.filter_bank(tm_f_l, fs, 4, ild_bands)
                        tm_f_bp_r = psylab.signal.filter_bank(tm_f_r, fs, 4, ild_bands)

                        c_bp_g_l = np.zeros(shape=c_f_bp_l.shape)
                        c_bp_g_r = np.zeros(shape=c_f_bp_r.shape)
                        t_m_bp_g_l = np.zeros(shape=t_m_bp_l.shape)
                        t_m_bp_g_r = np.zeros(shape=t_m_bp_r.shape)
                        #m1_m_bp_g_l = np.zeros(shape=m1_m_bp_l.shape)
                        #m1_m_bp_g_r = np.zeros(shape=m1_m_bp_r.shape)
                        m2_m_bp_g_l = np.zeros(shape=m2_m_bp_l.shape)
                        m2_m_bp_g_r = np.zeros(shape=m2_m_bp_r.shape)

                        # Loop through bands, apply magnification
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
                            
                            all_itd_tracks.append(itd_track)
            
                            # Apply gain track not to hrtf'd data, but to mag-spec'd data
                            # This allows magnified ILD cues to be applied without itd cues present
                            # Which is needed for unvocoded stimuli
                            c_bp_g_l[:,i] = ild.apply_gain_track(c_m_bp_l[:,i], -1*c_gain_l)
                            c_bp_g_r[:,i] = ild.apply_gain_track(c_m_bp_r[:,i], -1*c_gain_r)
                            t_m_bp_g_l[:,i] = ild.apply_gain_track(t_m_bp_l[:,i], -1*gain_l)
                            t_m_bp_g_r[:,i] = ild.apply_gain_track(t_m_bp_r[:,i], -1*gain_r)
                            #m1_m_bp_g_l[:,i] = ild.apply_gain_track(m1_m_bp_l[:,i], gain_l)
                            #m1_m_bp_g_r[:,i] = ild.apply_gain_track(m1_m_bp_r[:,i], gain_l)
                            m2_m_bp_g_l[:,i] = ild.apply_gain_track(m2_m_bp_l[:,i], -1*gain_l)
                            m2_m_bp_g_r[:,i] = ild.apply_gain_track(m2_m_bp_r[:,i], -1*gain_r)

                        c_data_l = c_bp_g_l.sum(axis=1)
                        c_data_r = c_bp_g_r.sum(axis=1)
                        t_data_l = t_m_bp_g_l.sum(axis=1)
                        t_data_r = t_m_bp_g_r.sum(axis=1)
                        #m1_data_l = m1_m_bp_g_l.sum(axis=1)
                        #m1_data_r = m1_m_bp_g_r.sum(axis=1)
                        m2_data_l = m2_m_bp_g_l.sum(axis=1)
                        m2_data_r = m2_m_bp_g_r.sum(axis=1)
                from mpl_toolkits.axes_grid1 import make_axes_locatable

                def plot_acoustic_analysis(data, curr_ax, bands, time, plot_min, plot_max, curr_cmap):
                    # Plot snr as a function of time and frequency band
                    fs = 44100
                    snr_map = curr_ax.imshow(data, cmap=curr_cmap, interpolation='nearest', aspect='auto', origin = 'lower', extent = (time[0],time[-1],-0.5, np.size(data,axis=0)-0.5), vmin=plot_min, vmax = plot_max)
                    
                    plt.colorbar(snr_map)
                    divider = make_axes_locatable(curr_ax) 
                    # Add frequency labels
                    freq_labels = []
                    for iband in range(0, len(bands) - 1):
                        freq_labels.append(str(round(bands[iband])) + "-" + str(round(bands[iband+1])) + " Hz")
                        # Add horizontal lines at boundaries
                        curr_ax.axhline(y = iband + 0.5, xmin = time[0], xmax = time[-1], color = 'k',linestyle = '-')
                        
                        # plot data
                        # this_data = data[iband,:]
                        # this_data = (this_data - np.min(this_data))/(np.max(this_data) - np.min(this_data))
                        # curr_ax.plot(this_data + iband*2)
                        
                        
                    curr_ax.set_yticks(range(0, len(bands) - 1))
                    curr_ax.get_yaxis().set_ticklabels(freq_labels)
                fig, axes = plt.subplots(nrows=1,ncols=1)
                curr_ax = axes
                time = np.arange(np.shape(all_itd_tracks)[1])
                import matplotlib
                cmap_rwb = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","white","black"])
                plot_acoustic_analysis(np.array(all_itd_tracks), curr_ax, ild_bands, time, -100, 100, cmap_rwb)

            ## Add location cue; combine target and masker
            #sig_l = np.concatenate((c_data_l, c_isi_data, t_data_l + m2_data_l))
            #sig_r = np.concatenate((c_data_r, c_isi_data, t_data_r + m2_data_r))
            #if levels['side'] == 'r':
            #    sig_out = np.column_stack((sig_r, sig_l))
            #else:
            #    sig_out = np.column_stack((sig_l, sig_r))

            if levels['side'] == 'r':
                sig_l = np.concatenate((c_data_r, c_isi_data, t_data_l + m2_data_l))
                sig_r = np.concatenate((c_data_l, c_isi_data, t_data_r + m2_data_r))
            else:
                sig_l = np.concatenate((c_data_l, c_isi_data, t_data_l + m2_data_l))
                sig_r = np.concatenate((c_data_r, c_isi_data, t_data_r + m2_data_r))
            
            sig_out = np.column_stack((sig_l, sig_r))

            if pract:
                sig_out = psylab.signal.atten(sig_out, pract_atten)

            #sfwrite(this_sf_filepath_out, sig_out, fs, format='wav')
            write(this_sf_filepath_out, sig_out, fs)
            t_str += ','.join([','.join(map(str, i)) for i in zip(t_words, t_times)]) + "\n"
            m_str += ','.join([','.join(map(str, i)) for i in zip(m2_words, m2_times)]) + "\n"

            fid_all_times = open(os.path.join(out_path_data, f'{exp}__s_{subjectID}__Word_Times.csv'), 'a')
            fid_all_times.write(t_str)
            fid_all_times.write(m_str)
            fid_all_times.close()

        if pract:
            pract = False

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
