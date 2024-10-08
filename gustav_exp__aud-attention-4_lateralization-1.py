# -*- coding: utf-8 -*-

# A Gustav settings file!

import os, sys
import numpy as np
import time
import psylab
import gustav
import gamepad
from gustav.forms import lateralization as theForm
import medussa as m

def setup(exp):

    machine = psylab.config.local_settings(conf_file='config/psylab.conf')
    exp.user.machine_name = machine.get_str('name')
    exp.user.audio_id = machine.get_list_int('audio_id')
    workdir = machine.get_path('workdir')
    exp.stim.stimdir = machine.get_path('stimdir')
    dev_name = machine.get_str('audiodev_name')
    dev_ch = machine.get_int('audiodev_ch')
    #dev_name = "Fireface800"
    #dev_ch = 8
    devs = m.get_available_devices()
    dev_id = None
    for i,di in enumerate(devs):
        name = psylab.string.as_str(di.name)
        if name.startswith(dev_name):
            dev_id = i,i,dev_ch
    if dev_id:
        exp.stim.audiodev = m.open_device(*dev_id)
    else:
        raise Exception(f"The audio device {dev_name} was not found")

    exp.stim.gamepad = gamepad.gamepad()
    exp.stim.range_gamepad = machine.get_list_float('gamepad_cal_x')  # Gamepad range
    exp.stim.range_binaural = {'ild': [-30,30],                       # ILD range to use
                               'itd': [-1500,1500],                   # ITD range to use
                              }
    exp.stim.units = {'ild': 'dB',
                      'itd': 'us',
                     }
    exp.stim.locations = ['-30', '30', '-90', '90']
    exp.stim.fs = 44100

    this_syll_1,fs = m.read_file("stim/aud-attention-4_hrtf/syllables/ba_M2.wav")
    this_syll_2,fs = m.read_file("stim/aud-attention-4_hrtf/syllables/da_M2.wav")
    this_syll_3,fs = m.read_file("stim/aud-attention-4_hrtf/syllables/ga_M2.wav")
    exp.stim.isi = np.zeros(int(psylab.signal.ms2samp(250, exp.stim.fs)))
    exp.stim.syllables_clean = [this_syll_1, this_syll_2, this_syll_3]

    # General Experimental Variables
    exp.name = 'aud-attention-4_lat-1'
    exp.method = 'constant' # 'constant' for constant stimuli, or 'adaptive' for a staircase procedure (SRT, etc)
    exp.prompt = 'Match the location of the second word to the first' # Prompt for subject
    exp.frontend = 'tk'
    exp.logFile = os.path.join(workdir,'logs','$name.log')
    exp.logConsole = False
    exp.logConsoleDelay = True
    exp.debug = False
    exp.recordData = True
    exp.dataFile = os.path.join(workdir,'data','$name.csv')
    exp.dataString_header = "# A datafile created by Gustav!\n# \n# Experiment: $name\n# \n\nS,Trial,Date,Block,Condition,@currentvars[],Stim,Resp\n"
    exp.dataString_post_trial = "$subj,$trial,$date,$block,$condition,$currentvars[],$user[az],$response\n"
    exp.logString_pre_block = "\n Block $block of $blocks started at $time; Condition: $condition ; $currentvarsvals[' ; ']\n"
    exp.logString_post_trial = "   Trial $trial, target: $user[az], response: $response\n"
    exp.logString_post_block = " Block $block of $blocks ended at $time; Condition: $condition ; $currentvarsvals[' ; ']\n"
    exp.cacheTrials = False
    exp.validKeys = '1,2';  # Not used for qt_lateralization form, which takes mouseclicks as responses
    exp.quitKey = '/'
    exp.note = "Lateralization task to generate ILDs and ITDs that elicit equivalent lateral positions of hrtfs"
    exp.comments = '''\
    '''

    if not exp.subjID:
        exp.subjID = exp.term.get_input(parent=None, title = "Gustav!", prompt = 'Enter a Subject ID:')

    hrtf_path = os.path.join(workdir,'stim','aud-attention-4_hrtf','hrtfs', f"hrtf_s_{exp.subjID}.npy")
    exp.stim.s_hrtfs = np.load(hrtf_path)
    exp.stim.display = f"\nSubject: {exp.subjID}\n" # To hold lat data for display at end of experiment
    exp.stim.display_dict = {}
    """EXPERIMENT VARIABLES
        There are 2 kinds of variables: factorial and covariable

        Levels added as 'factorial' variables will be factorialized with each
        other. So, if you have 2 fact variables A & B, each with 3 levels, you
        will end up with 9 conditions: A1B1, A1B2, A1B3, A2B1 etc..

        Levels added as 'covariable' variables will simply be listed (in parallel
        with the corresponding levels from the other variables) in the order
        specified. So, if you have 2 'covariable' variables A & B, each with 3
        levels, you will end up with 3 conditions: A1B1, A2B2, and A3B3. All
        'covariable' variables must have either the same number of levels, or
        exactly one level. When only one level is specified, that level will
        be used in all 'covariable' conditions. Eg., A1B1, A2B1, A3B1, etc.

        You can use both types of variables in the same experiment, but both
        factorial and covariable must contain exactly the same set of variable
        names. factorial levels are processed first, covariable levels are added
        at the end.

        Both factorial and covariable are Python ordered dicts, where the keys 
        are variable names, and the values are lists of levels. During the 
        experiment, you have access to the current level of each variable. For 
        example, if you have the following variable:
        
        exp.var.factorial['target'] = ['Male', 'Female']
        
        Then, you can find out what the level is at any point in the experiment 
        with exp.var.current['target'], which would return either 'Male' or 
        'Female' depending on what the condition happened to be. This is 
        probably most useful to generate your stimuli, eg., in the pre_trial 
        function. 
    """
    exp.var.factorial['cue'] =  [
                                        'itd',
                                        'ild',
                                     ]
    exp.var.factorial['location'] = [
                                        '-30',
                                        '30',
                                        '-90',
                                        '90',
                                     ]

    """CONSTANT METHOD VARIABLES
        The method of constant limits requires three variables to be set.
            trialsperblock
            startblock [crash recovery]
            starttrial [crash recovery]
    """
    exp.var.constant = {
        'trialsperblock' : 1,  
        'startblock' : 1,
        'starttrial' : 1,
        }

    """CONDITION PRESENTATION ORDER
        Use 'prompt' to prompt for condition on each block, 'random' to randomize
        condition order, 'menu' to be able to choose from a list of conditions at
        the start of a run, 'natural' to use natural order (1:end), or a
        print-range style string to specify the order ('1-10, 12, 15'). You can
        make the first item in the print range 'random' to randomize the specified
        range.
    """
    exp.var.order = 'x3, random'

    """IGNORE CONDITIONS
        A list of condition numbers to ignore. These conditions will not be
        reflected in the total number of conditions to be run, etc. They will
        simply be skipped as they are encountered during a session.
    """
    exp.var.ignore = []

    '''USER VARIABLES
        Add any additional variables you need here
    '''

"""CUSTOM PROMPT
    If you want a custom response prompt, define a function for it
    here. exp.run.response should receive the response as a string, and
    if you want to cancel the experiment, set both exp.run.block_on and
    exp.run.pylab_is_go to False
"""
def prompt_response(exp):
    pass
    #resp = None
    #while True:
    #    ax,pos = exp.stim.gamepad.listen(timeout=.2)
    #    if ax=='x':
    #        resp = pos

    #        ((pos - exp.stim.range_gamepad[0]) / (exp.stim.range_gamepad[1]-exp.stim.range_gamepad[0]) * ((exp.stim.range_bin[1]-exp.stim.range_bin[0])) + exp.stim.range_[0])
    #    elif ax=='trigger' and pos==0 and resp:
    #        exp.run.response = resp
    #        break

        #r = exp.interface.get_resp(exp.prompt)
        #if r == exp.quitKey:
        #    exp.run.block_on = False
        #    exp.run.gustav_is_go = False
        #    break
        #elif isinstance(r, (float,int)) and int(r)*5 in exp.stim.locations:
        #    exp.run.response = str(int(r)*5)
        #    break

"""PRE_TRIAL
    This function gets called on every trial to generate the stimulus, and
    do any other processing you need. All settings and variables are
    available. For the current level of a variable, use
    exp.var.current['varname']. 
"""
def pre_trial(exp):
    try:
        responding = True # Set to False when user presses button
        i = 0
        while responding:

            stan_l = exp.stim.stan_syllables_l[exp.stim.stan_order[i]]
            stan_r = exp.stim.stan_syllables_r[exp.stim.stan_order[i]]

            compar = exp.stim.syllables_clean[exp.stim.comp_order[i]]
     
            if exp.var.current['cue'] == 'itd':
                delay = np.zeros(psylab.signal.ms2samp(np.abs(exp.stim.this_cue)/1000, exp.stim.fs))
                if exp.stim.this_cue < 0:
                    comp_l = np.concatenate((compar, delay))
                    comp_r = np.concatenate((delay, compar))
                else:
                    comp_l = np.concatenate((delay, compar))
                    comp_r = np.concatenate((compar, delay))
            else:
                if exp.stim.this_cue > 0:
                    comp_l = psylab.signal.atten(compar, np.abs(exp.stim.this_cue))
                    comp_r = compar
                else:
                    comp_l = compar
                    comp_r = psylab.signal.atten(compar, np.abs(exp.stim.this_cue))
            
            exp.interface.update_Title_Center(f"{exp.stim.this_cue} {exp.stim.unit}", redraw=True)
            sig_l = psylab.signal.atten(np.concatenate((stan_l, exp.stim.isi, comp_l, exp.stim.isi, exp.stim.isi)), 10)
            sig_r = psylab.signal.atten(np.concatenate((stan_r, exp.stim.isi, comp_r, exp.stim.isi, exp.stim.isi)), 10)
            sig_out = np.column_stack((sig_l, sig_r))

            s = exp.stim.audiodev.open_array(sig_out,exp.stim.fs)
            s.play()
            exp.interface.show_Notify_Right(show=True, redraw=True)
            while s.is_playing:
                ax,pos = exp.stim.gamepad.listen(timeout=.1)
                if ax=='x':
                    # convert axis position (-1>=1) to binaural cue (ild: -20>=20; itd: -1000>=1000)
                    # range_gamepad is the calibrated physical limits of the joystick (eg, move all the way right may not return 1.0)
                    exp.stim.pos_0_1 = pos
                    exp.stim.this_cue = exp.stim.pos_0_1 * (exp.stim.range_bin[1]-exp.stim.range_bin[0]) + exp.stim.range_bin[0]
                    exp.interface.set_marker_pos(pos, show=True)

                elif ax=='trigger' and pos==0:
                    exp.stim.display += f"{exp.var.current['location']}: {np.round(exp.stim.this_cue,2)} {exp.stim.unit}\n"
#                    if exp.var.current['location'] in exp.stim.display_dict.keys():
#                        exp.stim.display_dict[exp.var.current['location']].append(np.round(exp.stim.this_cue,2))

                    responding = False
                    break
            exp.interface.show_Notify_Right(show=False, redraw=True)
            i += 1
            if i == 30:
                i = 0
    except Exception as e:
        exp.interface.destroy()
        raise e


def present_trial(exp):
    pass
    #m.play_array(exp.stim.out,exp.stim.fs)
    #s = exp.audiodev.open_array(exp.stim.out,exp.stim.fs)
    #s.play()
    #while s.is_playing:
    #    time.sleep(0.1)

def post_trial(exp):
    pass

def pre_exp(exp):
    exp.interface = theForm.Interface()
    exp.interface.update_Notify_Left("Hit a key\nto begin", show=False, redraw=False)
    exp.interface.update_Notify_Right('Playing', show=False, redraw=False)
    exp.interface.update()

def post_exp(exp):
    exp.interface.destroy()
    print(exp.stim.display)

def pre_block(exp):
    try:
        exp.stim.stan_order = []
        exp.stim.comp_order = []
        this_ord = np.arange(3)
        for i in np.arange(10):
            np.random.shuffle(this_ord)
            exp.stim.stan_order.extend(this_ord)
            np.random.shuffle(this_ord)
            exp.stim.comp_order.extend(this_ord)

        exp.stim.range_bin = exp.stim.range_binaural[exp.var.current['cue']] # use the right range (itd or ild) for this block
        exp.stim.unit = exp.stim.units[exp.var.current['cue']]
        exp.stim.this_cue = 0 # Start cue at zero
        loc_i = exp.stim.locations.index(exp.var.current['location']) # Deg az -> index
        hrtf_l_data = exp.stim.s_hrtfs[:, 0, loc_i]                   # Grab left ear for current loc
        hrtf_r_data = exp.stim.s_hrtfs[:, 1, loc_i]                   # Grab right ear for current loc

        ba_l = np.convolve(exp.stim.syllables_clean[0],hrtf_l_data)
        ba_r = np.convolve(exp.stim.syllables_clean[0],hrtf_r_data)
        da_l = np.convolve(exp.stim.syllables_clean[1],hrtf_l_data)
        da_r = np.convolve(exp.stim.syllables_clean[1],hrtf_r_data)
        ga_l = np.convolve(exp.stim.syllables_clean[2],hrtf_l_data)
        ga_r = np.convolve(exp.stim.syllables_clean[2],hrtf_r_data)

        # Equate to unproc. But since there are ILDs, use the average of the spatialized versions
        ba_l *= psylab.signal.rms(exp.stim.syllables_clean[0]) / np.mean((psylab.signal.rms(ba_l), psylab.signal.rms(ba_r)))
        ba_r *= psylab.signal.rms(exp.stim.syllables_clean[0]) / np.mean((psylab.signal.rms(ba_l), psylab.signal.rms(ba_r)))
        da_l *= psylab.signal.rms(exp.stim.syllables_clean[1]) / np.mean((psylab.signal.rms(da_l), psylab.signal.rms(da_r)))
        da_r *= psylab.signal.rms(exp.stim.syllables_clean[1]) / np.mean((psylab.signal.rms(da_l), psylab.signal.rms(da_r)))
        ga_l *= psylab.signal.rms(exp.stim.syllables_clean[2]) / np.mean((psylab.signal.rms(ga_l), psylab.signal.rms(ga_r)))
        ga_r *= psylab.signal.rms(exp.stim.syllables_clean[2]) / np.mean((psylab.signal.rms(ga_l), psylab.signal.rms(ga_r)))

        # Make left-ear and right-ear lists of each spatialized syllable
        exp.stim.stan_syllables_l = [ba_l, da_l, ga_l]
        exp.stim.stan_syllables_r = [ba_r, da_r, ga_r]

        exp.interface.update_Status_Right("%g of %g" % (exp.run.block+1, exp.run.nblocks), redraw=False)
        exp.interface.update_Status_Center("Cue: %s" % exp.var.current['cue'], redraw=False)
        exp.interface.show_Notify_Left(show=True, redraw=True)
        while True:
            ax,pos = exp.stim.gamepad.listen(timeout=.1)
            if ax=='trigger' and pos==0:
                break
        exp.interface.show_Notify_Left(show=False, redraw=True)


        #exp.interface.set_text_block("Block %g of %g" % (exp.run.block+1, exp.run.nblocks))
        #exp.interface.update()
    except Exception as e:
        exp.interface.destroy()
        raise e


if __name__ == '__main__':
    argv = sys.argv[1:]
    argv.append("--experimentFile=%s" % os.path.realpath(__file__))
    gustav.main(argv)
#    import inspect
#    fname = inspect.getfile( inspect.currentframe() )
#    psylab.gustav.run(settingsFile=fname)
