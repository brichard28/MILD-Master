## vocode_stimuli.py

# Created by: Benjamin Richardson
# May 16th, 2026

# This code imports stimuli for Sahil and Abby's grant proposal and vocodes them.

# FROM THE PROPOSAL: 8-channel noise-vocoder, where the eight channels comprise a set of band-pass filters
# (cutoff frequencies logarithmically spaced between 125 and 5500 Hz). Within each channel,
# the amplitude envelope will be extracted using half-wave rectification, low-pass filtered
# (at either 400 Hz or half of the channel’s bandwidth, whichever is lesser),
# and used to modulate band-passed noise; the output will be summed across the eight channels.


import soundfile as sf

# Load the audio data and sampling rate
data, samplerate = sf.read('audio.flac')

