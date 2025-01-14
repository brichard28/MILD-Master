% Primary Authors: Victoria Figarola, Benjamin Richardson 7/21/23
% Secondary Authors: Emaya Anand, Maanasa Guru Adimurthy
% EPOCHING

% 'fullpilot1','fullpilot2','fullpilot3','eeg_pilot_1','eeg_pilot_4','eeg_pilot_5','eeg_pilot_6','eeg_pilot_7','eeg_pilot_8'); % char();
subID = ['button_press_pilot_2']; % set current subject ID

% Set directories
whos_using = 'Bon';

if whos_using == 'Ben'
    addpath('/home/ben/Documents/MATLAB/eeglab2023.1/')
    pre_pro_epoched_data_folder = '/home/ben/Documents/GitHub/fNIRSandGerbils/prepro_epoched_data/';
elseif whos_using == 'Bon' % Ben Laptop
    addpath('\Users\benri\Documents\eeglab2023.0\')
    pre_pro_epoched_data_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';

end

% Load dataset
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; 
EEG = pop_loadset('filename', [subID, '_ICAdone.set'], 'filepath', pre_pro_epoched_data_folder);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
EEG = eeg_checkset( EEG );


% remove extraneous triggers

EEG.event(~ismember(string({EEG.event(:).type}), {'31231','63999'})) = [];
EEG.urevent(~ismember([EEG.urevent(:).type],[31231,63999])) = [];



%check trigger latency distances, remove double triggers
distance_threshold = 500;
all_latencies = [EEG.urevent(:).latency];
all_types = [EEG.urevent(:).type];
all_distances = diff(all_latencies);
num_dist_below_threshold = sum(all_distances < distance_threshold);
disp('Below is the number of instances where triggers are too close together');
disp(num_dist_below_threshold);
figure; xline(all_latencies);

i = 1;
urevents_to_remove = [];
while i < numel(all_latencies)
    if all_latencies(i+1) - all_latencies(i) < distance_threshold
        urevents_to_remove = [urevents_to_remove, i+1];
    end
    i = i+1;
end
where_below_threshold = find(all_distances < distance_threshold);
% find corresponding events to urevents and remove triggers from EEG.event
% and EEG.urevent
z = {EEG.event(:).urevent};
z(find(cellfun(@isempty,z))) = {0};
events_to_remove = [];
for j = 1:length(urevents_to_remove)
    this_urevent = urevents_to_remove(j);
    events_to_remove = [events_to_remove, find(string(z) == string(this_urevent))];
end
EEG.event(events_to_remove) = [];
EEG.urevent(urevents_to_remove) = [];

%all epochs

EEG = pop_epoch( EEG, {"31231" , "63999"}, [0  16], 'newname', [subID, 'all epochs'], 'epochinfo', 'yes');

EEG = eeg_checkset( EEG );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'gui', 'off');
EEG = eeg_checkset( EEG );
save([pre_pro_epoched_data_folder ,subID, 'all_epoch.mat'], "EEG")

eeglab redraw;