%% mildmaster_button_press_postprocessing.m
% Benjamin Richardson


% Set directories
whos_using = 'Bon';
if all(whos_using == 'Ben')
    addpath('/home/ben/Documents/MATLAB/eeglab2023.1');
    dir = '/home/ben/Documents/GitHub/MILD-Master/';
    dir_mildmaster = '/home/ben/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER Behavior Files/button-press-pilot.xlsx';
elseif all(whos_using == 'Bon')
    addpath('C:\Users\benri\Documents\eeglab2023.1');
    dir = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
    dir_mildmaster = 'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\button-press-pilot.xlsx';
    prepro_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';
end

curr_subject_ID = char('button_press_pilot_1','button_press_pilot_2'); % char();
% Set analysis parameters
erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 750; % 750 ms after onset of word
nsubjects = size(curr_subject_ID,1);
word_length = 0.3;
frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
fs = 256;

for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID)
    data_by_button_press = [];

    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == strtrim(string(curr_subject_ID(isubject,:)))); % find the rows in the spreadsheet which belong to this subject
    conditions = BehaviorTable.Condition(rows_this_subject); % conditions by trial for this subject
    condition_names = {'start_press'};
    this_subject_table = BehaviorTable(rows_this_subject,:);


    epochs_filename = join([prepro_folder,strtrim(curr_subject_ID(isubject,:)),'all_epoch.mat'],'');
    this_EEG = load(epochs_filename);
    eeg_struct_name = fieldnames(this_EEG);
    this_EEG = getfield(this_EEG,string(eeg_struct_name(1)));
    these_epochs = this_EEG.data; % 32 channels x Time x 240 trials
    num_tot_trials = size(these_epochs,3); % look into this

    ERP_info_button_press = struct('SubID',{},'Trial',{});

    eeg_time = this_EEG.times; % in milliseconds

    noise_thresh = 50; % 80;

    for itrial = 1:size(this_EEG.data,3)% for each trial (should be 144)
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        if mod(itrial,10) == 0
            disp(itrial)
        end
        icondition = conditions(itrial);

        which_soundfile_this_trial = BehaviorTable.Soundfile(rows_this_subject(itrial));
        which_soundfile_this_trial = cell2mat(which_soundfile_this_trial);
        slash_indices = find(which_soundfile_this_trial == '/');
        slash_index = max(slash_indices);
        which_soundfile_this_trial = which_soundfile_this_trial(slash_index + 1:slash_index + 3);
        if contains(which_soundfile_this_trial, '__')
            which_soundfile_this_trial = str2num(which_soundfile_this_trial(1));
        elseif contains(which_soundfile_this_trial,'_')
            which_soundfile_this_trial = str2num(which_soundfile_this_trial(1:2));
        else
            which_soundfile_this_trial = str2num(which_soundfile_this_trial);
        end

        %% ISOLATE BUTTON PRESSES
        % Find this trial button presses
        this_trial_click_times = table2array(this_subject_table(itrial,9:end));
        this_trial_click_times(isnan(this_trial_click_times)) = [];
        for iclick = 1:length(this_trial_click_times) % for each target word onset...
            resampled_search_time = this_trial_click_times(iclick);
            button_press_delay = 0; % ms
            [~,start_time] = min(abs(eeg_time - (resampled_search_time + erp_window_start_time + button_press_delay))); %
            [~,end_time] = min(abs(eeg_time - (resampled_search_time + erp_window_end_time)));%


            if end_time - start_time == 204
                end_time = end_time + 1;
            end

            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(these_epochs(:,start_time:end_time,itrial)) > noise_thresh,'all')
                disp('ERP rejected')
                continue


            end

            % Isolate ERP
            this_erp = these_epochs(:,start_time:end_time,itrial);
            data_by_button_press = cat(3, data_by_button_press,this_erp);

            % Append Info
            if ~isempty(ERP_info_button_press)
                ERP_info_button_press.SubID = [ERP_info_button_press.SubID; curr_subject_ID(isubject,:)];
                ERP_info_button_press.Trial = [ERP_info_button_press.Trial, itrial];
            else
                ERP_info_button_press(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_button_press(1).Trial = itrial;
            end



        end

    end

    data_by_button_press_baselined = nan(size(data_by_button_press));

    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press,2));
    [~,baseline_start_index_buttonpress] = min(abs(single_onset_time_buttonpress - erp_window_start_time));
    [~,baseline_end_index_buttonpress] = min(abs(single_onset_time_buttonpress - 0));

    for ichannel = 1:32
        % all data
        data_by_button_press_baselined(ichannel,:,:) = data_by_button_press(ichannel,:,:) - mean(data_by_button_press(ichannel,baseline_start_index_buttonpress:baseline_end_index_buttonpress,:),'all');
    end

    %% SAVE INFO FOR THIS SUBBY
    save(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'),'data_by_button_press_baselined','ERP_info_button_press','-v7.3')

end
