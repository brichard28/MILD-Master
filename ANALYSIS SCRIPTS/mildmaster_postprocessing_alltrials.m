%% mildmaster_postprocessing_alltrials.m

% Set directories
whos_using = 'Bon';
if all(whos_using == 'Ben')
    addpath('/home/ben/Documents/MATLAB/eeglab2023.1');
    dir = '/home/ben/Documents/GitHub/fNIRSandGerbils/';
    dir_mildmaster = '/home/ben/Documents/GitHub/fNIRSandGerbils/data/fNIRSandGerbils.xlsx';
elseif all(whos_using == 'Bon')
    addpath('C:\Users\benri\Documents\eeglab2023.1');
    dir = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
    dir_mildmaster = 'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx';
    prepro_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';
end

% 'mild_master_1',...
%     'mild_master_3',...
%     'mild_master_4',...
%     'mild_master_5',...
%     'mild_master_6',...
%     'mild_master_8',...
%     'mild_master_9',...
%     'mild_master_10',...
%     'mild_master_11',...
%     'mild_master_12',...
%     'mild_master_14',...
%     'mild_master_15',...
%     'mild_master_16',...
%     'mild_master_17',...
%     'mild_master_18',...
%     'mild_master_19',...
%     'mild_master_20',...
%     'mild_master_22',...
%     'mild_master_23',...
%     'mild_master_24',...

curr_subject_ID = char('mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40'); % char();

% Set analysis parameters
erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 950; % 750 ms after onset of word
nsubjects = size(curr_subject_ID,1);
word_length = 0.3;
frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
fs = 256;

mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master';


%% RUN ANALYSIS WITHOUT BUTTON PRESS SUBTRACTION
%% For each subject.....
for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID + " no button press")

    cue_dur = 2.0;

    % Load word onset times data
    WordTimesTable = readtable(append(mild_master_root,"\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(subID) + "__Word_Times.csv"));


    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == strtrim(string(curr_subject_ID(isubject,:)))); % find the rows in the spreadsheet which belong to this subject
    conditions = BehaviorTable.Condition(rows_this_subject); % conditions by trial for this subject
    condition_names = {'side=l_itd=0_az_itd=0_az=5_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=0_az=15_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=5_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=0_az=5_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=0_az=15_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=15_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=5_az=0_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=15_az=0_mag=0_lpf=0'};
    small_itd_cond = [3,7];
    large_itd_cond = [6,8];
    small_ild_cond = [1,4];
    large_ild_cond = [2,5];
    this_subject_table = BehaviorTable(rows_this_subject,:);

    % Create empty arrays for ERPs
    data_by_pair_onset = [];
    data_by_button_press = [];


    % Create empty arrays for info for each ERP
    % Will contain subID, trial, and word (if target)
    ERP_info = struct('SubID',{},'Trial',{},'Lead_Stream',{},'Lag_Stream',{},'Lead_Word',{},'Lag_Word',{},'Condition',{});
    ERP_info_button_press = struct('SubID',{},'Trial',{});

    % Load EEG for this subject
    epochs_filename = join([prepro_folder,strtrim(curr_subject_ID(isubject,:)),'all_epoch_no_button_press.mat'],'');
    this_EEG = load(epochs_filename);
    eeg_struct_name = fieldnames(this_EEG);
    this_EEG = getfield(this_EEG,string(eeg_struct_name(1)));
    these_epochs = this_EEG.data; % 32 channels x Time x 240 trials
    num_tot_trials = size(these_epochs,3); % look into this


    num_trials = size(this_EEG.data,3);
    if subID == "mild_master_2"
        num_trials = 78;
    elseif ismember(subID,["mild_master_9 ","mild_master_10","mild_master_16","mild_master_18","mild_master_28"])
        num_trials = 119;
    end
    % baseline all epochs to 1 second



    % Define time vector for extracting target ERPs
    eeg_time = this_EEG.times; % in milliseconds
    resampled_audio_time = -1:1/fs:16;
    resampled_audio_time = resampled_audio_time.*1000;


    %     [~,baseline_start_index] = min(abs(resampled_audio_time - -1));
    %     [~,baseline_end_index] = min(abs(resampled_audio_time - 0));
    %     for ichannel = 1:32
    %         these_epochs(ichannel,:,:) = these_epochs(ichannel,:,:) - mean(these_epochs(ichannel,baseline_start_index:baseline_end_index,:),'all');
    %     end

    % Define time vector for extracting masker ERPs

    noise_thresh = 50; % 80;
    
    for itrial = 1:num_trials% for each trial (should be 120)
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        if mod(itrial,20) == 0
            disp(itrial)
        end
        icondition = conditions(itrial);

        %this_trial_target_onsets = all_target_onsets(itrial).onsets;
        % actually this_trial_target_onsets should be the one that matches
        % the SOUNDFILE, not the trial number. They were not played in that
        % order
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

        this_trial_target_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = (table2array(this_trial_target_all(:,2:2:end))./44100) + cue_dur; % have to add cue back in

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = (table2array(this_trial_masker_all(:,2:2:end))./44100) + cue_dur;

        % Find which token leads (to sort ERPs that way as well)
        this_trial_whether_target_lead = sign(this_trial_masker_times - this_trial_target_times);
        this_trial_whether_target_lead(this_trial_whether_target_lead == -1) = 0; % 0 when masker leads, 1 when target leads

        this_trial_target_bash_times = this_trial_target_times(ismember(this_trial_target_words,{'bash'}));
        this_trial_masker_bash_times = this_trial_masker_times(ismember(this_trial_masker_words,{'bash'}));
        %% ISOLATE BUTTON PRESSES
        % Find this trial button presses
        this_trial_click_times = table2array(this_subject_table(itrial,9:end));
        this_trial_click_times(isnan(this_trial_click_times)) = [];
        for iclick = 1:length(this_trial_click_times) % for each target word onset...
            this_click_time = this_trial_click_times(iclick) + 0.702*fs;
            resampled_search_time = this_click_time;
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

        %% ISOLATE WORD PAIR ONSETS

        this_trial_pair_onset_times = 2:0.6:11;
        for ionset = 2:length(this_trial_pair_onset_times) % for each word pair onset, excluding the first
            if this_trial_whether_target_lead(ionset) == 1
                resampled_search_time = this_trial_pair_onset_times(ionset)*1000;
            elseif this_trial_whether_target_lead(ionset) == 0
                resampled_search_time = this_trial_pair_onset_times(ionset)*1000;
            end

            [~,start_time] = min(abs(eeg_time - (resampled_search_time + erp_window_start_time))); %
            [~,end_time] = min(abs(eeg_time - (resampled_search_time + erp_window_end_time)));%

            if end_time - start_time == 204
                end_time = end_time + 1;
            end


            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(these_epochs(:,start_time:end_time,itrial)) > noise_thresh,'all')
                %disp('ERP rejected')
                continue
            end

            % Isolate ERP

            this_erp = these_epochs(:,start_time:end_time,itrial);
            data_by_pair_onset = cat(3,data_by_pair_onset,this_erp);

            % add to lead or lag ERPs, depending


            % Append Info
            if ~isempty(ERP_info)
                ERP_info.SubID = [ERP_info.SubID; curr_subject_ID(isubject,:)];
                ERP_info.Trial = [ERP_info.Trial, itrial];
                if this_trial_whether_target_lead(ionset) == 1
                    ERP_info.Lead_Stream = [ERP_info.Lead_Stream; 'Target'];
                    ERP_info.Lag_Stream = [ERP_info.Lag_Stream; 'Masker'];
                    ERP_info.Lead_Word = [ERP_info.Lead_Word; this_trial_target_words(ionset)];
                    ERP_info.Lag_Word = [ERP_info.Lag_Word; this_trial_masker_words(ionset)];
                elseif this_trial_whether_target_lead(ionset) == 0
                    ERP_info.Lead_Stream = [ERP_info.Lead_Stream; 'Masker'];
                    ERP_info.Lag_Stream = [ERP_info.Lag_Stream; 'Target'];
                    ERP_info.Lead_Word = [ERP_info.Lead_Word; this_trial_masker_words(ionset)];
                    ERP_info.Lag_Word = [ERP_info.Lag_Word; this_trial_target_words(ionset)];
                end
                ERP_info.Condition = [ERP_info.Condition, conditions(itrial)];
            else
                ERP_info(1).SubID = curr_subject_ID(isubject,:);
                ERP_info(1).Trial = itrial;
                if this_trial_whether_target_lead(ionset) == 1
                    ERP_info(1).Lead_Stream = 'Target';
                    ERP_info(1).Lag_Stream = 'Masker';
                    ERP_info(1).Lead_Word = this_trial_target_words(ionset);
                    ERP_info(1).Lag_Word = this_trial_masker_words(ionset);
                elseif this_trial_whether_target_lead(ionset) == 0
                    ERP_info(1).Lead_Stream = 'Masker';
                    ERP_info(1).Lag_Stream = 'Target';
                    ERP_info(1).Lead_Word = this_trial_masker_words(ionset);
                    ERP_info(1).Lag_Word = this_trial_target_words(ionset);
                end
                ERP_info(1).Condition = conditions(itrial);
            end


        end



    end

    %% Concatenate and baseline within each channel for this subject
    % Baseline to the mean voltage during the baseline period over ALL
    % trials
    % Masker and target will be baselined TOGETHER
    % Lead and Lag will be baselined separately, for those analyses

    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset,2));
    [~,baseline_start_index] = min(abs(single_onset_time - erp_window_start_time));
    [~,baseline_end_index] = min(abs(single_onset_time - 0));

    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press,2));
    [~,baseline_start_index_buttonpress] = min(abs(single_onset_time_buttonpress - erp_window_start_time));
    [~,baseline_end_index_buttonpress] = min(abs(single_onset_time_buttonpress - 0));

    % all data
    data_by_button_press_baselined = nan(size(data_by_button_press));
    data_by_pair_onset_baselined = nan(size(data_by_pair_onset));

    for ichannel = 1:32
        % all data
        data_by_button_press_baselined(ichannel,:,:) = data_by_button_press(ichannel,:,:) - mean(data_by_button_press(ichannel,baseline_start_index_buttonpress:baseline_end_index_buttonpress,:),'all');
        data_by_pair_onset_baselined(ichannel,:,:) = data_by_pair_onset(ichannel,:,:) - mean(data_by_pair_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');

    end

    all_data_button_press(isubject,:,:) = squeeze(mean(data_by_button_press_baselined,3));
    all_data_pair_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined,3));

    %all_info_button_press(isubject).info = ERP_info_button_press;
    %all_info_target(isubject).info = ERP_info_target;

    % sort all data into conditions
    itd5_by_pair_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond))),3));

    itd15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond))),3));

    ild5_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond))),3));

    ild15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond))),3));

    %     colorbar


    %% SAVE INFO FOR THIS SUBBY
    save(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_no_button_press.mat'),'data_by_pair_onset_baselined','data_by_button_press_baselined','ERP_info_button_press','ERP_info','-v7.3')

end











%% RUN ANALYSIS WITH BUTTON PRESS SUBTRACTION
%% For each subject.....
for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID + " WITH button press")

    cue_dur = 2.0;

    % Load word onset times data
    WordTimesTable = readtable(append(mild_master_root,"\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(subID) + "__Word_Times.csv"));


    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == strtrim(string(curr_subject_ID(isubject,:)))); % find the rows in the spreadsheet which belong to this subject
    conditions = BehaviorTable.Condition(rows_this_subject); % conditions by trial for this subject
    condition_names = {'side=l_itd=0_az_itd=0_az=5_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=0_az=15_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=5_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=0_az=5_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=0_az=15_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=15_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az_itd=5_az=0_mag=0_lpf=0',...
        'side=l_itd=0_az_itd=15_az=0_mag=0_lpf=0'};
    small_itd_cond = [3,7];
    large_itd_cond = [6,8];
    small_ild_cond = [1,4];
    large_ild_cond = [2,5];
    this_subject_table = BehaviorTable(rows_this_subject,:);

    % Create empty arrays for ERPs
    data_by_pair_onset = [];
    data_by_button_press = [];


    % Create empty arrays for info for each ERP
    % Will contain subID, trial, and word (if target)
    ERP_info = struct('SubID',{},'Trial',{},'Lead_Stream',{},'Lag_Stream',{},'Lead_Word',{},'Lag_Word',{},'Condition',{});
    ERP_info_button_press = struct('SubID',{},'Trial',{});

    % Load EEG for this subject
    epochs_filename = join([prepro_folder,strtrim(curr_subject_ID(isubject,:)),'all_epoch_yes_button_press.mat'],'');
    this_EEG = load(epochs_filename);
    eeg_struct_name = fieldnames(this_EEG);
    this_EEG = getfield(this_EEG,string(eeg_struct_name(1)));
    these_epochs = this_EEG.data; % 32 channels x Time x 240 trials
    num_tot_trials = size(these_epochs,3); % look into this


    num_trials = size(this_EEG.data,3);
    if subID == "mild_master_2"
        num_trials = 78;
    elseif ismember(subID,["mild_master_9 ","mild_master_10","mild_master_16","mild_master_18","mild_master_28"])
        num_trials = 119;
    end
    % baseline all epochs to 1 second



    % Define time vector for extracting target ERPs
    eeg_time = this_EEG.times; % in milliseconds
    resampled_audio_time = -1:1/fs:16;
    resampled_audio_time = resampled_audio_time.*1000;


    %     [~,baseline_start_index] = min(abs(resampled_audio_time - -1));
    %     [~,baseline_end_index] = min(abs(resampled_audio_time - 0));
    %     for ichannel = 1:32
    %         these_epochs(ichannel,:,:) = these_epochs(ichannel,:,:) - mean(these_epochs(ichannel,baseline_start_index:baseline_end_index,:),'all');
    %     end

    % Define time vector for extracting masker ERPs

    noise_thresh = 50; % 80;
    
    for itrial = 1:num_trials% for each trial (should be 120)
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        if mod(itrial,20) == 0
            disp(itrial)
        end
        icondition = conditions(itrial);

        %this_trial_target_onsets = all_target_onsets(itrial).onsets;
        % actually this_trial_target_onsets should be the one that matches
        % the SOUNDFILE, not the trial number. They were not played in that
        % order
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

        this_trial_target_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = (table2array(this_trial_target_all(:,2:2:end))./44100) + cue_dur; % have to add cue back in

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = (table2array(this_trial_masker_all(:,2:2:end))./44100) + cue_dur;

        % Find which token leads (to sort ERPs that way as well)
        this_trial_whether_target_lead = sign(this_trial_masker_times - this_trial_target_times);
        this_trial_whether_target_lead(this_trial_whether_target_lead == -1) = 0; % 0 when masker leads, 1 when target leads

        this_trial_target_bash_times = this_trial_target_times(ismember(this_trial_target_words,{'bash'}));
        this_trial_masker_bash_times = this_trial_masker_times(ismember(this_trial_masker_words,{'bash'}));
        %% ISOLATE BUTTON PRESSES
        % Find this trial button presses
        this_trial_click_times = table2array(this_subject_table(itrial,9:end));
        this_trial_click_times(isnan(this_trial_click_times)) = [];
        for iclick = 1:length(this_trial_click_times) % for each target word onset...
            this_click_time = this_trial_click_times(iclick) + 0.702*fs;
            resampled_search_time = this_click_time;
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

        %% ISOLATE WORD PAIR ONSETS

        this_trial_pair_onset_times = 2:0.6:11;
        for ionset = 2:length(this_trial_pair_onset_times) % for each word pair onset, excluding the first
            if this_trial_whether_target_lead(ionset) == 1
                resampled_search_time = this_trial_pair_onset_times(ionset)*1000;
            elseif this_trial_whether_target_lead(ionset) == 0
                resampled_search_time = this_trial_pair_onset_times(ionset)*1000;
            end

            [~,start_time] = min(abs(eeg_time - (resampled_search_time + erp_window_start_time))); %
            [~,end_time] = min(abs(eeg_time - (resampled_search_time + erp_window_end_time)));%

            if end_time - start_time == 204
                end_time = end_time + 1;
            end


            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(these_epochs(:,start_time:end_time,itrial)) > noise_thresh,'all')
                %disp('ERP rejected')
                continue
            end

            % Isolate ERP

            this_erp = these_epochs(:,start_time:end_time,itrial);
            data_by_pair_onset = cat(3,data_by_pair_onset,this_erp);

            % add to lead or lag ERPs, depending


            % Append Info
            if ~isempty(ERP_info)
                ERP_info.SubID = [ERP_info.SubID; curr_subject_ID(isubject,:)];
                ERP_info.Trial = [ERP_info.Trial, itrial];
                if this_trial_whether_target_lead(ionset) == 1
                    ERP_info.Lead_Stream = [ERP_info.Lead_Stream; 'Target'];
                    ERP_info.Lag_Stream = [ERP_info.Lag_Stream; 'Masker'];
                    ERP_info.Lead_Word = [ERP_info.Lead_Word; this_trial_target_words(ionset)];
                    ERP_info.Lag_Word = [ERP_info.Lag_Word; this_trial_masker_words(ionset)];
                elseif this_trial_whether_target_lead(ionset) == 0
                    ERP_info.Lead_Stream = [ERP_info.Lead_Stream; 'Masker'];
                    ERP_info.Lag_Stream = [ERP_info.Lag_Stream; 'Target'];
                    ERP_info.Lead_Word = [ERP_info.Lead_Word; this_trial_masker_words(ionset)];
                    ERP_info.Lag_Word = [ERP_info.Lag_Word; this_trial_target_words(ionset)];
                end
                ERP_info.Condition = [ERP_info.Condition, conditions(itrial)];
            else
                ERP_info(1).SubID = curr_subject_ID(isubject,:);
                ERP_info(1).Trial = itrial;
                if this_trial_whether_target_lead(ionset) == 1
                    ERP_info(1).Lead_Stream = 'Target';
                    ERP_info(1).Lag_Stream = 'Masker';
                    ERP_info(1).Lead_Word = this_trial_target_words(ionset);
                    ERP_info(1).Lag_Word = this_trial_masker_words(ionset);
                elseif this_trial_whether_target_lead(ionset) == 0
                    ERP_info(1).Lead_Stream = 'Masker';
                    ERP_info(1).Lag_Stream = 'Target';
                    ERP_info(1).Lead_Word = this_trial_masker_words(ionset);
                    ERP_info(1).Lag_Word = this_trial_target_words(ionset);
                end
                ERP_info(1).Condition = conditions(itrial);
            end


        end



    end

    %% Concatenate and baseline within each channel for this subject
    % Baseline to the mean voltage during the baseline period over ALL
    % trials
    % Masker and target will be baselined TOGETHER
    % Lead and Lag will be baselined separately, for those analyses

    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset,2));
    [~,baseline_start_index] = min(abs(single_onset_time - erp_window_start_time));
    [~,baseline_end_index] = min(abs(single_onset_time - 0));

    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press,2));
    [~,baseline_start_index_buttonpress] = min(abs(single_onset_time_buttonpress - erp_window_start_time));
    [~,baseline_end_index_buttonpress] = min(abs(single_onset_time_buttonpress - 0));

    % all data
    data_by_button_press_baselined = nan(size(data_by_button_press));
    data_by_pair_onset_baselined = nan(size(data_by_pair_onset));

    for ichannel = 1:32
        % all data
        data_by_button_press_baselined(ichannel,:,:) = data_by_button_press(ichannel,:,:) - mean(data_by_button_press(ichannel,baseline_start_index_buttonpress:baseline_end_index_buttonpress,:),'all');
        data_by_pair_onset_baselined(ichannel,:,:) = data_by_pair_onset(ichannel,:,:) - mean(data_by_pair_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');

    end

    all_data_button_press(isubject,:,:) = squeeze(mean(data_by_button_press_baselined,3));
    all_data_pair_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined,3));

    %all_info_button_press(isubject).info = ERP_info_button_press;
    %all_info_target(isubject).info = ERP_info_target;

    % sort all data into conditions
    itd5_by_pair_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond))),3));

    itd15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond))),3));

    ild5_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond))),3));

    ild15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond))),3));

    %     colorbar


    %% SAVE INFO FOR THIS SUBBY
    save(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_yes_button_press.mat'),'data_by_pair_onset_baselined','data_by_button_press_baselined','ERP_info_button_press','ERP_info','-v7.3')

end


