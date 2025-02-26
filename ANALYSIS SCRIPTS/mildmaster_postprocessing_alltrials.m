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

curr_subject_ID = char('mild_master_1',...
    'mild_master_3',...
    'mild_master_4',...
    'mild_master_6',...
    'mild_master_8',...
    'mild_master_9',...
    'mild_master_10',...
    'mild_master_11',...
    'mild_master_12',...
    'mild_master_14',...
    'mild_master_15',...
    'mild_master_16',...
    'mild_master_17',...
    'mild_master_18',...
    'mild_master_19',...
    'mild_master_22',...
    'mild_master_23'); % char();
% Set analysis parameters
erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 950; % 750 ms after onset of word
nsubjects = size(curr_subject_ID,1);
word_length = 0.3;
frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
fs = 256;

%% For each subject.....
for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID)

    cue_dur = 2.0;

    % Load word onset times data
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(subID) + "__Word_Times.csv");


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
    data_by_target_onset = [];
    data_by_masker_onset = [];
    data_by_button_press = [];

    data_by_lead_target_onset = [];
    data_by_lag_target_onset = [];
    data_by_lead_masker_onset = [];
    data_by_lag_masker_onset = [];

    % Create empty arrays for info for each ERP
    % Will contain subID, trial, and word (if target)
    ERP_info_masker = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_target = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_button_press = struct('SubID',{},'Trial',{});

    ERP_info_lead_target = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_lag_target = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});

    ERP_info_lead_masker = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_lag_masker = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});

    % Load EEG for this subject
    epochs_filename = join([prepro_folder,strtrim(curr_subject_ID(isubject,:)),'all_epoch.mat'],'');
    this_EEG = load(epochs_filename);
    eeg_struct_name = fieldnames(this_EEG);
    this_EEG = getfield(this_EEG,string(eeg_struct_name(1)));
    these_epochs = this_EEG.data; % 32 channels x Time x 240 trials
    num_tot_trials = size(these_epochs,3); % look into this


    num_trials = size(this_EEG.data,3);
    if subID == "mild_master_2"
        num_trials = 78;
    elseif ismember(subID,["mild_master_9 ","mild_master_10","mild_master_16","mild_master_18"])
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
        if ismember(subID,{'fullpilot1','fullpilot2','fullpilot3'})
            audio_onset_delay = 30;
        else
            audio_onset_delay = 0;
        end

        for ionset = 1:length(this_trial_target_times) % for each word pair onset, excluding the first
            if this_trial_whether_target_lead(ionset) == 1
                resampled_search_time = this_trial_target_times(ionset)*1000 + audio_onset_delay;
            elseif this_trial_whether_target_lead(ionset) == 0
                resampled_search_time = this_trial_masker_times(ionset)*1000 + audio_onset_delay;
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
            data_by_target_onset = cat(3,data_by_target_onset,this_erp);

            % add to lead or lag ERPs, depending
            if this_trial_whether_target_lead(ionset) == 1
                data_by_lead_target_onset = cat(3,data_by_lead_target_onset,this_erp);
            elseif this_trial_whether_target_lead(ionset) == 0
                data_by_lag_target_onset = cat(3,data_by_lag_target_onset,this_erp);
            end



            % Append Info
            if ~isempty(ERP_info_target)
                ERP_info_target.SubID = [ERP_info_target.SubID; curr_subject_ID(isubject,:)];
                ERP_info_target.Trial = [ERP_info_target.Trial, itrial];
                ERP_info_target.Word = [ERP_info_target.Word; this_trial_target_words(ionset)];
                ERP_info_target.Condition = [ERP_info_target.Condition, conditions(itrial)];
            else
                ERP_info_target(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_target(1).Trial = itrial;
                ERP_info_target(1).Word = this_trial_target_words(ionset);
                ERP_info_target(1).Condition = conditions(itrial);
            end

            % lead data

            if ~isempty(ERP_info_lead_target) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lead_target.SubID = [ERP_info_lead_target.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lead_target.Trial = [ERP_info_lead_target.Trial, itrial];
                ERP_info_lead_target.Word = [ERP_info_lead_target.Word; this_trial_target_words(ionset)];
                ERP_info_lead_target.OtherWord = [ERP_info_lead_target.OtherWord; this_trial_masker_words(ionset)];
                ERP_info_lead_target.Condition = [ERP_info_lead_target.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lead_target) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lead_target(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lead_target(1).Trial = itrial;
                ERP_info_lead_target(1).Word = this_trial_target_words(ionset);
                ERP_info_lead_target.OtherWord = this_trial_masker_words(ionset);
                ERP_info_lead_target(1).Condition = conditions(itrial);
            end

            % lag data
            if ~isempty(ERP_info_lag_target) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lag_target.SubID = [ERP_info_lag_target.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lag_target.Trial = [ERP_info_lag_target.Trial, itrial];
                ERP_info_lag_target.Word = [ERP_info_lag_target.Word; this_trial_target_words(ionset)];
                ERP_info_lag_target.OtherWord = [ERP_info_lag_target.OtherWord; this_trial_masker_words(ionset)];
                ERP_info_lag_target.Condition = [ERP_info_lag_target.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lag_target) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lag_target(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lag_target(1).Trial = itrial;
                ERP_info_lag_target(1).Word = this_trial_target_words(ionset);
                ERP_info_lag_target(1).OtherWord = this_trial_masker_words(ionset);

                ERP_info_lag_target(1).Condition = conditions(itrial);
            end
        end




        %% Same thing but for Masker
        for ionset = 1:length(this_trial_masker_times) % for each word pair onset, excluding the first
            if this_trial_whether_target_lead(ionset) == 1
                resampled_search_time = this_trial_masker_times(ionset)*1000 + audio_onset_delay;
            elseif this_trial_whether_target_lead(ionset) == 0
                resampled_search_time = this_trial_masker_times(ionset)*1000 + audio_onset_delay;
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
            data_by_masker_onset = cat(3,data_by_masker_onset,this_erp);

            % add to lead or lag ERPs, depending
            if this_trial_whether_target_lead(ionset) == 1
                data_by_lag_masker_onset = cat(3,data_by_lag_masker_onset,this_erp);
            elseif this_trial_whether_target_lead(ionset) == 0
                data_by_lead_masker_onset = cat(3,data_by_lead_masker_onset,this_erp);
            end



            % Append Info
            if ~isempty(ERP_info_masker)
                ERP_info_masker.SubID = [ERP_info_masker.SubID; curr_subject_ID(isubject,:)];
                ERP_info_masker.Trial = [ERP_info_masker.Trial, itrial];
                ERP_info_masker.Word = [ERP_info_masker.Word; this_trial_masker_words(ionset)];
                ERP_info_masker.Condition = [ERP_info_masker.Condition, conditions(itrial)];
            else
                ERP_info_masker(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_masker(1).Trial = itrial;
                ERP_info_masker(1).Word = this_trial_masker_words(ionset);
                ERP_info_masker(1).Condition = conditions(itrial);
            end

            % lead data

            if ~isempty(ERP_info_lead_masker) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lead_masker.SubID = [ERP_info_lead_masker.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lead_masker.Trial = [ERP_info_lead_masker.Trial, itrial];
                ERP_info_lead_masker.Word = [ERP_info_lead_masker.Word; this_trial_masker_words(ionset)];
                ERP_info_lead_masker.OtherWord = [ERP_info_lead_masker.OtherWord; (ionset)];
                ERP_info_lead_masker.Condition = [ERP_info_lead_masker.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lead_masker) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lead_masker(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lead_masker(1).Trial = itrial;
                ERP_info_lead_masker(1).Word = this_trial_masker_words(ionset);
                ERP_info_lead_masker.OtherWord = this_trial_masker_words(ionset);
                ERP_info_lead_masker(1).Condition = conditions(itrial);
            end

            % lag data
            if ~isempty(ERP_info_lag_masker) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lag_masker.SubID = [ERP_info_lag_masker.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lag_masker.Trial = [ERP_info_lag_masker.Trial, itrial];
                ERP_info_lag_masker.Word = [ERP_info_lag_masker.Word; this_trial_masker_words(ionset)];
                ERP_info_lag_masker.OtherWord = [ERP_info_lag_masker.OtherWord; this_trial_target_words(ionset)];
                ERP_info_lag_masker.Condition = [ERP_info_lag_masker.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lag_masker) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lag_masker(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lag_masker(1).Trial = itrial;
                ERP_info_lag_masker(1).Word = this_trial_masker_words(ionset);
                ERP_info_lag_masker(1).OtherWord = this_trial_target_words(ionset);

                ERP_info_lag_masker(1).Condition = conditions(itrial);
            end
        end

    end

    %% Concatenate and baseline within each channel for this subject
    % Baseline to the mean voltage during the baseline period over ALL
    % trials
    % Masker and target will be baselined TOGETHER
    % Lead and Lag will be baselined separately, for those analyses

    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset,2));
    [~,baseline_start_index] = min(abs(single_onset_time - erp_window_start_time));
    [~,baseline_end_index] = min(abs(single_onset_time - 0));

    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press,2));
    [~,baseline_start_index_buttonpress] = min(abs(single_onset_time_buttonpress - erp_window_start_time));
    [~,baseline_end_index_buttonpress] = min(abs(single_onset_time_buttonpress - 0));

    % all data
    data_by_button_press_baselined = nan(size(data_by_button_press));
    data_by_target_onset_baselined = nan(size(data_by_target_onset));
    data_by_masker_onset_baselined = nan(size(data_by_masker_onset));

    % lead data
    data_by_lead_target_onset_baselined = nan(size(data_by_lead_target_onset));
    data_by_lead_masker_onset_baselined = nan(size(data_by_lead_masker_onset));

    % lag data
    data_by_lag_target_onset_baselined = nan(size(data_by_lag_target_onset));
    data_by_lag_masker_onset_baselined = nan(size(data_by_lag_masker_onset));
    % BOFA
    data_by_both_lead = cat(3,data_by_lead_target_onset,data_by_lead_masker_onset);
    data_by_both_lag = cat(3,data_by_lag_target_onset,data_by_lag_masker_onset);
    for ichannel = 1:32
        % all data
        data_by_button_press_baselined(ichannel,:,:) = data_by_button_press(ichannel,:,:) - mean(data_by_button_press(ichannel,baseline_start_index_buttonpress:baseline_end_index_buttonpress,:),'all');
        data_by_target_onset_baselined(ichannel,:,:) = data_by_target_onset(ichannel,:,:) - mean(data_by_target_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_masker_onset_baselined(ichannel,:,:) = data_by_masker_onset(ichannel,:,:) - mean(data_by_masker_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');

        % lead data
        data_by_lead_target_onset_baselined(ichannel,:,:) = data_by_lead_target_onset(ichannel,:,:) - mean(data_by_lead_target_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_lead_masker_onset_baselined(ichannel,:,:) = data_by_lead_masker_onset(ichannel,:,:) - mean(data_by_lead_masker_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');

        % lag data
        data_by_lag_target_onset_baselined(ichannel,:,:) = data_by_lag_target_onset(ichannel,:,:) - mean(data_by_lag_target_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_lag_masker_onset_baselined(ichannel,:,:) = data_by_lag_masker_onset(ichannel,:,:) - mean(data_by_lag_masker_onset(ichannel,baseline_start_index:baseline_end_index,:),'all');

    end

    all_data_button_press(isubject,:,:) = squeeze(mean(data_by_button_press_baselined,3));
    all_data_target(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined,3));

    all_data_lead_target(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined,3));

    all_data_lag_target(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined,3));

    %all_info_button_press(isubject).info = ERP_info_button_press;
    %all_info_target(isubject).info = ERP_info_target;

    % sort all data into conditions
    itd5_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,small_itd_cond))),3));

    itd15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,large_itd_cond))),3));

    ild5_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,small_ild_cond))),3));

    ild15_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,large_ild_cond))),3));

    %% ALL WORDS
    % sort lead data into conditions
    itd5_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond))),3));
    itd15_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond))),3));
    ild5_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond))),3));
    ild15_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond))),3));

    % sort lag data into conditions
    itd5_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond))),3));
    itd15_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond))),3));
    ild5_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond))),3));
    ild15_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond))),3));

    %% JUST BASH

    itd5_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond).*ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    itd15_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond).*ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    ild5_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond).*ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    ild15_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond).*ismember(ERP_info_lead_target(:).Word,"bash")')),3));

    % sort lag data into conditions
    itd5_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond).*ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    itd15_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond).*ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    ild5_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond).*ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    ild15_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond).*ismember(ERP_info_lag_target(:).Word,"bash")')),3));


    %% JUST DASH, GASH
    itd5_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond).*~ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    itd15_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond).*~ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    ild5_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond).*~ismember(ERP_info_lead_target(:).Word,"bash")')),3));
    ild15_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond).*~ismember(ERP_info_lead_target(:).Word,"bash")')),3));

    % sort lag data into conditions
    itd5_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond).*~ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    itd15_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond).*~ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    ild5_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond).*~ismember(ERP_info_lag_target(:).Word,"bash")')),3));
    ild15_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond).*~ismember(ERP_info_lag_target(:).Word,"bash")')),3));

    %% Plot all data frontocentral erp for each subject
    %     figure;
    %     subplot(1,4,1)
    %     hold on
    %     p1 = plot(single_onset_time,squeeze(mean(itd5_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    %     p2 = plot(single_onset_time,squeeze(mean(itd5_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    %     title('5 deg ITDs','FontSize',18)
    %     xline(0)
    %     xline(250)
    %     legend([p1(1),p2(1)],{'Attend','Ignore'})
    %     ylim([-4,3])
    %     xlim([-100,750])
    %     xlabel('Time (ms)','FontSize',18)
    %     ylabel('Voltage (uV)','FontSize',18)
    %
    %     subplot(1,4,2)
    %     hold on
    %     p1 = plot(single_onset_time,squeeze(mean(itd15_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    %     p2 = plot(single_onset_time,squeeze(mean(itd15_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    %     title('15 deg ITDs','FontSize',18)
    %     xline(0)
    %     xline(250)
    %     legend([p1(1),p2(1)],{'Attend','Ignore'})
    %     ylim([-4,3])
    %     xlim([-100,750])
    %     xlabel('Time (ms)','FontSize',18)
    %
    %     subplot(1,4,3)
    %     hold on
    %     p1 = plot(single_onset_time,squeeze(mean(ild5_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    %     p2 = plot(single_onset_time,squeeze(mean(ild5_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    %     title('5 deg ILDs','FontSize',18)
    %     xline(0)
    %     xline(250)
    %     legend([p1(1),p2(1)],{'Attend','Ignore'})
    %     ylim([-4,3])
    %     xlim([-100,750])
    %     xlabel('Time (ms)','FontSize',18)
    %
    %
    %     subplot(1,4,4)
    %     hold on
    %     p1 = plot(single_onset_time,squeeze(mean(ild15_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    %     p2 = plot(single_onset_time,squeeze(mean(ild15_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    %     title('15 deg ILDs','FontSize',18)
    %     xline(0)
    %     xline(250)
    %     legend([p1(1),p2(1)],{'Attend','Ignore'})
    %     ylim([-4,3])
    %     xlim([-100,750])
    %     xlabel('Time (ms)','FontSize',18)
    %
    %
    %     sgtitle(subID)

    %% Plot frontocentral ERPs
%     y_min = -4;
%     y_max = 5;
%     figure;
%     subplot(1,4,1)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd5_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd5_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
%     ylabel('Voltage (uV)','FontSize',18)
% 
%     subplot(1,4,2)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd15_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd15_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
%     subplot(1,4,3)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild5_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild5_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     subplot(1,4,4)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild15_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild15_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     sgtitle([subID,' Frontocentral ERPs'])
% 
%     %% Plot Bash separately
%     y_min = -4;
%     y_max = 5;
%     figure;
%     subplot(1,4,1)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd5_by_lead_target_onset_bash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd5_by_lag_target_onset_bash(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
%     ylabel('Voltage (uV)','FontSize',18)
% 
%     subplot(1,4,2)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd15_by_lead_target_onset_bash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd15_by_lag_target_onset_bash(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
%     subplot(1,4,3)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild5_by_lead_target_onset_bash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild5_by_lag_target_onset_bash(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     subplot(1,4,4)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild15_by_lead_target_onset_bash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild15_by_lag_target_onset_bash(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     sgtitle([subID,' Frontocentral ERPs BASH ONLY'])
% 
% 
%     %% Plot dash, gash
% 
%     y_min = -4;
%     y_max = 5;
%     figure;
%     subplot(1,4,1)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd5_by_lead_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd5_by_lag_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
%     ylabel('Voltage (uV)','FontSize',18)
% 
%     subplot(1,4,2)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(itd15_by_lead_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(itd15_by_lag_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ITDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
%     subplot(1,4,3)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild5_by_lead_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild5_by_lag_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-b');
%     title('5 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     subplot(1,4,4)
%     hold on
%     p1 = plot(single_onset_time,squeeze(mean(ild15_by_lead_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-r');
%     p2 = plot(single_onset_time,squeeze(mean(ild15_by_lag_target_onset_dashgash(isubject,frontocentral_channels,:),2)),'-b');
%     title('15 deg ILDs','FontSize',18)
%     xline(0,'--r','LineWidth',2)
%     xline(250,'--b','LineWidth',2)
%     xline(600)
%     legend([p1(1),p2(1)],{'Attend Lead','Attend Lag'})
%     ylim([y_min,y_max])
%     xlim([-100,750])
%     xlabel('Time (ms)','FontSize',18)
% 
% 
%     sgtitle([subID,' Frontocentral ERPs {DASH,GASH} ONLY'])
% 
% 

    %% Plot whole trial EEG at frontocentral channels in each condition
%     [~,baseline_start_index_wholetrial] = min(abs(eeg_time + 1000));
%     [~,baseline_end_index_wholetrial] = min(abs(eeg_time + 0));
%     these_epochs = these_epochs - mean(these_epochs(:,baseline_start_index_wholetrial:baseline_end_index_wholetrial,:),[2,3]);
% 
%     figure;
%     %subplot(size(curr_subject_ID,1),1,isubject)
%     hold on;
%     p1 = plot(eeg_time,mean(these_epochs(frontocentral_channels,:,ismember(conditions,small_itd_cond)),[1,3])); % itd 50
%     p2 = plot(eeg_time,mean(these_epochs(frontocentral_channels,:,ismember(conditions,large_itd_cond)),[1,3])); % itd 500
%     p3 = plot(eeg_time,mean(these_epochs(frontocentral_channels,:,ismember(conditions,small_ild_cond)),[1,3])); % ild 5
%     p4 = plot(eeg_time,mean(these_epochs(frontocentral_channels,:,ismember(conditions,large_ild_cond)),[1,3])); % itd 5 + mag
% 
%     %legend({'Scrambled','Unscrambled'})
%     if isubject == size(curr_subject_ID,1)
%         xlabel('Time (ms)','FontSize',18)
%     elseif isubject == round(size(curr_subject_ID,1)/2)
%         ylabel('Voltage (mV)','FontSize',18)
%     end
%     xline(0,'LineWidth',2)
%     xline(cue_dur*1000,'LineWidth',2)
%     xline([this_trial_target_times, this_trial_masker_times]*1000,'LineWidth',1)
%     title(subID)
%     legend([p1(1),p2(1),p3(1),p4(1)],{'itd5','itd15','ILD5deg','ILD5degMag'})

    %% Calculate CWT spectrograms at each trial and channel
%     parietooccipital_channels = 11:20;
%     left_parietooccipital_channels = [11,12,14,15];
%     right_parietooccipital_channels = [17,18,19,20];
%     frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
% 
%     M = 128;%fs*0.4;
%     window=hamming(M,'periodic');
%     noverlap=floor(0.75*M);
% 
%     Ndft=128;
%     epoch_start_time = 0;
%     epoch_end_time = 16;
%     time_window = linspace(epoch_start_time, epoch_end_time, size(these_epochs, 2));
% 
%     all_spectrograms_this_subject = nan(Ndft/2 + 1, 123, 32, size(these_epochs,3));
%     all_cwt_this_subject= cell(32,num_trials);
%     freq_range= [1 50];
    %     for ichannel=1:32
    %         for itrial=1:num_tot_trials
    %             spect_sub= these_epochs(ichannel,:,itrial);
    %             %[this_spectrogram,frequencies,time_window]=spectrogram(spect_sub,window,noverlap,Ndft,fs,'yaxis');
    %             %all_spectrograms_this_subject(:,:,ichannel, itrial) = abs(this_spectrogram);
    %             [cwt_coeffs, frequencies] = cwt(spect_sub, fs, 'FrequencyLimits', freq_range);
    %             all_cwt_this_subject{ichannel,itrial}=abs(cwt_coeffs);
    %         end
    %    end
    %% Plot  CWT spectrogram of the epoch at Frontocentral and Parietooccipital channels


    % Frontocentral
    %     fc_cwt_to_plot = [];
    %     for ichannel=1:length(frontocentral_channels)
    %         for itrial = 1:num_tot_trials
    %             fc_cwt_to_plot(:,:,ichannel,itrial) = all_cwt_this_subject{frontocentral_channels(ichannel),itrial};
    %         end
    %     end
%     cmin = 50;
%     cmax = 125;
%     plot_time_limits = [-1, 6];
    %     % Frontocentral channels
    %     figure;
    %     hold on
    %     imagesc(time_window,frequencies,squeeze(mean(fc_cwt_to_plot(:,:,:,:),[3,4])))
    %     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
    %     axis tight
    %     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
    %     %caxis([cmin,cmax])
    %     xline(0,'LineWidth',2)
    %     xline(2,'LineWidth',2)
    %     %yline(8,'LineWidth',2)
    %     %yline(12,'LineWidth',2)
    %     xlim([plot_time_limits(1),plot_time_limits(end)])
    %     title('Frontocentral Channels','FontSize',18)
    %     xlabel('Time (s)','FontSize',18)
    %     ylabel('Frequency (Hz)','FontSize',18)

    % Parietooccipital channels BROKEN UP BY ATTEND
    %     po_cwt_to_plot = [];
    %     for ichannel=1:length(parietooccipital_channels)
    %         for itrial = 1:num_tot_trials
    %             po_cwt_to_plot(:,:,ichannel,itrial) = all_cwt_this_subject{parietooccipital_channels(ichannel),itrial};
    %         end
    %     end
    %
    %
    %
    %
    %     attend_right_conditions = ismember(conditions,[1,2,3,7]);
    %     attend_left_conditions = ismember(conditions,[4,5,6,8]);
    %     figure;
    %     subplot(2,2,1) % left hemisphere, attend left
    %     imagesc(time_window,frequencies,squeeze(mean(po_cwt_to_plot(:,:,left_parietooccipital_channels - 10,ismember(conditions,attend_left_conditions)),[3,4])))
    %     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
    %     axis tight
    %     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
    %     %caxis([cmin,cmax])
    %     xline(0,'LineWidth',2)
    %     xline(2,'LineWidth',2)
    %     %yline(8,'LineWidth',2)
    %     %yline(12,'LineWidth',2)
    %     title('Left Hem. Attend Left','FontSize',18)
    %     xlabel('Time (s)','FontSize',18)
    %     ylabel('Frequency (Hz)','FontSize',18)
    %     xlim([plot_time_limits(1),plot_time_limits(end)])
    %
    %
    %     subplot(2,2,2) % right hemisphere, attend left
    %     imagesc(time_window,frequencies',squeeze(mean(po_cwt_to_plot(:,:,right_parietooccipital_channels - 10,attend_left_conditions),[3,4])))
    %     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
    %     axis tight
    %     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
    %     %caxis([cmin,cmax])
    %     xline(0,'LineWidth',2)
    %     xline(2,'LineWidth',2)
    %     %yline(8,'LineWidth',2)
    %     %yline(12,'LineWidth',2)
    %     title('Right Hem. Attend Left','FontSize',18)
    %     xlabel('Time (s)','FontSize',18)
    %     ylabel('Frequency (Hz)','FontSize',18)
    %     xlim([plot_time_limits(1),plot_time_limits(end)])
    %
    %
    %     subplot(2,2,3) % left hemisphere, attend right
    %     imagesc(time_window,frequencies',squeeze(mean(po_cwt_to_plot(:,:,left_parietooccipital_channels - 10,attend_right_conditions),[3,4])))
    %     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
    %     axis tight
    %     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
    %     %caxis([cmin,cmax])
    %     xline(0,'LineWidth',2)
    %     xline(2,'LineWidth',2)
    %     %yline(8,'LineWidth',2)
    %     %yline(12,'LineWidth',2)
    %     title('Left Hem. Attend Right','FontSize',18)
    %     xlabel('Time (s)','FontSize',18)
    %     ylabel('Frequency (Hz)','FontSize',18)
    %     xlim([plot_time_limits(1),plot_time_limits(end)])
    %
    %
    %     subplot(2,2,4) % right hemisphere, attend right
    %     imagesc(time_window,frequencies',squeeze(mean(po_cwt_to_plot(:,:,right_parietooccipital_channels - 10,attend_right_conditions),[3,4])))
    %     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
    %     axis tight
    %     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
    %     %caxis([cmin,cmax])
    %     xline(0,'LineWidth',2)
    %     xline(2,'LineWidth',2)
    %     %yline(8,'LineWidth',2)
    %     %yline(12,'LineWidth',2)
    %     title('Right Hem. Attend Right','FontSize',18)
    %     xlabel('Time (s)','FontSize',18)
    %     ylabel('Frequency (Hz)','FontSize',18)
    %     xlim([plot_time_limits(1),plot_time_limits(end)])
    %
    %
    %
    %     sgtitle('Parietooccipital Channels','FontSize',18)

    % what if I subtract the mean over the WHOLE trial (get rid of the
    % ERP?)
    %      parietooccipital_channels = 11:20;
    %     figure;
    %     M = fs/10;
    %     window=hamming(M,'periodic');
    %     noverlap=floor(0.75*M);
    %     Ndft=128;
    %     this_data = nan(size(these_epochs));
    %     for ichannel = 1:32
    %         this_data(ichannel,:,:) = these_epochs(ichannel,:,:) - mean(these_epochs(ichannel,:,:),3);
    %     end
    %     this_data = mean(this_data(parietooccipital_channels,:,:),[1,3]);
    %     spectrogram(this_data*10e7,window,noverlap,Ndft,fs,'yaxis')
    %     ylim([0,30])
    %     xlim([0,6])
    %     xticklabels(0:1:16);
    %     curr_tick_labels = xticklabels;
    %     xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    %     caxis([0,10])
    %     xline(1,'LineWidth',2)
    %     xline(2.8,'LineWidth',2)
    %     yline(8,'LineWidth',2)
    %     yline(12,'LineWidth',2)
    %     title('Parietooccipital Channels WITH MEAN SUBTRACTED','FontSize',18)


    %% Whole Trial topoplot
    %     figure;
    %     hold on
    %     cmin = -3;
    %     cmax = 3;
    %     fs = 256;
    %
    %     topoplot_indices = round(0:0.25*fs:(((13000)/1000)*fs));
    %     topoplot_indices(1) = 1;
    %     topoplot_times = -1000:250:13000;
    %
    %     iplot = 1;
    %
    %     itime = 1;
    %     for itopo = topoplot_indices
    %         subplot(4,round((length(topoplot_indices))/4),iplot);
    %         this_data = mean(these_epochs(:,itopo,:), [2,3]);
    %         topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
    %         title([num2str(topoplot_times(itime)),' ms'])
    %         iplot = iplot + 1;
    %         itime = itime + 1;
    %     end
    %
    %     sgtitle(subID)
    %
    %
    %
    %     %button press time trace at Cz, Fz, and Pz
    %     figure;
    %     subplot(3,1,1)
    %     this_data = mean(all_data_button_press(:,32,:),2);
    %     plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
    %     title('Button Press Cz','FontSize',18)
    %     ylim([-3,4])
    %     subplot(3,1,2)
    %     this_data = mean(all_data_button_press(:,31,:),2);
    %     plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
    %     title('Button Press Fz','FontSize',18)
    %     ylim([-3,4])
    %     subplot(3,1,3)
    %     this_data = mean(all_data_button_press(:,13,:),2);
    %     plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
    %     title('Button Press Pz','FontSize',18)
    %     ylim([-3,4])



    % button press topoplot
    %     figure;
    %     topoplot_indices = round(0:0.1*fs:(((erp_window_end_time - (erp_window_start_time + button_press_delay))/1000)*fs));
    %     topoplot_indices(1) = 1;
    %     topoplot_times = -500:100:800;
    %
    %     iplot = 1;
    %     itime = 1;
    %     for itopo = topoplot_indices
    %         subplot(1,length(topoplot_indices),iplot);
    %         this_data = squeeze(mean(all_data_button_press(:,:,itopo),[1,3]));
    %         topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
    %         title([num2str(topoplot_times(itime)),' ms'])
    %         iplot = iplot + 1;
    %         itime = itime + 1;
    %     end
    %     colorbar


    %% SAVE INFO FOR THIS SUBBY
    save(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'),'data_by_target_onset_baselined','data_by_lead_target_onset_baselined','data_by_lag_target_onset_baselined','data_by_masker_onset_baselined','data_by_lead_masker_onset_baselined','data_by_lag_masker_onset_baselined','data_by_button_press_baselined','ERP_info_button_press','ERP_info_lead_target','ERP_info_lag_target','ERP_info_target','ERP_info_lead_masker','ERP_info_lag_masker','ERP_info_masker','-v7.3')

end

