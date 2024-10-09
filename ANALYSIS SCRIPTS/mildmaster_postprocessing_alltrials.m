%% gerbilmaster_postprocessing_alltrials.m

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

curr_subject_ID =  char('fullpilot1');

% Set analysis parameters
erp_window_start_time = -100; % 100 ms before onset of word
erp_window_end_time = 750; % 750 ms after onset of word
nsubjects = size(curr_subject_ID,1);
word_length = 0.3;
num_tot_trials = 240; % look into this
frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
fs = 256;
cue_dur = 1.8;
%% For each subject.....
for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID)

    % Load word onset times data
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(subID) + "__Word_Times.csv");

    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == string(curr_subject_ID(isubject,:))); % find the rows in the spreadsheet which belong to this subject
    conditions = BehaviorTable.Condition(rows_this_subject); % conditions by trial for this subject
    condition_names = {'side=r_itd=500_az=0_mag=0_lpf=0',...
        'side=r_itd=50_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az=5_mag=0_lpf=0',...
        'side=l_itd=0_az=5_mag=1_lpf=0',...
        'side=l_itd=50_az=0_mag=0_lpf=0',...
        'side=l_itd=500_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az= 5_mag=1_lpf=0',...
        'side=l_itd=0_az=5_mag=0_lpf=0'};
    this_subject_table = BehaviorTable(rows_this_subject,:);

    % Create empty arrays for ERPs
    data_by_masker_onset = [];
    data_by_target_onset = [];
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

    ERP_info_lead_masker = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_lead_target = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_lag_masker = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});
    ERP_info_lag_target = struct('SubID',{},'Trial',{},'Word',{},'Condition',{});

    % Load EEG for this subject
    epochs_filename = join([prepro_folder,strtrim(curr_subject_ID(isubject,:)),'all_epoch.mat'],'');
    this_EEG = load(epochs_filename);
    eeg_struct_name = fieldnames(this_EEG);
    this_EEG = getfield(this_EEG,string(eeg_struct_name(1)));
    these_epochs = this_EEG.data; % 32 channels x Time x 240 trials


    % Define time vector for extracting target ERPs
    eeg_time = this_EEG.times; % in milliseconds
    resampled_audio_time = -1:1/fs:16;
    resampled_audio_time = resampled_audio_time.*1000;

    % Define time vector for extracting masker ERPs
    stimulus_length = 12; % seconds
    word_length = 0.3; % seconds

    noise_thresh = 100; % 80;

    for itrial = 1:size(this_EEG.data,3)% for each trial (should be 144)
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        if mod(itrial,10) == 0
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

        %% ISOLATE BUTTON PRESSES
        % Find this trial button presses
        this_trial_click_times = table2array(this_subject_table(itrial,9:end));
        this_trial_click_times(isnan(this_trial_click_times)) = [];
        for iclick = 1:length(this_trial_click_times) % for each target word onset...
            resampled_search_time = this_trial_click_times(iclick);
            button_press_delay = -500; % ms
            [~,start_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_start_time + button_press_delay))); %
            [~,end_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_end_time)));%


            if end_time - start_time == 345
                end_time = end_time + 1;
            end

            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(detrend(these_epochs(:,start_time:end_time,itrial))) > noise_thresh,'all')
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

        %% ISOLATE TARGET WORD ONSETS
        % Within Target Onsets

        for ionset = 1:length(this_trial_target_times) % for each target word onset...
            resampled_search_time = this_trial_target_times(ionset)*1000;
            [~,start_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_start_time))); %
            [~,end_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_end_time)));%

            if end_time - start_time == 218
                end_time = end_time -1;
            end


            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(detrend(these_epochs(:,start_time:end_time,itrial))) > noise_thresh,'all')
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
                ERP_info_lead_target.Condition = [ERP_info_lead_target.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lead_target) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lead_target(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lead_target(1).Trial = itrial;
                ERP_info_lead_target(1).Word = this_trial_target_words(ionset);
                ERP_info_lead_target(1).Condition = conditions(itrial);
            end

            % lag data
            if ~isempty(ERP_info_lag_target) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lag_target.SubID = [ERP_info_lag_target.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lag_target.Trial = [ERP_info_lag_target.Trial, itrial];
                ERP_info_lag_target.Word = [ERP_info_lag_target.Word; this_trial_target_words(ionset)];
                ERP_info_lag_target.Condition = [ERP_info_lag_target.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lag_target) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lag_target(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lag_target(1).Trial = itrial;
                ERP_info_lag_target(1).Word = this_trial_target_words(ionset);
                ERP_info_lag_target(1).Condition = conditions(itrial);
            end
        end

        %% ISOLATE MASKER WORD ONSETS
        % Background Onsets (masker onsets)
        for ionset = 1:length(this_trial_masker_times)

            resampled_search_time = this_trial_masker_times(ionset)*1000;
            [~,start_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_start_time))); %
            [~,end_time] = min(abs(resampled_audio_time - (resampled_search_time + erp_window_end_time)));%

            %resampled_search_index = (masker_time(ionset) + (1*fs) - (0.1*fs));
            %start_time = round(resampled_search_index);
            %end_time = round(start_time + floor(((erp_window_end_time - erp_window_start_time)/1000)*fs));
            if end_time - start_time == 218
                end_time = end_time -1;
            end

            % Reject epochs with amplitude above +/- 100 uV
            if any(abs(detrend(these_epochs(:,start_time:end_time,itrial))) > noise_thresh,'all')
                disp('ERP rejected')
                continue
                %add variance here
            end

            this_erp = these_epochs(:,start_time:end_time,itrial);
            data_by_masker_onset = cat(3, data_by_masker_onset,this_erp);

            % add to lead or lag ERPs, depending
            if this_trial_whether_target_lead(ionset) == 0
                data_by_lead_masker_onset = cat(3,data_by_lead_masker_onset,this_erp);
            elseif this_trial_whether_target_lead(ionset) == 1
                data_by_lag_masker_onset = cat(3,data_by_lag_masker_onset,this_erp);
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
                ERP_info_lead_masker.Condition = [ERP_info_lead_masker.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lead_masker) && this_trial_whether_target_lead(ionset) == 0
                ERP_info_lead_masker(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lead_masker(1).Trial = itrial;
                ERP_info_lead_masker(1).Word = this_trial_masker_words(ionset);
                ERP_info_lead_masker(1).Condition = conditions(itrial);
            end

            % lag data
            if ~isempty(ERP_info_lag_masker) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lag_masker.SubID = [ERP_info_lag_masker.SubID; curr_subject_ID(isubject,:)];
                ERP_info_lag_masker.Trial = [ERP_info_lag_masker.Trial, itrial];
                ERP_info_lag_masker.Word = [ERP_info_lag_masker.Word; this_trial_masker_words(ionset)];
                ERP_info_lag_masker.Condition = [ERP_info_lag_masker.Condition, conditions(itrial)];
            elseif isempty(ERP_info_lag_masker) && this_trial_whether_target_lead(ionset) == 1
                ERP_info_lag_masker(1).SubID = curr_subject_ID(isubject,:);
                ERP_info_lag_masker(1).Trial = itrial;
                ERP_info_lag_masker(1).Word = this_trial_masker_words(ionset);
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
    data_by_target_and_masker = cat(3,data_by_target_onset,data_by_masker_onset);

    % lead data
    data_by_lead_target_onset_baselined = nan(size(data_by_lead_target_onset));
    data_by_lead_masker_onset_baselined = nan(size(data_by_lead_masker_onset));
    data_by_lead_target_and_masker = cat(3,data_by_lead_target_onset,data_by_lead_masker_onset);

    % lag data
    data_by_lag_target_onset_baselined = nan(size(data_by_lag_target_onset));
    data_by_lag_masker_onset_baselined = nan(size(data_by_lag_masker_onset));
    data_by_lag_target_and_masker = cat(3,data_by_lag_target_onset,data_by_lag_masker_onset);

    for ichannel = 1:32
        % all data
        data_by_button_press_baselined(ichannel,:,:) = data_by_button_press(ichannel,:,:) - mean(data_by_button_press(ichannel,baseline_start_index_buttonpress:baseline_end_index_buttonpress,:),'all');
        data_by_target_onset_baselined(ichannel,:,:) = data_by_target_onset(ichannel,:,:) - mean(data_by_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_masker_onset_baselined(ichannel,:,:) = data_by_masker_onset(ichannel,:,:) - mean(data_by_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');

        % lead data
        data_by_lead_target_onset_baselined(ichannel,:,:) = data_by_lead_target_onset(ichannel,:,:) - mean(data_by_lead_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_lead_masker_onset_baselined(ichannel,:,:) = data_by_lead_masker_onset(ichannel,:,:) - mean(data_by_lead_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');

        % lag data
        data_by_lag_target_onset_baselined(ichannel,:,:) = data_by_lag_target_onset(ichannel,:,:) - mean(data_by_lag_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');
        data_by_lag_masker_onset_baselined(ichannel,:,:) = data_by_lag_masker_onset(ichannel,:,:) - mean(data_by_lag_target_and_masker(ichannel,baseline_start_index:baseline_end_index,:),'all');

    end

    all_data_button_press(isubject,:,:) = squeeze(mean(data_by_button_press_baselined,3));
    all_data_target(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined,3));
    all_data_masker(isubject,:,:) = squeeze(mean(data_by_masker_onset_baselined,3));

    all_data_lead_target(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined,3));
    all_data_lead_masker(isubject,:,:) = squeeze(mean(data_by_lead_masker_onset_baselined,3));

    all_data_lag_target(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined,3));
    all_data_lag_masker(isubject,:,:) = squeeze(mean(data_by_lag_masker_onset_baselined,3));

    %all_info_button_press(isubject).info = ERP_info_button_press;
    %all_info_target(isubject).info = ERP_info_target;
    %all_info_masker(isubject).info = ERP_info_masker;

    % sort all data into conditions
    itd50_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,[2,5]))),3));
    itd50_by_masker_onset(isubject,:,:) = squeeze(mean(data_by_masker_onset_baselined(:,:,logical(ismember(ERP_info_masker(:).Condition,[2,5]))),3));

    itd500_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,[1,6]))),3));
    itd500_by_masker_onset(isubject,:,:) = squeeze(mean(data_by_masker_onset_baselined(:,:,logical(ismember(ERP_info_masker(:).Condition,[1,6]))),3));

    ild5_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,[3,8]))),3));
    ild5_by_masker_onset(isubject,:,:) = squeeze(mean(data_by_masker_onset_baselined(:,:,logical(ismember(ERP_info_masker(:).Condition,[3,8]))),3));

    ild5mag_by_target_onset(isubject,:,:) = squeeze(mean(data_by_target_onset_baselined(:,:,logical(ismember(ERP_info_target(:).Condition,[4,7]))),3));
    ild5mag_by_masker_onset(isubject,:,:) = squeeze(mean(data_by_masker_onset_baselined(:,:,logical(ismember(ERP_info_masker(:).Condition,[4,7]))),3));

    % sort lead data into conditions
    itd50_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]))),3));
    itd50_by_lead_masker_onset(isubject,:,:) = squeeze(mean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,[2,5]))),3));

    itd500_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]))),3));
    itd500_by_lead_masker_onset(isubject,:,:) = squeeze(mean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,[1,6]))),3));

    ild5_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]))),3));
    ild5_by_lead_masker_onset(isubject,:,:) = squeeze(mean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,[3,8]))),3));

    ild5mag_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]))),3));
    ild5mag_by_lead_masker_onset(isubject,:,:) = squeeze(mean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,[4,7]))),3));

    % sort lag data into conditions
    itd50_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]))),3));
    itd50_by_lag_masker_onset(isubject,:,:) = squeeze(mean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,[2,5]))),3));

    itd500_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]))),3));
    itd500_by_lag_masker_onset(isubject,:,:) = squeeze(mean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,[1,6]))),3));

    ild5_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]))),3));
    ild5_by_lag_masker_onset(isubject,:,:) = squeeze(mean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,[3,8]))),3));

    ild5mag_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]))),3));
    ild5mag_by_lag_masker_onset(isubject,:,:) = squeeze(mean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,[4,7]))),3));



    %% Plot all data frontocentral erp for each subject
    figure;
    subplot(1,4,1)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd50_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd50_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('50 us ITD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)
    ylabel('Voltage (uV)','FontSize',18)

    subplot(1,4,2)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd500_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd500_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('500 us ITD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)

    subplot(1,4,3)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    subplot(1,4,4)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5mag_by_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5mag_by_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD + MAG','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    sgtitle(subID)

    %% Plot lead frontocentral ERP
    figure;
    subplot(1,4,1)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd50_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd50_by_lead_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('50 us ITD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)
    ylabel('Voltage (uV)','FontSize',18)

    subplot(1,4,2)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd500_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd500_by_lead_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('500 us ITD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)

    subplot(1,4,3)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5_by_lead_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    subplot(1,4,4)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5mag_by_lead_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5mag_by_lead_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD + MAG','FontSize',18)
    xline(0)
    xline(250)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    sgtitle([subID,' Lead ERPs'])


    %% Plot lag frontocentral ERP
    figure;
    subplot(1,4,1)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd50_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd50_by_lag_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('50 us ITD','FontSize',18)
    xline(0)
    xline(350)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)
    ylabel('Voltage (uV)','FontSize',18)

    subplot(1,4,2)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(itd500_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(itd500_by_lag_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('500 us ITD','FontSize',18)
    xline(0)
    xline(350)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)

    subplot(1,4,3)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5_by_lag_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD','FontSize',18)
    xline(0)
    xline(350)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    subplot(1,4,4)
    hold on
    p1 = plot(single_onset_time,squeeze(mean(ild5mag_by_lag_target_onset(isubject,frontocentral_channels,:),2)),'-r');
    p2 = plot(single_onset_time,squeeze(mean(ild5mag_by_lag_masker_onset(isubject,frontocentral_channels,:),2)),'-b');
    title('5 deg ILD + MAG','FontSize',18)
    xline(0)
    xline(350)
    legend([p1(1),p2(1)],{'Attend','Ignore'})
    ylim([-4,3])
    xlim([-100,500])
    xlabel('Time (ms)','FontSize',18)


    sgtitle([subID,' lag ERPs'])


    %% Plot whole trial EEG at Cz
    [~,baseline_start_index_wholetrial] = min(abs(eeg_time + 1000));
    [~,baseline_end_index_wholetrial] = min(abs(eeg_time + 0));
    these_epochs = these_epochs - mean(these_epochs(:,baseline_start_index_wholetrial:baseline_end_index_wholetrial,:),[2,3]);
    
    figure;
    subplot(size(curr_subject_ID,1),1,isubject)
    plot(eeg_time,mean(these_epochs(32,:,:),[1,3]))
    hold on;plot(eeg_time,mean(these_epochs(32,:,:),[1,3]))
    %legend({'Scrambled','Unscrambled'})
    if isubject == size(curr_subject_ID,1)
        xlabel('Time (ms)','FontSize',18)
    elseif isubject == round(size(curr_subject_ID,1)/2)
        ylabel('Voltage (mV)','FontSize',18)
    end
    xline(0,'LineWidth',2)
    xline(1800,'LineWidth',2)
    xline([1.8 2.05 2.4 2.65 3 3.25 3.6 3.85 4.2 4.45 4.8 5.05 5.4 5.65 6 6.25]*1000,'LineWidth',1)
    title(subID)

    %% Plot spectrogram of the epoch at Frontocentral and Parietooccipital channels
    % Frontocentral channels
    figure;
    M = fs/10;
    window=hamming(M,'periodic');
    noverlap=floor(0.75*M);
    Ndft=128;
    this_data = mean(these_epochs(frontocentral_channels,:,:),[1,3]);
    spectrogram(this_data,window,noverlap,Ndft,fs,'yaxis')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([-40,5])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    %yline(8,'LineWidth',2)
    %yline(12,'LineWidth',2)
    title('Frontocentral Channels','FontSize',18)

    % Parietooccipital channels ATTTEND LEFT
    parietooccipital_channels = 11:20;
    figure;
    M = fs/10;
    window=hamming(M,'periodic');
    noverlap=floor(0.75*M);
    Ndft=128;
    this_data = mean(these_epochs(parietooccipital_channels,:,:),[1,3]);
    spectrogram(this_data,window,noverlap,Ndft,fs,'yaxis')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([-40,5])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Parietooccipital Channels','FontSize',18)

     % what if I subtract the mean over the WHOLE trial (get rid of the
     % ERP?)
     parietooccipital_channels = 11:20;
    figure;
    M = fs/10;
    window=hamming(M,'periodic');
    noverlap=floor(0.75*M);
    Ndft=128;
    this_data = nan(size(these_epochs));
    for ichannel = 1:32
        this_data(ichannel,:,:) = these_epochs(ichannel,:,:) - mean(these_epochs(ichannel,:,:),3);
    end
    this_data = mean(this_data(parietooccipital_channels,:,:),[1,3]);
    spectrogram(this_data*10e7,window,noverlap,Ndft,fs,'yaxis')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([0,10])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Parietooccipital Channels WITH MEAN SUBTRACTED','FontSize',18)


    %% Topoplots in each condition
    figure;
    hold on
    cmin = -3;
    cmax = 3;
    fs = 256;

    topoplot_indices = round(0:0.05*fs:(((600)/1000)*fs));
    topoplot_indices(1) = 1;
    topoplot_times = -100:50:500;

    iplot = 1;

    % ITD50 Attend
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ITD50\newlineAttend','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(itd50_by_target_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ITD50 Ignore
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ITD50\newlineIgnore','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(itd50_by_masker_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ITD500 Attend
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ITD500\newlineAttend','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(itd500_by_target_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ITD500 Ignore
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ITD500\newlineIgnore','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(itd500_by_masker_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ILD5 Attend
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ILD5\newlineAttend','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(ild5_by_target_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ILD5 Ignore
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ILD5\newlineIgnore','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(ild5_by_masker_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ILD5Mag Attend
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ILD5 Mag\newlineAttend','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(ild5mag_by_target_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % ILD5 Ignore
    subplot(8,length(topoplot_indices)+ 1,iplot);
    text(-1,0.5,'ILD5 Mag\newlineIgnore','Interpreter','tex','FontSize',18);
    axis off
    iplot = iplot+1;

    itime = 1;
    for itopo = topoplot_indices
        subplot(8,length(topoplot_indices)+ 1,iplot);
        this_data = mean(ild5mag_by_masker_onset(:,:,itopo), [1,3]);
        topoplot(this_data,this_EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end

    % button press time trace at Cz, Fz, and Pz
        figure;
        subplot(3,1,1)
        this_data = mean(all_data_button_press(:,32,:),2);
        plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
        title('Button Press Cz','FontSize',18)
        ylim([-3,4])
        subplot(3,1,2)
        this_data = mean(all_data_button_press(:,31,:),2);
        plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
        title('Button Press Fz','FontSize',18)
        ylim([-3,4])
        subplot(3,1,3)
        this_data = mean(all_data_button_press(:,13,:),2);
        plot(single_onset_time_buttonpress,squeeze(mean(this_data,1)),'-k')
        title('Button Press Pz','FontSize',18)
        ylim([-3,4])



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
    save(append('Results_Subject_',string(curr_subject_ID(isubject,:)),'.mat'),'data_by_masker_onset_baselined','data_by_target_onset_baselined','data_by_button_press_baselined','ERP_info_button_press','ERP_info_masker','ERP_info_target','-v7.3')

end

