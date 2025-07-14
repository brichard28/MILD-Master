% Primary Authors: Victoria Figarola, Benjamin Richardson 7/21/23
% Secondary Authors: Emaya Anand, Maanasa Guru Adimurthy
% EPOCHING

% ,...
curr_subject_ID = char('mild_master_29'); % char();
mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
dir_mildmaster = append(mild_master_root,'RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx');

button_fig = figure;
for isubject = 1:size(curr_subject_ID,1)

    subID = strtrim(curr_subject_ID(isubject,:));
    % Set directories
    whos_using = 'Bon';

    addpath('/home/ben/Documents/MATLAB/eeglab2023.1');
    pre_pro_epoched_data_folder = append(mild_master_root,'/prepro_epoched_data/');

    % Load dataset
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', [subID, '_ICAdone.set'], 'filepath', pre_pro_epoched_data_folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
    EEG = eeg_checkset( EEG );


    % remove extraneous triggers

    EEG.event(~ismember(string({EEG.event(:).type}), {'31231','63999'})) = [];
    EEG.urevent(~ismember([EEG.urevent(:).type],[31231,63999])) = [];



    %check trigger latency distances, remove double triggers
    distance_threshold =256;
    all_latencies = [EEG.urevent(:).latency];
    all_types = [EEG.urevent(:).type];
    all_distances = diff(all_latencies);
    num_dist_below_threshold = sum(all_distances < distance_threshold);
    disp('Below is the number of instances where triggers are too close together');
    disp(num_dist_below_threshold);
    figure; xline(all_latencies);

    
    if ismember(subID,{'mild_master_25'})
        urevents_to_remove = [2,3];
        i = 4;
    else
        urevents_to_remove = [];
        i = 1;
    end
    while i < numel(all_latencies)
        % find any latencies which are within distance_threshold of this
        % latency
        z = all_latencies - all_latencies(i);
        z(z <= 0) = inf;
        zz = find(z < distance_threshold);
        if ~isempty(zz)
            urevents_to_remove = [urevents_to_remove, zz];
        end
        i = i + 1;

    end
    urevents_to_remove = unique(urevents_to_remove);
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


    %% epoch trigger for button press, save data
    if ~ismember(subID,{'mild_master_2'})

        
        % take the first trigger
        EEG.event(1).type = "7070";
        index_button_press_start = EEG.event(1).latency;

        EEG_button_press = pop_epoch(EEG, {"7070"}, [0 600], 'newname', [subID, 'button_press'], 'epochinfo', 'yes');

        button_press_data_raw = EEG_button_press.data;





        % Isolate average button press during buttonpress
        PressTimesTable = readtable(append(mild_master_root,'/RESULTS DATA/MILD-MASTER Behavior Files/mild-master__s_',strtrim(string(subID)),'__button_press_times.csv'),'ReadVariableNames',false,'FileType','text','Range','A2:ET2');

        fs = 256;
        eeg_time = 0:1/fs:((length(button_press_data_raw) - 1)/fs);
        button_window_start_time = -1000; % 100 ms before onset of word
        button_window_end_time = 1000; % 750 ms after onset of word

        button_press_data_this_subject = [];
        for ipress = 1:width(PressTimesTable)
            this_time=strtrim(string(PressTimesTable{:,ipress}));
            this_time=erase(this_time,"'");
            this_time=erase(this_time,"[");
            this_time=erase(this_time,"]");
            this_time=str2double(this_time);

            this_time = this_time + 702;

            this_search_time = this_time/1000;
            [~,this_eeg_index_start] = min(abs(eeg_time - (this_search_time + (button_window_start_time/1000))));
            [~,this_eeg_index_end] = min(abs(eeg_time - (this_search_time + (button_window_end_time/1000))));

            if this_eeg_index_end - this_eeg_index_start == 1537
                this_eeg_index_end = this_eeg_index_end + 1;
            end

            button_press_data_this_subject = cat(3,button_press_data_this_subject,button_press_data_raw(:,this_eeg_index_start:this_eeg_index_end));
        end

        single_onset_time = linspace(button_window_start_time,button_window_end_time,size(button_press_data_this_subject,2));
        button_press_crop_start = -300;
        button_press_crop_end = 300;
        [~,index_window_start] = min(abs(single_onset_time - button_press_crop_start));
        [~,index_window_end] = min(abs(single_onset_time - button_press_crop_end));
        mean_button_press= squeeze(mean(button_press_data_this_subject(:,index_window_start:index_window_end,:),3));

        % Generate button press kernel (raised cosine, hanning window)
        button_press_window = hann(length(mean_button_press));

        button_press_kernel = [];
        for ichannel = 1:32
            button_press_kernel(ichannel,:) = button_press_window.*squeeze(mean_button_press(ichannel,:))';
        end
        figure(button_fig)
        if size(curr_subject_ID,1) > 5
            subplot(round(size(curr_subject_ID,1)/5),8,isubject)
        else
            subplot(round(size(curr_subject_ID,1)),4,isubject)
        end
        plot(linspace(button_press_crop_start,button_press_crop_end,size(button_press_kernel,2)),squeeze(mean(button_press_kernel,1)),'k')

        save(append(subID,"_button_press_data.mat"),'button_press_data_raw','button_press_window','mean_button_press','index_button_press_start','PressTimesTable');

        % remove triggers
        EEG.event(1) = [];
        EEG.urevent(1) = [];
    end

    if subID == "mild_master_29"
        EEG.event(1) = [];
        EEG.urevent(1) = [];
    end

    %% Epoch without button press subtraction
    if subID ~= "mild_master_40"
        EEG_no_button_press = pop_epoch( EEG, {"31231" , "63999"}, [-1  16], 'newname', [subID, 'all epochs'], 'epochinfo', 'yes');
    else
        EEG_no_button_press = pop_epoch( EEG, {31231 , 63999}, [-1  16], 'newname', [subID, 'all epochs'], 'epochinfo', 'yes');
    end

    EEG_no_button_press = eeg_checkset(EEG_no_button_press );
    [ALLEEG EEG_no_button_press CURRENTSET] = pop_newset(ALLEEG, EEG_no_button_press, 2, 'gui', 'off');
    EEG_no_button_press = eeg_checkset(EEG_no_button_press );
    save([pre_pro_epoched_data_folder ,subID, 'all_epoch_no_button_press.mat'], "EEG_no_button_press")

    eeglab redraw;

    %% Perform button Press subtraction

    trigger_times = [EEG.event(:).latency];
    if length(trigger_times) ~= 120
        disp('Not 120 triggeers, checking something...')
        %pause
    end

    this_subject_EEG_data = EEG.data;
    this_subject_button_press_data = zeros(size(this_subject_EEG_data));

    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == strtrim(string(curr_subject_ID(isubject,:)))); % find the rows in the spreadsheet which belong to this subject

    if length(rows_this_subject) == 119 && length(trigger_times) ==120
       trigger_times(end) = [];
       disp('119 recorded click lines')
    end
    this_subject_table = BehaviorTable(rows_this_subject,:);
    WordTimesTable = readtable(append(mild_master_root,"RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(string(curr_subject_ID(isubject,:))) + "__Word_Times.csv"));

    for itrial = 1:length(trigger_times)
        start_index_in_eeg_this_trial = EEG.event(itrial).latency;
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial

        %% ISOLATE BUTTON PRESSES
        % Find this trial button presses
        this_trial_click_times = table2array(this_subject_table(itrial,9:end));
        this_trial_click_times(isnan(this_trial_click_times)) = [];

        this_trial_click_times = this_trial_click_times/1000;
        % remove double clicks
        click_distances = diff(this_trial_click_times);
        click_distances_to_remove = find(click_distances < 0.2);
        this_trial_click_times(click_distances_to_remove + 1) = [];

        this_trial_click_times = this_trial_click_times + 0.702;

        for iclick = 1:length(this_trial_click_times) % for each target word onset...
            this_click_time = this_trial_click_times(iclick);
            resampled_search_time = round(this_click_time * fs);

            kernel_start_index = start_index_in_eeg_this_trial + resampled_search_time + round((button_press_crop_start/1000)*fs);
            kernel_end_index = start_index_in_eeg_this_trial + resampled_search_time + round((button_press_crop_end/1000)*fs);

            this_subject_button_press_data(:,kernel_start_index:kernel_end_index) = button_press_kernel;
       end


    end

    %% Epoch with button press subtraction

    EEG_yes_button_press = EEG;
    EEG_yes_button_press.data = this_subject_EEG_data - this_subject_button_press_data;
    if subID ~= "mild_master_40"
        EEG_yes_button_press = pop_epoch( EEG_yes_button_press, {"31231" , "63999"}, [-1  16], 'newname', [subID, 'all epochs'], 'epochinfo', 'yes');
    else
        EEG_yes_button_press = pop_epoch( EEG_yes_button_press, {31231 , 63999}, [-1  16], 'newname', [subID, 'all epochs'], 'epochinfo', 'yes');
    end
    EEG_yes_button_press = eeg_checkset(EEG_yes_button_press);
    [ALLEEG EEG_yes_button_press CURRENTSET] = pop_newset(ALLEEG, EEG_yes_button_press, 2, 'gui', 'off');
    EEG_yes_button_press = eeg_checkset(EEG_yes_button_press);
    save([pre_pro_epoched_data_folder ,subID, 'all_epoch_yes_button_press.mat'], "EEG_yes_button_press")

    eeglab redraw;
    %% all epochs


end