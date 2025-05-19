%% mildmaster_eeg_calculations.m

% Calculate P1, N1, and P300 statistics for mild master
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
    'mild_master_20',...
    'mild_master_22',...
    'mild_master_23',...
    'mild_master_24',...
    'mild_master_25',...
    'mild_master_26',...
    'mild_master_27',...
    'mild_master_28',...
    'mild_master_29',...
    'mild_master_30',...
    'mild_master_31',...
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40', ...
    'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_46','mild_master_48'); % char();
num_subjects = size(curr_subject_ID,1);

mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
erp_window_start_time = -25; % 100 ms before onset of word
erp_window_end_time = 1200; % 750 ms after onset of word
fs = 256;

all_subs_fig_p1n1 = figure();
all_subs_fig_lead_p3 = figure();
all_subs_fig_lag_p3 = figure();

all_subs_p1 = [];
all_subs_n1 = [];
all_subs_p3 = [];


all_subs_p1_target_left_onset = [];
all_subs_p1_target_right_onset = [];
all_subs_p1_masker_left_onset = [];
all_subs_p1_masker_right_onset = [];

all_subs_n1_target_left_onset = [];
all_subs_n1_target_right_onset = [];
all_subs_n1_masker_left_onset = [];
all_subs_n1_masker_right_onset = [];

all_subs_p3_target_left_onset = [];
all_subs_p3_target_right_onset = [];
all_subs_p3_masker_left_onset = [];
all_subs_p3_masker_right_onset = [];
    
for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_yes_button_press.mat'))

    %% Plot all channels, remove noisy ones in time domain

    % Plotting parameters
    button_press_delay = 0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    parietooccipital_channels = 11:20;
    cz_index = 32;
    fz_index = 31;
    pz_index = 13;


    %% Find this subjects' P1, N1, and P300 latency using the mean over ALL lead lag tokens
    % Do this separately for lead, lag P1 and N1
    % We will use Cz and Fz average to calculate P1 and N1, and Pz to
    % calculate P300

    % For P300, only basing on trials that have bash!

    % Lead P1: local maximum betweeen 50 and 150 ms
    % Lead N1: local minimum between 100 and 200
    % Lead P3: local maximum between 400 and 600 ms
    [~,lead_p1_start_index] = min(abs(single_onset_time - (50)));
    [~,lead_p1_end_index] = min(abs(single_onset_time - (150)));
    [~,lead_n1_start_index] = min(abs(single_onset_time - (100)));
    [~,lead_n1_end_index] = min(abs(single_onset_time - (200)));
    [~,lead_p3_start_index] = min(abs(single_onset_time - (450)));
    [~,lead_p3_end_index] = min(abs(single_onset_time - (550)));

    % Lag P1: local maximum betweeen 300 and 400 ms
    % Lag N1: local minimum between 350 and 450 ms
    % Lag P3: local maximum between 650 and 850 ms
    [~,lag_p1_start_index] = min(abs(single_onset_time - (300)));
    [~,lag_p1_end_index] = min(abs(single_onset_time - (400)));
    [~,lag_n1_start_index] = min(abs(single_onset_time - (350)));
    [~,lag_n1_end_index] = min(abs(single_onset_time - (450)));
    [~,lag_p3_start_index] = min(abs(single_onset_time - (700)));
    [~,lag_p3_end_index] = min(abs(single_onset_time - (800)));

    % take the average
    this_sub_cz_fz_average = squeeze(mean(data_by_pair_onset_baselined([frontocentral_channels],:,:),[1,3]));
    this_sub_lead_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lead_Word,{'bash'})),[1,3]));
    this_sub_lag_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lag_Word,{'bash'})),[1,3]));

    % lead p1
    [~,p1_locs, ~, p1_proms] = findpeaks(this_sub_cz_fz_average(lead_p1_start_index:lead_p1_end_index));
    [~,which_peak] = max(p1_proms);
    this_sub_lead_p1_index = p1_locs(which_peak) + lead_p1_start_index - 1;
    this_sub_lead_p1_time = single_onset_time(this_sub_lead_p1_index);

    % lead n1
    [~,n1_locs, ~, n1_proms] = findpeaks(-1*this_sub_cz_fz_average(lead_n1_start_index:lead_n1_end_index));
    [~,which_peak] = max(n1_proms);
    this_sub_lead_n1_index = n1_locs(which_peak) + lead_n1_start_index - 1;
    this_sub_lead_n1_time = single_onset_time(this_sub_lead_n1_index);

    % lead p3
    [~,p3_locs, ~, p3_proms] = findpeaks(this_sub_lead_bash_pz_average(lead_p3_start_index:lead_p3_end_index));
    [~,which_peak] = max(p3_proms);
    this_sub_lead_p3_index = p3_locs(which_peak) + lead_p3_start_index - 1;
    this_sub_lead_p3_time = single_onset_time(this_sub_lead_p3_index);

    % lag p1
    [~,p1_locs, ~, p1_proms] = findpeaks(this_sub_cz_fz_average(lag_p1_start_index:lag_p1_end_index));
    [~,which_peak] = max(p1_proms);
    this_sub_lag_p1_index = p1_locs(which_peak) + lag_p1_start_index - 1;
    this_sub_lag_p1_time = single_onset_time(this_sub_lag_p1_index);

    % lag n1
    [~,n1_locs, ~, n1_proms] = findpeaks(-1*this_sub_cz_fz_average(lag_n1_start_index:lag_n1_end_index));
    [~,which_peak] = max(n1_proms);
    this_sub_lag_n1_index = n1_locs(which_peak) + lag_n1_start_index - 1;
    this_sub_lag_n1_time = single_onset_time(this_sub_lag_n1_index);

    % lag p3
    [~,p3_locs, ~, p3_proms] = findpeaks(this_sub_lag_bash_pz_average(lag_p3_start_index:lag_p3_end_index));
    [~,which_peak] = max(p3_proms);
    this_sub_lag_p3_index = p3_locs(which_peak) + lag_p3_start_index - 1;
    this_sub_lag_p3_time = single_onset_time(this_sub_lag_p3_index);

    figure(all_subs_fig_p1n1)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_cz_fz_average,'k')
    scatter(this_sub_lead_p1_time,this_sub_cz_fz_average(this_sub_lead_p1_index),'or','filled');
    scatter(this_sub_lead_n1_time,this_sub_cz_fz_average(this_sub_lead_n1_index),'ob','filled');
    scatter(this_sub_lag_p1_time,this_sub_cz_fz_average(this_sub_lag_p1_index),'og','filled');
    scatter(this_sub_lag_n1_time,this_sub_cz_fz_average(this_sub_lag_n1_index),'om','filled');

    figure(all_subs_fig_lead_p3)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_lead_bash_pz_average,'k')
    scatter(this_sub_lead_p3_time,this_sub_lead_bash_pz_average(this_sub_lead_p3_index),'or','filled');
    title('Lead P3 (bash in leading position ONLY)')

    figure(all_subs_fig_lag_p3)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_lag_bash_pz_average,'k')
    scatter(this_sub_lag_p3_time,this_sub_lag_bash_pz_average(this_sub_lag_p3_index),'og','filled');
    title('lag P3 (bash in lagging position ONLY)')


    all_subs_lead_p1_times(isubject,:) = this_sub_lead_p1_time;
    all_subs_lag_p1_times(isubject,:) = this_sub_lag_p1_time;

    all_subs_lead_n1_times(isubject,:) = this_sub_lead_n1_time;
    all_subs_lag_n1_times(isubject,:) = this_sub_lag_n1_time;

    all_subs_lead_p3_times(isubject,:) = this_sub_lead_p3_time;
    all_subs_lag_p3_times(isubject,:) = this_sub_lag_p3_time;

    all_subs_lead_p1_indices(isubject,:) = this_sub_lead_p1_index;
    all_subs_lag_p1_indices(isubject,:) = this_sub_lag_p1_index;

    all_subs_lead_n1_indices(isubject,:) = this_sub_lead_n1_index;
    all_subs_lag_n1_indices(isubject,:) = this_sub_lag_n1_index;

    all_subs_lead_p3_indices(isubject,:) = this_sub_lead_p3_index;
    all_subs_lag_p3_indices(isubject,:) = this_sub_lag_p3_index;

    % Build data arrays to plot later
    % Num subjects x num conditions x num situations x num channels (32)
    % Condition order: small ITDs, large ITDs, small ILDs, large ILDs
    % ORGANIZATION:
    % (1) lead response, bash leads, target leads
    % (2) lead response, bash lags, target leads
    % (3) lead response, no bash, target leads
    % (4) lag response, bash leads, target leads
    % (5) lag response, bash lags, target leads
    % (6) lag response, no bash, target leads
    % (7) lead response, bash leads, masker leads
    % (8) lead response, bash lags, masker leads
    % (9) lead response, no bash, masker leads
    % (10) lag response, bash leads, masker leads
    % (11) lag response, bash lags, masker leads
    % (12) lag response, no bash, masker leads

    small_itd_cond = [3,7];
    large_itd_cond = [6,8];
    small_ild_cond = [1,4];
    large_ild_cond = [2,5];

    all_conditions = [small_itd_cond;large_itd_cond;small_ild_cond;large_ild_cond];

    % p1
    peak_integration_time = 0.010; % s
    for icondition = 1:4
        these_conditions = all_conditions(icondition,:);

        % should these be trapz instead?? to integreate over time?
        all_subs_p1(isubject,icondition,1,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,2,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p1(isubject,icondition,3,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,4,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,5,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p1(isubject,icondition,6,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));

        all_subs_p1(isubject,icondition,7,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,8,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p1(isubject,icondition,9,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,10,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p1(isubject,icondition,11,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p1(isubject,icondition,12,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
    
    
        % Break up data by direction
        target_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        target_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));   
        masker_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        masker_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));
        

        all_subs_p1_target_left_onset(isubject,icondition,:) = squeeze(nanmean(target_left_onset_data,[2,3]));
        all_subs_p1_target_right_onset(isubject,icondition,:) = squeeze(nanmean(target_right_onset_data,[2,3]));
        all_subs_p1_masker_left_onset(isubject,icondition,:) = squeeze(nanmean(masker_left_onset_data,[2,3]));
        all_subs_p1_masker_right_onset(isubject,icondition,:) = squeeze(nanmean(masker_right_onset_data,[2,3]));
        
    end

    % n1
    for icondition = 1:4
        these_conditions = all_conditions(icondition,:);
        all_subs_n1(isubject,icondition,1,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,2,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_n1(isubject,icondition,3,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,4,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,5,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_n1(isubject,icondition,6,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));

        all_subs_n1(isubject,icondition,7,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,8,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_n1(isubject,icondition,9,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,10,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_n1(isubject,icondition,11,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_n1(isubject,icondition,12,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));

        % Break up data by direction
        target_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        target_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));   
        masker_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        masker_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));
        

        all_subs_n1_target_left_onset(isubject,icondition,:) = squeeze(nanmean(target_left_onset_data,[2,3]));
        all_subs_n1_target_right_onset(isubject,icondition,:) = squeeze(nanmean(target_right_onset_data,[2,3]));
        all_subs_n1_masker_left_onset(isubject,icondition,:) = squeeze(nanmean(masker_left_onset_data,[2,3]));
        all_subs_n1_masker_right_onset(isubject,icondition,:) = squeeze(nanmean(masker_right_onset_data,[2,3]));
    
    
    end
    % p3
    for icondition = 1:4
        these_conditions = all_conditions(icondition,:);
        all_subs_p3(isubject,icondition,1,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,2,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p3(isubject,icondition,3,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,4,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,5,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p3(isubject,icondition,6,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));

        all_subs_p3(isubject,icondition,7,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,8,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p3(isubject,icondition,9,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p3_index  - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,10,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
        all_subs_p3(isubject,icondition,11,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),[2,3]));
        all_subs_p3(isubject,icondition,12,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),[2,3]));
    
        % Break up data by direction
        target_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        target_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));   
        masker_left_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),data_by_pair_onset_baselined(:,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'L'}))));
        masker_right_onset_data = cat(3,data_by_pair_onset_baselined(:,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),data_by_pair_onset_baselined(:,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Target_Direction,{'R'}))));
        

        all_subs_p3_target_left_onset(isubject,icondition,:) = squeeze(nanmean(target_left_onset_data,[2,3]));
        all_subs_p3_target_right_onset(isubject,icondition,:) = squeeze(nanmean(target_right_onset_data,[2,3]));
        all_subs_p3_masker_left_onset(isubject,icondition,:) = squeeze(nanmean(masker_left_onset_data,[2,3]));
        all_subs_p3_masker_right_onset(isubject,icondition,:) = squeeze(nanmean(masker_right_onset_data,[2,3]));
    end

end

save('ERP_individual_indices.mat',"all_subs_lead_p1_times","all_subs_lag_p1_times","all_subs_lag_p3_indices","all_subs_lead_p3_indices","all_subs_lag_n1_indices","all_subs_lead_n1_indices")

%% PLOT p1 and n1 topoplots bash lead vs. bash lag vs. no bash (left panel lead, right panel lag)
% eeglab
% figure;
% subplot(1,2,1) % lead responses
% topoplot()
%
%
% subplot(1,2,2) % lag responses

%% PLOT p3 topoplots  bash lead vs. bash lag vs. no bash (left panel lead, right panel lag)


%% PLOT errorbar P1-N1 difference bash lead vs. bash lag vs. no bash (left panel lead, right panel lag)

all_subs_p1_minus_n1 = all_subs_p1 - all_subs_n1;


condition_titles = {'Small ITD','Large ITD','Small ILD','Large ILD'};

for icondition = 1:4
    figure;
    this_data = squeeze(mean(all_subs_p1_minus_n1(:,icondition,:,[frontocentral_channels]),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);


    subplot(1,2,1)
    hold on
    errorbar(0.9:1:3.1,this_mean(1:3),this_sem(1:3),'k','LineStyle','none','Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b'); % Target Leads
    errorbar(1.1:1:3.1,this_mean(7:9),this_sem(7:9),'k','LineStyle','none','Marker','o','MarkerEdgeColor','g','MarkerFaceColor','g'); % Masker Leads
    title('Lead Word Response')
    xticks(1:3)
    xticklabels({'Bash Leads','Bash Lags','No Bash'})
    legend({'Target Leads (Masker Lags)','Masker Leads (Target Lags)'})
    ylabel('P1-N1 (uV)')
    ylim([0,4])

    subplot(1,2,2)
    hold on
    errorbar(0.9:1:3.1,this_mean(4:6),this_sem(4:6),'k','LineStyle','none','Marker','o','MarkerEdgeColor','g','MarkerFaceColor','g'); % Target Leads (Masker Lags)
    errorbar(1.1:1:3.1,this_mean(10:12),this_sem(10:12),'k','LineStyle','none','Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b'); % Masker Leads (Target Lags)
    title('Lag Word Response')
    xticks(1:3)
    xticklabels({'Bash Leads','Bash Lags','No Bash'})
    legend({'Target Leads (Masker Lags)','Masker Leads (Target Lags)'})
    ylabel('P1-N1 (uV)')
    ylim([0,4])

    sgtitle(append(condition_titles(icondition),' P1 N1 Differences'))
end
%% PLOT errorbar P3 bash lead vs. bash lag vs. no bash (left panel lead, right panel lag)



for icondition = 1:4
    figure;
    this_data = squeeze(mean(all_subs_p3(:,icondition,:,[parietooccipital_channels]),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);


    subplot(1,2,1)
    hold on
    errorbar(0.9:1:3.1,this_mean(1:3),this_sem(1:3),'k','LineStyle','none','Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b'); % Target Leads
    errorbar(1.1:1:3.1,this_mean(7:9),this_sem(7:9),'k','LineStyle','none','Marker','o','MarkerEdgeColor','g','MarkerFaceColor','g'); % Masker Leads
    title('Lead Word Response')
    xticks(1:3)
    xticklabels({'Bash Leads','Bash Lags','No Bash'})
    legend({'Target Leads (Masker Lags)','Masker Leads (Target Lags)'})
    ylabel('P300 (uV)')
    ylim([-2,7])

    subplot(1,2,2)
    hold on
    errorbar(0.9:1:3.1,this_mean(4:6),this_sem(4:6),'k','LineStyle','none','Marker','o','MarkerEdgeColor','g','MarkerFaceColor','g'); % Target Leads (Masker Lags)
    errorbar(1.1:1:3.1,this_mean(10:12),this_sem(10:12),'k','LineStyle','none','Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b'); % Masker Leads (Target Lags)
    title('Lag Word Response')
    xticks(1:3)
    xticklabels({'Bash Leads','Bash Lags','No Bash'})
    legend({'Target Leads (Masker Lags)','Masker Leads (Target Lags)'})
    ylabel('P300 (uV)')
    ylim([-2,7])

    sgtitle(append(condition_titles(icondition),' P300 Height'))
end


%% ATTEND LEFT VS. ATTEND RIGHT


% TOPOPLOTS P1
eeglab;
EEG_struct_for_topographies = load(append(mild_master_root,'prepro_epoched_data\epochs_for_topographies.mat'));
EEG_struct_for_topographies = EEG_struct_for_topographies.EEG_yes_button_press;
cmin = -4;
cmax = 4;

for icondition = 1:4
    figure;
    subplot(2,2,1)
    this_data = squeeze(mean(all_subs_p1_target_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Left Onset")

    subplot(2,2,2)
    this_data = squeeze(mean(all_subs_p1_target_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Right Onset")

    subplot(2,2,3)
    this_data = squeeze(mean(all_subs_p1_masker_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Left Onset")

    subplot(2,2,4)
    this_data = squeeze(mean(all_subs_p1_masker_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Right Onset")


    this_colorbar=colorbar;
    this_colorbar.Position = [0.93 0.168 0.022 0.7];

    sgtitle(append(condition_titles(icondition), ' P1 topoplots'),'FontSize',18)
end

% TOPOPLOTS N1
cmin = -2;
cmax = 2;
for icondition = 1:4
    figure;
    subplot(2,2,1)
    this_data = squeeze(mean(all_subs_n1_target_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Left Onset")

    subplot(2,2,2)
    this_data = squeeze(mean(all_subs_n1_target_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Right Onset")

    subplot(2,2,3)
    this_data = squeeze(mean(all_subs_n1_masker_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Left Onset")

    subplot(2,2,4)
    this_data = squeeze(mean(all_subs_n1_masker_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Right Onset")


    this_colorbar=colorbar;
    this_colorbar.Position = [0.93 0.168 0.022 0.7];

    sgtitle(append(condition_titles(icondition), ' N1 topoplots'),'FontSize',18)
end

% TOPOPLOTS P3
cmin = 0;
cmax = 4;
for icondition = 1:4
    figure;
    subplot(2,2,1)
    this_data = squeeze(mean(all_subs_p3_target_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Left Onset")

    subplot(2,2,2)
    this_data = squeeze(mean(all_subs_p3_target_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Right Onset")

    subplot(2,2,3)
    this_data = squeeze(mean(all_subs_p3_masker_left_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Right Left Onset")

    subplot(2,2,4)
    this_data = squeeze(mean(all_subs_p3_masker_right_onset(:,icondition,:,:),[2,4]));
    this_mean = mean(this_data,1);
    this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
    topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
    title("Attend Left Right Onset")


    this_colorbar=colorbar;
    this_colorbar.Position = [0.93 0.168 0.022 0.7];

    sgtitle(append(condition_titles(icondition), ' P3 topoplots'),'FontSize',18)
end
