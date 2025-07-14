%% mildmaster_n1_lateralization

% Create N1 lateralization topoplots and calculations for mild master
% experiment.

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

all_subs_p1 = [];
all_subs_n1 = [];
all_subs_p3 = [];

all_subs_lead_p1 = [];
all_subs_lag_p1 = [];
all_subs_lead_n1 = [];
all_subs_lag_n1 = [];

all_subs_lead_p1_attend_left = [];
all_subs_lead_p1_attend_right = [];
all_subs_lag_p1_attend_left = [];
all_subs_lag_p1_attend_right = [];
all_subs_lead_n1_attend_left = [];
all_subs_lead_n1_attend_right = [];
all_subs_lag_n1_attend_left = [];
all_subs_lag_n1_attend_right = [];

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
    [~,lead_p1_start_index] = min(abs(single_onset_time - (75)));
    [~,lead_p1_end_index] = min(abs(single_onset_time - (125)));
    [~,lead_n1_start_index] = min(abs(single_onset_time - (125)));
    [~,lead_n1_end_index] = min(abs(single_onset_time - (175)));
    [~,lead_p3_start_index] = min(abs(single_onset_time - (450)));
    [~,lead_p3_end_index] = min(abs(single_onset_time - (550)));

    % Lag P1: local maximum betweeen 300 and 400 ms
    % Lag N1: local minimum between 350 and 450 ms
    % Lag P3: local maximum between 650 and 850 ms
    [~,lag_p1_start_index] = min(abs(single_onset_time - (325)));
    [~,lag_p1_end_index] = min(abs(single_onset_time - (375)));
    [~,lag_n1_start_index] = min(abs(single_onset_time - (375)));
    [~,lag_n1_end_index] = min(abs(single_onset_time - (425)));
    [~,lag_p3_start_index] = min(abs(single_onset_time - (700)));
    [~,lag_p3_end_index] = min(abs(single_onset_time - (800)));

    % take the average
    this_sub_cz_fz_average = squeeze(mean(data_by_pair_onset_baselined([frontocentral_channels],:,:),[1,3]));
    this_sub_lead_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lead_Word,{'bash'})),[1,3]));
    this_sub_lag_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lag_Word,{'bash'})),[1,3]));

    % lead p1
    [~,temp_lead_p1_index] = max(this_sub_cz_fz_average(lead_p1_start_index:lead_p1_end_index));
    this_sub_lead_p1_index = lead_p1_start_index + temp_lead_p1_index;
    this_sub_lead_p1_time = single_onset_time(this_sub_lead_p1_index);

    % lead n1
    [~,temp_lead_n1_index] = min(this_sub_cz_fz_average(lead_n1_start_index:lead_n1_end_index));
    this_sub_lead_n1_index = lead_n1_start_index + temp_lead_n1_index;
    this_sub_lead_n1_time = single_onset_time(this_sub_lead_n1_index);

    % lag p1
    [~,temp_lag_p1_index] = max(this_sub_cz_fz_average(lag_p1_start_index:lag_p1_end_index));
    this_sub_lag_p1_index = lag_p1_start_index + temp_lag_p1_index;
    this_sub_lag_p1_time = single_onset_time(this_sub_lag_p1_index);

    % lag n1
    [~,temp_lag_n1_index] = min(this_sub_cz_fz_average(lag_n1_start_index:lag_n1_end_index));
    this_sub_lag_n1_index = lag_n1_start_index + temp_lag_n1_index;
    this_sub_lag_n1_time = single_onset_time(this_sub_lag_n1_index);

    figure(all_subs_fig_p1n1)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_cz_fz_average,'k')
    scatter(this_sub_lead_p1_time,this_sub_cz_fz_average(this_sub_lead_p1_index),'or','filled');
    scatter(this_sub_lead_n1_time,this_sub_cz_fz_average(this_sub_lead_n1_index),'ob','filled');
    scatter(this_sub_lag_p1_time,this_sub_cz_fz_average(this_sub_lag_p1_index),'og','filled');
    scatter(this_sub_lag_n1_time,this_sub_cz_fz_average(this_sub_lag_n1_index),'om','filled');

    all_subs_lead_p1_times(isubject,:) = this_sub_lead_p1_time;
    all_subs_lag_p1_times(isubject,:) = this_sub_lag_p1_time;

    all_subs_lead_n1_times(isubject,:) = this_sub_lead_n1_time;
    all_subs_lag_n1_times(isubject,:) = this_sub_lag_n1_time;

    all_subs_lead_p1_indices(isubject,:) = this_sub_lead_p1_index;
    all_subs_lag_p1_indices(isubject,:) = this_sub_lag_p1_index;

    all_subs_lead_n1_indices(isubject,:) = this_sub_lead_n1_index;
    all_subs_lag_n1_indices(isubject,:) = this_sub_lag_n1_index;

    % Build data arrays to plot later
    % Num subjects x num conditions x num situations x num channels (32)
    % Condition order: small ITDs, large ITDs, small ILDs, large ILDs
    % ORGANIZATION:
    % (1) lead response, bash leads, target leads, attend right
    % (2) lead response, bash lags, target leads, attend right
    % (3) lead response, no bash, target leads, attend right
    % (4) lag response, bash leads, target leads, attend right
    % (5) lag response, bash lags, target leads, attend right
    % (6) lag response, no bash, target leads, attend right
    % (7) lead response, bash leads, masker leads, attend right
    % (8) lead response, bash lags, masker leads, attend right
    % (9) lead response, no bash, masker leads, attend right
    % (10) lag response, bash leads, masker leads, attend right
    % (11) lag response, bash lags, masker leads, attend right
    % (12) lag response, no bash, masker leads, attend right

    % (13) lead response, bash leads, target leads, attend left
    % (14) lead response, bash lags, target leads, attend left
    % (15) lead response, no bash, target leads, attend left
    % (16) lag response, bash leads, target leads, attend left
    % (17) lag response, bash lags, target leads, attend left
    % (18) lag response, no bash, target leads, attend left
    % (19) lead response, bash leads, masker leads, attend left
    % (20) lead response, bash lags, masker leads, attend left
    % (21) lead response, no bash, masker leads, attend left
    % (22) lag response, bash leads, masker leads, attend left
    % (23) lag response, bash lags, masker leads, attend left
    % (24) lag response, no bash, masker leads, attend left

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
        all_subs_p1(isubject,icondition,1,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,2,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,3,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,4,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,5,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,6,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));

        all_subs_p1(isubject,icondition,7,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,8,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,9,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,10,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,11,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_p1(isubject,icondition,12,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));

        all_subs_p1(isubject,icondition,13,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,14,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,15,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,16,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,17,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,18,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));

        all_subs_p1(isubject,icondition,19,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,20,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,21,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,22,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,23,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_p1(isubject,icondition,24,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));




    end

    % n1
    for icondition = 1:4
        these_conditions = all_conditions(icondition,:);
        all_subs_n1(isubject,icondition,1,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,2,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,3,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,4,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,5,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,6,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));

        all_subs_n1(isubject,icondition,7,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,8,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,9,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,10,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,11,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));
        all_subs_n1(isubject,icondition,12,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'R'}))),[2,3]));

        all_subs_n1(isubject,icondition,13,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,14,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,15,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,16,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,17,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,18,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));

        all_subs_n1(isubject,icondition,19,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,20,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,21,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,22,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,23,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));
        all_subs_n1(isubject,icondition,24,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}).*ismember(ERP_info(:).Target_Direction,{'L'}))),[2,3]));


    end

    all_subs_lead_p1(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,[small_itd_cond,large_itd_cond,small_ild_cond,large_ild_cond])')),[2,3]));
    all_subs_lag_p1(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,[small_itd_cond,large_itd_cond,small_ild_cond,large_ild_cond])')),[2,3]));

    all_subs_lead_n1(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,[small_itd_cond,large_itd_cond,small_ild_cond,large_ild_cond])')),[2,3]));
    all_subs_lag_n1(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,[small_itd_cond,large_itd_cond,small_ild_cond,large_ild_cond])')),[2,3]));



    all_subs_lead_p1_attend_left(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'L'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Target'}))),[2,3]));
    all_subs_lead_p1_attend_right(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'R'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Target'}))),[2,3]));
    all_subs_lag_p1_attend_left(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'L'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}))),[2,3]));
    all_subs_lag_p1_attend_right(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'R'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}))),[2,3]));


    all_subs_lead_n1_attend_left(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'L'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Target'}))),[2,3]));
    all_subs_lead_n1_attend_right(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'R'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Target'}))),[2,3]));
    all_subs_lag_n1_attend_left(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'L'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}))),[2,3]));
    all_subs_lag_n1_attend_right(isubject,:) = squeeze(nanmean(data_by_pair_onset_baselined(:,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Target_Direction,{'R'}).*ismember(ERP_info(:).Condition,[large_itd_cond])'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}))),[2,3]));



end

%% ATTEND LEFT Vs. ATTEND RIGHT
eeglab;
EEG_struct_for_topographies = load(append(mild_master_root,'prepro_epoched_data\epochs_for_topographies.mat'));
EEG_struct_for_topographies = EEG_struct_for_topographies.EEG_yes_button_press;

cmin = -3;
cmax = 3;

figure;
subplot(1,4,1)
topoplot(squeeze(mean(all_subs_lead_n1_attend_left,1)), EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
title('Lead N1 Attend Left')
subplot(1,4,2)
topoplot(squeeze(mean(all_subs_lead_n1_attend_right,1)), EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
title('Lead N1 Attend Right')
subplot(1,4,3)
topoplot(squeeze(mean(all_subs_lag_n1_attend_left,1)), EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
title('Lag N1 Attend Left')
subplot(1,4,4)
topoplot(squeeze(mean(all_subs_lag_n1_attend_right,1)), EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
title('Lag N1 Attend Right')
this_colorbar=colorbar;
    this_colorbar.Position = [0.93 0.168 0.022 0.7];

%% ALL CONDITIONS AND SITUATIONS
condition_titles = {'Small ITD','Large ITD','Small ILD','Large ILD'};
situation_titles = {'Bash Leads, Lead, Target, R','Bash Lags, Lead, Target, R','No Bash, Lead, Target, R',...
                    'Bash Leads, Lag, Masker, R','Bash Lags, Lag, Masker, R','No Bash, Lag, Masker, R',...
                    'Bash Leads, Lead, Masker, R','Bash Lags, Lead, Masker, R','No Bash, Lead, Masker, R',...
                    'Bash Leads, Lag, Target, R','Bash Lags, Lag, Target, R','No Bash, Lag, Target, R',...
                    'Bash Leads, Lead, Target, L','Bash Lags, Lead, Target, L','No Bash, Lead, Target, L',...
                    'Bash Leads, Lag, Masker, L','Bash Lags, Lag, Masker, L','No Bash, Lag, Masker, L',...
                    'Bash Leads, Lead, Masker, L','Bash Lags, Lead, Masker, L','No Bash, Lead, Masker, L',...
                    'Bash Leads, Lag, Target, L','Bash Lags, Lag, Target, L','No Bash, Lag, Target, L'};

% TOPOPLOTS P1


% TOPOPLOTS N1
cmin = -3;
cmax = 3;
for icondition = 1:4
    figure;
    for isituation = 1:24
        subplot(4,6,isituation)
        this_data = squeeze(all_subs_n1(:,icondition,isituation,:)); % -  squeeze(mean(all_subs_p1,[2,3]));
        this_mean = mean(this_data,1);
        this_sem = std(this_data,[],1)/(sqrt(num_subjects)-1);
        topoplot(this_mean, EEG_struct_for_topographies.chanlocs, 'maplimits',[cmin,cmax]);
        title(situation_titles(isituation))
    end
    this_colorbar=colorbar;
    this_colorbar.Position = [0.93 0.168 0.022 0.7];

end