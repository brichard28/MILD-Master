%% mildmaster_eeg_plotting.m

% Author: Benjamin Richardson
% 10/14/2024

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
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40',...
    'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_46','mild_master_48'); % char();
num_subjects = size(curr_subject_ID,1);
mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master\';

itd5_lead_bash_target_lag_nonbash_masker = [];
itd5_lead_nonbash_target_lag_bash_masker = [];
itd5_lead_nonbash_target_lag_nonbash_masker= [];

itd5_lead_bash_masker_lag_nonbash_target = [];
itd5_lead_nonbash_masker_lag_bash_target = [];
itd5_lead_nonbash_masker_lag_nonbash_target= [];

itd15_lead_bash_target_lag_nonbash_masker = [];
itd15_lead_nonbash_target_lag_bash_masker = [];
itd15_lead_nonbash_target_lag_nonbash_masker= [];

itd15_lead_bash_masker_lag_nonbash_target = [];
itd15_lead_nonbash_masker_lag_bash_target = [];
itd15_lead_nonbash_masker_lag_nonbash_target= [];


ild5_lead_bash_target_lag_nonbash_masker = [];
ild5_lead_nonbash_target_lag_bash_masker = [];
ild5_lead_nonbash_target_lag_nonbash_masker= [];

ild5_lead_bash_masker_lag_nonbash_target = [];
ild5_lead_nonbash_masker_lag_bash_target = [];
ild5_lead_nonbash_masker_lag_nonbash_target= [];

ild15_lead_bash_target_lag_nonbash_masker = [];
ild15_lead_nonbash_target_lag_bash_masker = [];
ild15_lead_nonbash_target_lag_nonbash_masker= [];

ild15_lead_bash_masker_lag_nonbash_target = [];
ild15_lead_nonbash_masker_lag_bash_target = [];
ild15_lead_nonbash_masker_lag_nonbash_target= [];

erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 1500; % 750 ms after onset of word

small_itd_cond = [3,7];
large_itd_cond = [6,8];
small_ild_cond = [1,4];
large_ild_cond = [2,5];

all_subs_p1 = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p1.csv');
all_subs_n1 = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_n1.csv');
all_subs_p2 = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p2.csv');
all_subs_p3 = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p3.csv');




%% PLOT WITH BUTTON PRESS SUBTRACTION
for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_no_button_press.mat'))

    %% Plot all channels, remove noisy ones in time domain

    % Plotting parameters
    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [31, 5, 26, 8, 32, 23, 9, 22]; % Fz, FC1, FC2, C3, Cz, C4, CP1, and CP2
    frontocentral_channel_names = {'Fz', 'FC1', 'FC2', 'C3', 'Cz', 'C4', 'CP1', 'CP2'};
    parietooccipital_channels = [12, 13, 14, 15, 16, 17, 18, 19] ;%  P3, Pz, PO3, O1, Oz, O2, PO4, and P4
    parietooccipital_channel_names = {'P3', 'Pz', 'PO3', 'O1', 'Oz', 'O2', 'PO4', 'P4'};
    cz_index = 32;

    %% Break up data by position and word type
    % ITD5
    itd5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));


    itd15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,1).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Responded,0).*ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));


    %% Button Press

    all_subjects_button_press(isubject,:,:) = squeeze(nanmean(data_by_button_press_baselined,3));
end


%% PLOT TIME TRACES: Target (solid line) vs. Masker (dashed line)
% Color: Bash lead (blue), bash lag (orange), no bash (green)



colors = [0, 0.4470, 0.7410;   % Blue
    0.8500, 0.3250, 0.0980;   % Orange
    0.4660, 0.6740, 0.1880];  % Green
curr_channel_indices = {[frontocentral_channels],[parietooccipital_channels]}; % Fz, Cz, and Pz
for i = 1:length(curr_channel_indices)
    curr_channel_index = cell2mat(curr_channel_indices(i));


    figure;
    if curr_channel_index == 31
        sgtitle('Fz Button Press Subtracted','FontSize',18)
        xlim_min = -25;
        xlim_max = 1000;
    elseif curr_channel_index == 32
        sgtitle('Cz Button Press Subtracted','FontSize',18)
        xlim_min = -25;
        xlim_max = 1000;
    elseif curr_channel_index == 13
        sgtitle('Pz Button Press Subtracted','FontSize',18)
        xlim_min = -25;
        xlim_max = 1000;
    elseif sum(ismember(curr_channel_index,frontocentral_channels)) == length(frontocentral_channels)
        sgtitle('Frontocentral Channels','FontSize',18)
        xlim_min = -25;
        xlim_max = 450;
    elseif sum(ismember(curr_channel_index,parietooccipital_channels)) == length(parietooccipital_channels)
        sgtitle('Parietooccipital Channels','FontSize',18)
        xlim_min = -25;
        xlim_max = 1000;
    elseif sum(ismember(curr_channel_index,[31,32])) == 2
        xlim_min = -25;
        xlim_max = 450;
    end
    ymin = -3;
    ymax = 4.5;

    subplot(2,2,1) % Small ITD Time Trace
    hold on
    this_target_lead_bash_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_lead_bash_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    this_target_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

    this_target_no_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_no_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    plot1 = plot(single_onset_time,nanmean(this_target_lead_bash_data,1),'Color',colors(1,:),'LineStyle','-','LineWidth',1.5);
    plot2 = plot(single_onset_time,nanmean(this_masker_lead_bash_data,1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
    plot3 = plot(single_onset_time,nanmean(this_target_lag_bash_data,1),'Color',colors(2,:),'LineStyle','-','LineWidth',1.5);
    plot4 = plot(single_onset_time,nanmean(this_masker_lag_bash_data,1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
    plot5 = plot(single_onset_time,nanmean(this_target_no_bash_data,1),'Color',colors(3,:),'LineStyle','-','LineWidth',1.5);
    plot6 = plot(single_onset_time,nanmean(this_masker_no_bash_data,1),'Color',colors(3,:),'LineStyle','--','LineWidth',1.5);
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Small ITD','FontSize',18)
    xline(0,':k','LineWidth',2)
    xline(250,':k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,2,2) % Large ITD Time Trace
    hold on
    this_target_lead_bash_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_lead_bash_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    this_target_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

    this_target_no_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_no_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    plot1 = plot(single_onset_time,nanmean(this_target_lead_bash_data,1),'Color',colors(1,:),'LineStyle','-','LineWidth',1.5);
    plot2 = plot(single_onset_time,nanmean(this_masker_lead_bash_data,1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
    plot3 = plot(single_onset_time,nanmean(this_target_lag_bash_data,1),'Color',colors(2,:),'LineStyle','-','LineWidth',1.5);
    plot4 = plot(single_onset_time,nanmean(this_masker_lag_bash_data,1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
    plot5 = plot(single_onset_time,nanmean(this_target_no_bash_data,1),'Color',colors(3,:),'LineStyle','-','LineWidth',1.5);
    plot6 = plot(single_onset_time,nanmean(this_masker_no_bash_data,1),'Color',colors(3,:),'LineStyle','--','LineWidth',1.5);
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Large ITD','FontSize',18)
    xline(0,':k','LineWidth',2)
    xline(250,':k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,2,3) % Small ILD Time Trace
    hold on
    this_target_lead_bash_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_lead_bash_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    this_target_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

    this_target_no_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_no_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    plot1 = plot(single_onset_time,nanmean(this_target_lead_bash_data,1),'Color',colors(1,:),'LineStyle','-','LineWidth',1.5);
    plot2 = plot(single_onset_time,nanmean(this_masker_lead_bash_data,1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
    plot3 = plot(single_onset_time,nanmean(this_target_lag_bash_data,1),'Color',colors(2,:),'LineStyle','-','LineWidth',1.5);
    plot4 = plot(single_onset_time,nanmean(this_masker_lag_bash_data,1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
    plot5 = plot(single_onset_time,nanmean(this_target_no_bash_data,1),'Color',colors(3,:),'LineStyle','-','LineWidth',1.5);
    plot6 = plot(single_onset_time,nanmean(this_masker_no_bash_data,1),'Color',colors(3,:),'LineStyle','--','LineWidth',1.5);
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Small ILD','FontSize',18)
    xline(0,':k','LineWidth',2)
    xline(250,':k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,2,4) % Large ILD Time Trace
    hold on
    this_target_lead_bash_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_lead_bash_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    this_target_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

    this_target_no_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_no_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

    plot1 = plot(single_onset_time,nanmean(this_target_lead_bash_data,1),'Color',colors(1,:),'LineStyle','-','LineWidth',1.5);
    plot2 = plot(single_onset_time,nanmean(this_masker_lead_bash_data,1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
    plot3 = plot(single_onset_time,nanmean(this_target_lag_bash_data,1),'Color',colors(2,:),'LineStyle','-','LineWidth',1.5);
    plot4 = plot(single_onset_time,nanmean(this_masker_lag_bash_data,1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
    plot5 = plot(single_onset_time,nanmean(this_target_no_bash_data,1),'Color',colors(3,:),'LineStyle','-','LineWidth',1.5);
    plot6 = plot(single_onset_time,nanmean(this_masker_no_bash_data,1),'Color',colors(3,:),'LineStyle','--','LineWidth',1.5);
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Large ILD','FontSize',18)
    xline(0,':k','LineWidth',2)
    xline(250,':k','LineWidth',2)
    xlim([xlim_min,xlim_max])




end





%% PLOT ERP Calculations
% Left column: Lead word response
% Right column: Lag word response
for i = 1:length(curr_channel_indices)
    curr_channel_index = cell2mat(curr_channel_indices(i));

    figure;
    hold on
    if sum(ismember(curr_channel_index,frontocentral_channels)) == length(frontocentral_channels) % plot p1-n1
        ymin = -1.5;
        ymax = 0.6;
        % Small ITD Measure

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_bash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_bash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        lead_p1n1_to_plot = [];
        lead_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p1n1_to_plot = [];
        lag_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,1)
        hold on
        h(1) = errorbar(1, mean(lead_p1n1_to_plot(1,:)), std(lead_p1n1_to_plot(1,:))./(sqrt(length(lead_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p1n1_to_plot(2,:)), std(lead_p1n1_to_plot(2,:))./(sqrt(length(lead_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p1n1_to_plot(3,:)), std(lead_p1n1_to_plot(3,:))./(sqrt(length(lead_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p1n1_to_plot(4,:)), std(lead_p1n1_to_plot(4,:))./(sqrt(length(lead_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p1n1_to_plot(5,:)), std(lead_p1n1_to_plot(5,:))./(sqrt(length(lead_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p1n1_to_plot(6,:)), std(lead_p1n1_to_plot(6,:))./(sqrt(length(lead_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,2)
        hold on
        h(1) = errorbar(1, mean(lag_p1n1_to_plot(1,:)), std(lag_p1n1_to_plot(1,:))./(sqrt(length(lag_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p1n1_to_plot(2,:)), std(lag_p1n1_to_plot(2,:))./(sqrt(length(lag_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p1n1_to_plot(3,:)), std(lag_p1n1_to_plot(3,:))./(sqrt(length(lag_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p1n1_to_plot(4,:)), std(lag_p1n1_to_plot(4,:))./(sqrt(length(lag_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p1n1_to_plot(5,:)), std(lag_p1n1_to_plot(5,:))./(sqrt(length(lag_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p1n1_to_plot(6,:)), std(lag_p1n1_to_plot(6,:))./(sqrt(length(lag_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])







        % Large ITD Measure

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_bash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_bash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_itd').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_itd').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        lead_p1n1_to_plot = [];
        lead_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p1n1_to_plot = [];
        lag_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,3)
        hold on
        h(1) = errorbar(1, mean(lead_p1n1_to_plot(1,:)), std(lead_p1n1_to_plot(1,:))./(sqrt(length(lead_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p1n1_to_plot(2,:)), std(lead_p1n1_to_plot(2,:))./(sqrt(length(lead_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p1n1_to_plot(3,:)), std(lead_p1n1_to_plot(3,:))./(sqrt(length(lead_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p1n1_to_plot(4,:)), std(lead_p1n1_to_plot(4,:))./(sqrt(length(lead_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p1n1_to_plot(5,:)), std(lead_p1n1_to_plot(5,:))./(sqrt(length(lead_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p1n1_to_plot(6,:)), std(lead_p1n1_to_plot(6,:))./(sqrt(length(lead_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,4)
        hold on
        h(1) = errorbar(1, mean(lag_p1n1_to_plot(1,:)), std(lag_p1n1_to_plot(1,:))./(sqrt(length(lag_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p1n1_to_plot(2,:)), std(lag_p1n1_to_plot(2,:))./(sqrt(length(lag_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p1n1_to_plot(3,:)), std(lag_p1n1_to_plot(3,:))./(sqrt(length(lag_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p1n1_to_plot(4,:)), std(lag_p1n1_to_plot(4,:))./(sqrt(length(lag_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p1n1_to_plot(5,:)), std(lag_p1n1_to_plot(5,:))./(sqrt(length(lag_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p1n1_to_plot(6,:)), std(lag_p1n1_to_plot(6,:))./(sqrt(length(lag_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])


        % Small ILD

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_bash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_bash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'small_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'small_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        lead_p1n1_to_plot = [];
        lead_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p1n1_to_plot = [];
        lag_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,5)
        hold on
        h(1) = errorbar(1, mean(lead_p1n1_to_plot(1,:)), std(lead_p1n1_to_plot(1,:))./(sqrt(length(lead_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p1n1_to_plot(2,:)), std(lead_p1n1_to_plot(2,:))./(sqrt(length(lead_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p1n1_to_plot(3,:)), std(lead_p1n1_to_plot(3,:))./(sqrt(length(lead_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p1n1_to_plot(4,:)), std(lead_p1n1_to_plot(4,:))./(sqrt(length(lead_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p1n1_to_plot(5,:)), std(lead_p1n1_to_plot(5,:))./(sqrt(length(lead_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p1n1_to_plot(6,:)), std(lead_p1n1_to_plot(6,:))./(sqrt(length(lead_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,6)
        hold on
        h(1) = errorbar(1, mean(lag_p1n1_to_plot(1,:)), std(lag_p1n1_to_plot(1,:))./(sqrt(length(lag_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p1n1_to_plot(2,:)), std(lag_p1n1_to_plot(2,:))./(sqrt(length(lag_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p1n1_to_plot(3,:)), std(lag_p1n1_to_plot(3,:))./(sqrt(length(lag_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p1n1_to_plot(4,:)), std(lag_p1n1_to_plot(4,:))./(sqrt(length(lag_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p1n1_to_plot(5,:)), std(lag_p1n1_to_plot(5,:))./(sqrt(length(lag_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p1n1_to_plot(6,:)), std(lag_p1n1_to_plot(6,:))./(sqrt(length(lag_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');
        ylim([ymin,ymax])


        % Large ILD

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_bash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_bash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_bash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Target').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_target_lag_nonbash_masker_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Target').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p1 = all_subs_p1(logical(ismember(string(all_subs_p1.Condition),'large_ild').*~ismember(string(all_subs_p1.Lead_Word),'bash').*~ismember(string(all_subs_p1.Lag_Word),'bash').*ismember(string(all_subs_p1.Lead_Stream),'Masker').*ismember(string(all_subs_p1.Electrode),frontocentral_channel_names)),:);
        lead_nonbash_masker_lag_nonbash_target_n1 = all_subs_n1(logical(ismember(string(all_subs_n1.Condition),'large_ild').*~ismember(string(all_subs_n1.Lead_Word),'bash').*~ismember(string(all_subs_n1.Lag_Word),'bash').*ismember(string(all_subs_n1.Lead_Stream),'Masker').*ismember(string(all_subs_n1.Electrode),frontocentral_channel_names)),:);


        lead_p1n1_to_plot = [];
        lead_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p1n1_to_plot = [];
        lag_p1n1_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_bash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_bash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_bash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_target_lag_nonbash_masker_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p1n1_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude - groupsummary(lead_nonbash_masker_lag_nonbash_target_n1, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,7)
        hold on
        h(1) = errorbar(1, mean(lead_p1n1_to_plot(1,:)), std(lead_p1n1_to_plot(1,:))./(sqrt(length(lead_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p1n1_to_plot(2,:)), std(lead_p1n1_to_plot(2,:))./(sqrt(length(lead_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p1n1_to_plot(3,:)), std(lead_p1n1_to_plot(3,:))./(sqrt(length(lead_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p1n1_to_plot(4,:)), std(lead_p1n1_to_plot(4,:))./(sqrt(length(lead_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p1n1_to_plot(5,:)), std(lead_p1n1_to_plot(5,:))./(sqrt(length(lead_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p1n1_to_plot(6,:)), std(lead_p1n1_to_plot(6,:))./(sqrt(length(lead_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,8)
        hold on
        h(1) = errorbar(1, mean(lag_p1n1_to_plot(1,:)), std(lag_p1n1_to_plot(1,:))./(sqrt(length(lag_p1n1_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p1n1_to_plot(2,:)), std(lag_p1n1_to_plot(2,:))./(sqrt(length(lag_p1n1_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p1n1_to_plot(3,:)), std(lag_p1n1_to_plot(3,:))./(sqrt(length(lag_p1n1_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p1n1_to_plot(4,:)), std(lag_p1n1_to_plot(4,:))./(sqrt(length(lag_p1n1_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p1n1_to_plot(5,:)), std(lag_p1n1_to_plot(5,:))./(sqrt(length(lag_p1n1_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p1n1_to_plot(6,:)), std(lag_p1n1_to_plot(6,:))./(sqrt(length(lag_p1n1_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');


        ylim([ymin,ymax])




    elseif sum(ismember(curr_channel_index,parietooccipital_channels)) == length(parietooccipital_channels) % plot p3-n1
        ymin = -1.5;
        ymax = 2;
        % Small ITD Measure

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        lead_p3_to_plot = [];
        lead_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p3_to_plot = [];
        lag_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,1)
        hold on
        h(1) = errorbar(1, mean(lead_p3_to_plot(1,:)), std(lead_p3_to_plot(1,:))./(sqrt(length(lead_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p3_to_plot(2,:)), std(lead_p3_to_plot(2,:))./(sqrt(length(lead_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p3_to_plot(3,:)), std(lead_p3_to_plot(3,:))./(sqrt(length(lead_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p3_to_plot(4,:)), std(lead_p3_to_plot(4,:))./(sqrt(length(lead_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p3_to_plot(5,:)), std(lead_p3_to_plot(5,:))./(sqrt(length(lead_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p3_to_plot(6,:)), std(lead_p3_to_plot(6,:))./(sqrt(length(lead_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,2)
        hold on
        h(1) = errorbar(1, mean(lag_p3_to_plot(1,:)), std(lag_p3_to_plot(1,:))./(sqrt(length(lag_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p3_to_plot(2,:)), std(lag_p3_to_plot(2,:))./(sqrt(length(lag_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p3_to_plot(3,:)), std(lag_p3_to_plot(3,:))./(sqrt(length(lag_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p3_to_plot(4,:)), std(lag_p3_to_plot(4,:))./(sqrt(length(lag_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p3_to_plot(5,:)), std(lag_p3_to_plot(5,:))./(sqrt(length(lag_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p3_to_plot(6,:)), std(lag_p3_to_plot(6,:))./(sqrt(length(lag_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');


        ylim([ymin,ymax])






        % Large ITD Measure

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_itd').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        lead_p3_to_plot = [];
        lead_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p3_to_plot = [];
        lag_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,3)
        hold on
        h(1) = errorbar(1, mean(lead_p3_to_plot(1,:)), std(lead_p3_to_plot(1,:))./(sqrt(length(lead_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p3_to_plot(2,:)), std(lead_p3_to_plot(2,:))./(sqrt(length(lead_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p3_to_plot(3,:)), std(lead_p3_to_plot(3,:))./(sqrt(length(lead_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p3_to_plot(4,:)), std(lead_p3_to_plot(4,:))./(sqrt(length(lead_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p3_to_plot(5,:)), std(lead_p3_to_plot(5,:))./(sqrt(length(lead_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p3_to_plot(6,:)), std(lead_p3_to_plot(6,:))./(sqrt(length(lead_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,4)
        hold on
        h(1) = errorbar(1, mean(lag_p3_to_plot(1,:)), std(lag_p3_to_plot(1,:))./(sqrt(length(lag_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p3_to_plot(2,:)), std(lag_p3_to_plot(2,:))./(sqrt(length(lag_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p3_to_plot(3,:)), std(lag_p3_to_plot(3,:))./(sqrt(length(lag_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p3_to_plot(4,:)), std(lag_p3_to_plot(4,:))./(sqrt(length(lag_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p3_to_plot(5,:)), std(lag_p3_to_plot(5,:))./(sqrt(length(lag_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p3_to_plot(6,:)), std(lag_p3_to_plot(6,:))./(sqrt(length(lag_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])


        % Small ILD

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'small_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        lead_p3_to_plot = [];
        lead_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p3_to_plot = [];
        lag_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,5)
        hold on
        h(1) = errorbar(1, mean(lead_p3_to_plot(1,:)), std(lead_p3_to_plot(1,:))./(sqrt(length(lead_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p3_to_plot(2,:)), std(lead_p3_to_plot(2,:))./(sqrt(length(lead_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p3_to_plot(3,:)), std(lead_p3_to_plot(3,:))./(sqrt(length(lead_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p3_to_plot(4,:)), std(lead_p3_to_plot(4,:))./(sqrt(length(lead_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p3_to_plot(5,:)), std(lead_p3_to_plot(5,:))./(sqrt(length(lead_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p3_to_plot(6,:)), std(lead_p3_to_plot(6,:))./(sqrt(length(lead_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,6)
        hold on
        h(1) = errorbar(1, mean(lag_p3_to_plot(1,:)), std(lag_p3_to_plot(1,:))./(sqrt(length(lag_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p3_to_plot(2,:)), std(lag_p3_to_plot(2,:))./(sqrt(length(lag_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p3_to_plot(3,:)), std(lag_p3_to_plot(3,:))./(sqrt(length(lag_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p3_to_plot(4,:)), std(lag_p3_to_plot(4,:))./(sqrt(length(lag_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p3_to_plot(5,:)), std(lag_p3_to_plot(5,:))./(sqrt(length(lag_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p3_to_plot(6,:)), std(lag_p3_to_plot(6,:))./(sqrt(length(lag_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        % Large ILD

        % LEAD RESPONSES (WordPosition = Lead)
        % Lead Bash Target, Lag Non-Bash Masker
        lead_bash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Bash Masker, Lag Non-Bash Target
        lead_bash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Target, Lag Bash Masker
        lead_nonbash_target_lag_bash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);

        % Lead Non-Bash Masker, Lag Bash Target
        lead_nonbash_masker_lag_bash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Target, Lag Non-Bash Masker
        lead_nonbash_target_lag_nonbash_masker_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Target').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        % Lead Non-Bash Masker, Lag Non-Bash Target
        lead_nonbash_masker_lag_nonbash_target_p3 = all_subs_p3(logical(ismember(string(all_subs_p3.Condition),'large_ild').*~ismember(string(all_subs_p3.Lead_Word),'bash').*~ismember(string(all_subs_p3.Lag_Word),'bash').*ismember(string(all_subs_p3.Lead_Stream),'Masker').*ismember(string(all_subs_p3.Electrode),parietooccipital_channel_names)),:);


        lead_p3_to_plot = [];
        lead_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;
        lead_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lead_Amplitude').mean_Lead_Amplitude;

        lag_p3_to_plot = [];
        lag_p3_to_plot(1,:) = groupsummary(lead_bash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(2,:) = groupsummary(lead_bash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(3,:) = groupsummary(lead_nonbash_masker_lag_bash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(4,:) = groupsummary(lead_nonbash_target_lag_bash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(5,:) = groupsummary(lead_nonbash_target_lag_nonbash_masker_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;
        lag_p3_to_plot(6,:)= groupsummary(lead_nonbash_masker_lag_nonbash_target_p3, 'S','mean','Lag_Amplitude').mean_Lag_Amplitude;



        subplot(4,2,7)
        hold on
        h(1) = errorbar(1, mean(lead_p3_to_plot(1,:)), std(lead_p3_to_plot(1,:))./(sqrt(length(lead_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lead_p3_to_plot(2,:)), std(lead_p3_to_plot(2,:))./(sqrt(length(lead_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lead_p3_to_plot(3,:)), std(lead_p3_to_plot(3,:))./(sqrt(length(lead_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lead_p3_to_plot(4,:)), std(lead_p3_to_plot(4,:))./(sqrt(length(lead_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lead_p3_to_plot(5,:)), std(lead_p3_to_plot(5,:))./(sqrt(length(lead_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lead_p3_to_plot(6,:)), std(lead_p3_to_plot(6,:))./(sqrt(length(lead_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])

        subplot(4,2,8)
        hold on
        h(1) = errorbar(1, mean(lag_p3_to_plot(1,:)), std(lag_p3_to_plot(1,:))./(sqrt(length(lag_p3_to_plot(1,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
        h(2) = errorbar(2, mean(lag_p3_to_plot(2,:)), std(lag_p3_to_plot(2,:))./(sqrt(length(lag_p3_to_plot(2,:))) - 1), ...
            'o', 'Color', colors(1,:), 'MarkerFaceColor', 'none');
        h(3) = errorbar(3, mean(lag_p3_to_plot(3,:)), std(lag_p3_to_plot(3,:))./(sqrt(length(lag_p3_to_plot(3,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
        h(4) = errorbar(4, mean(lag_p3_to_plot(4,:)), std(lag_p3_to_plot(4,:))./(sqrt(length(lag_p3_to_plot(4,:))) - 1), ...
            'o', 'Color', colors(2,:), 'MarkerFaceColor', 'none');
        h(5) = errorbar(5, mean(lag_p3_to_plot(5,:)), std(lag_p3_to_plot(5,:))./(sqrt(length(lag_p3_to_plot(5,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
        h(6) = errorbar(6, mean(lag_p3_to_plot(6,:)), std(lag_p3_to_plot(6,:))./(sqrt(length(lag_p3_to_plot(6,:))) - 1), ...
            'o', 'Color', colors(3,:), 'MarkerFaceColor', 'none');

        ylim([ymin,ymax])




    end


end

%% Plot all button press during experiment
% button_press_delay = 0;
% single_onset_time_buttonpress = linspace(erp_window_start_time,erp_window_end_time,size(all_subjects_button_press,3));
% curr_channel_index=frontocentral_channels;
% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% hold on
% this_data = squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2));
% plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'Color',colors(1,:))
% p1 = shadedErrorBar(single_onset_time_buttonpress,nanmean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('FC Button Press During Experiment','FontSize',18)
%
%
% single_onset_time_buttonpress = linspace(erp_window_start_time,erp_window_end_time,size(all_subjects_button_press,3));
% curr_channel_index=parietooccipital_channels;
% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% hold on
% this_data = squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2));
% plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'Color',colors(1,:))
% p1 = shadedErrorBar(single_onset_time_buttonpress,nanmean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('PO Button Press During Experiment','FontSize',18)

