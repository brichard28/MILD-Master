%% mildmaster_eeg_plotting.m

% Author: Benjamin Richardson
% 10/14/2024

curr_subject_ID =  char('fullpilot1','fullpilot2','fullpilot3');
itd50_by_lead_target_onset = [];
itd500_by_lead_target_onset = [];
ild5_by_lead_target_onset = [];
ild5mag_by_lead_target_onset = [];
itd50_by_lag_target_onset = [];
itd500_by_lag_target_onset = [];
ild5_by_lag_target_onset = [];
ild5mag_by_lag_target_onset = [];


for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',string(curr_subject_ID(isubject,:)),'.mat'))

    % Remove noisy ERPs
    
    % Plotting parameters
    erp_window_start_time = -100;
    erp_window_end_time = 750;
    button_press_delay = -500;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    cz_index = 32;
    

    % sort lead data into conditions
    itd50_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    itd500_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild5mag_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));

    % sort lag data into conditions
    itd50_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    itd500_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild5mag_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));

end

curr_channel_index = cz_index;

ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ITD50
hold on
this_lead_data = squeeze(mean(itd50_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(itd50_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time - 250])
ylabel('Voltage (uV)','FontSize',18)
title('50 us ITD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)

subplot(1,4,2) % ITD500
hold on
this_lead_data = squeeze(mean(itd500_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(itd500_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time - 250])
ylabel('Voltage (uV)','FontSize',18)
title('500 us ITD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,3) % ILD5
hold on
this_lead_data = squeeze(mean(ild5_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time - 250])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time - 250])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)

