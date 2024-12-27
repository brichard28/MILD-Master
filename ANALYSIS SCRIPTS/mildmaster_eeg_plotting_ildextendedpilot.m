%% mildmaster_eeg_plotting.m

% Author: Benjamin Richardson
% 10/14/2024

curr_subject_ID = char('eeg_pilot_2','eeg_pilot_3');%  char('fullpilot1','fullpilot2','fullpilot3'); % 
ild20_by_lead_target_onset = [];
ild20mag_by_lead_target_onset = [];
ild5_by_lead_target_onset = [];
ild5mag_by_lead_target_onset = [];
ild20_by_lag_target_onset = [];
ild20mag_by_lag_target_onset = [];
ild5_by_lag_target_onset = [];
ild5mag_by_lag_target_onset = [];


for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'))

    % Remove noisy ERPs
    
    % Plotting parameters
    erp_window_start_time = -100;
    erp_window_end_time = 750;
    button_press_delay = -500;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    cz_index = 32;
    
    conditions = {'side=r_itd=0_az=5_mag=0_lpf=0',...
'side=l_itd=0_az=20_mag=0_lpf=0',...
'side=r_itd=0_az=5_mag=1_lpf=0',...
'side=l_itd=0_az=5_mag=1_lpf=0',...
'side=r_itd=0_az=20_mag=1_lpf=0',...
'side=l_itd=0_az=20_mag=1_lpf=0',...
'side=l_itd=0_az=5_mag=0_lpf=0',...
'side=r_itd=0_az=20_mag=0_lpf=0',...
};

    %% ALL WORDS
    % sort lead data into conditions
    ild20_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,8]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild20mag_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[5,6]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,7]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild5mag_by_lead_target_onset(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,4]).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));

    % sort lag data into conditions
    ild20_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,8]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild20mag_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[5,6]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,7]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild5mag_by_lag_target_onset(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,4]).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));

    %% BASH ONLY
    % sort lead data into conditions
    ild20_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,8]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild20mag_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[5,6]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,7]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5mag_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,4]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));

    % sort lag data into conditions
    ild20_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,8]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild20mag_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[5,6]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,7]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5mag_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,4]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));



    %% {DASH,GASH} ONLY
    % sort lead data into conditions
    ild20_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,8]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild20mag_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[5,6]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,7]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,4]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));

    % sort lag data into conditions
    ild20_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,8]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild20mag_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[5,6]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,7]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,4]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));

end

curr_channel_index = frontocentral_channels; %cz_index;

%% Plot All Words
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ild20
hold on
this_lead_data = squeeze(mean(ild20_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20_by_lag_target_onset(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})

subplot(1,4,2) % ild20mag
hold on
this_lead_data = squeeze(mean(ild20mag_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20mag_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD + Mag','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,3) % ILD5
hold on
this_lead_data = squeeze(mean(ild5_by_lead_target_onset(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5_by_lag_target_onset(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
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
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)

sgtitle('All Words','FontSize',18)


%% Plot Bash Only
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ild20
hold on
this_lead_data = squeeze(mean(ild20_by_lead_target_onset_bash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20_by_lag_target_onset_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})

subplot(1,4,2) % ild20mag
hold on
this_lead_data = squeeze(mean(ild20mag_by_lead_target_onset_bash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20mag_by_lag_target_onset_bash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD + Mag','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,3) % ILD5
hold on
this_lead_data = squeeze(mean(ild5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset_bash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset_bash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)

sgtitle('BASH Only','FontSize',18)

%% Plot {Dash, Gash} Only
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ild20
hold on
this_lead_data = squeeze(mean(ild20_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})

subplot(1,4,2) % ild20mag
hold on
this_lead_data = squeeze(mean(ild20mag_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild20mag_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('20 deg ILD + Mag','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,3) % ILD5
hold on
this_lead_data = squeeze(mean(ild5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--r','LineWidth',2)
xline(250,'--b','LineWidth',2)

sgtitle('{DASH,GASH} Only','FontSize',18)