%% mildmaster_eeg_plotting.m

% Author: Benjamin Richardson
% 10/14/2024

% curr_subject_ID = char('mild_master_1',...
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
%     'mild_master_26',...
%     'mild_master_27',...
%     'mild_master_28',...
%     'mild_master_30',...
%     'mild_master_31',...
%     'mild_master_32'); % char();
curr_subject_ID = char('mild_master_1',...
    'mild_master_3',...
    'mild_master_4',...
    'mild_master_5',...
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
    'mild_master_23',...
    'mild_master_24',...
    'mild_master_25',...
    'mild_master_26',...
    'mild_master_27',...
    'mild_master_28',...
    'mild_master_29',...
    'mild_master_30',...
    'mild_master_31',...
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40'); % char();


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
erp_window_end_time = 2000; % 750 ms after onset of word

small_itd_cond = [3,7];
large_itd_cond = [6,8];
small_ild_cond = [1,4];
large_ild_cond = [2,5];

%% WITHOUT BUTTON PRESS SUBTRACTION

for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_no_button_press.mat'))

    %% Plot all channels, remove noisy ones in time domain

    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    parietooccipital_channels = 11:20;
    cz_index = 32;

    %% Break up data by position and word type
    % ITD5
    itd5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));


    itd15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    %% Button Press

    all_subjects_button_press(isubject,:,:) = squeeze(nanmean(data_by_button_press_baselined,3));
end

curr_channel_index = frontocentral_channels;

%% Comparing in first position
% Bash in target vs. Bash in masker vs. No bash at all
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 5 deg ITDs Button Press Included','FontSize',24)

% ITD15
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 15 deg ITDs Button Press Included','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 5 deg ILDs Button Press Included','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild15
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 15 deg ILDs Button Press Included','FontSize',24)

%% P300 plot: BASH in first pos. target minus BASH in first pos. non-target
curr_channel_index = parietooccipital_channels;

% Bash in target vs. Bash in masker vs. No bash at all
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 5 deg ITDs Button Press Included','FontSize',24)

% ITD15
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 15 deg ITDs Button Press Included','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 5 deg ILDs Button Press Included','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild15
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 15 deg ILDs Button Press Included','FontSize',24)


















%% PLOT WITH BUTTON PRESS SUBTRACTION
for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_yes_button_press.mat'))

    %% Plot all channels, remove noisy ones in time domain

    % Plotting parameters
    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    parietooccipital_channels = 11:20;
    cz_index = 32;

    %% Break up data by position and word type
    % ITD5
    itd5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));


    itd15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    itd15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    itd15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    itd15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    ild15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
    ild15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
    ild15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));

    %% Button Press

    all_subjects_button_press(isubject,:,:) = squeeze(nanmean(data_by_button_press_baselined,3));
end

curr_channel_index = frontocentral_channels;

%% Comparing in first position
% Bash in target vs. Bash in masker vs. No bash at all
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 5 deg ITDs Button Press Subtracted','FontSize',24)

% ITD15
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 15 deg ITDs Button Press Subtracted','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 5 deg ILDs Button Press Subtracted','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild15
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('FC 15 deg ILDs Button Press Subtracted','FontSize',24)

%% P300 plot: BASH in first pos. target minus BASH in first pos. non-target
curr_channel_index = parietooccipital_channels;

% Bash in target vs. Bash in masker vs. No bash at all
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 5 deg ITDs Button Press Subtracted','FontSize',24)

% ITD15
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 15 deg ITDs Button Press Subtracted','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 5 deg ILDs Button Press Subtracted','FontSize',24)


ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild15
subplot(1,2,1)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Target Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})

subplot(1,2,2)
hold on
this_lead_bash_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
this_lag_bash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));


p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_bash_data,1),std(this_lead_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_bash_data,1),std(this_lag_bash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
title('Masker Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash Leads','Bash Lags','No Bash'})
xlabel('Time re: lead word onset (ms)','FontSize',18)

sgtitle('PO 15 deg ILDs Button Press Subtracted','FontSize',24)














%% Plot all button press during experiment
button_press_delay = 0;
single_onset_time_buttonpress = linspace(erp_window_start_time,erp_window_end_time,size(all_subjects_button_press,3));
curr_channel_index=frontocentral_channels;
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
hold on
this_data = squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2));
plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'-k')
p1 = shadedErrorBar(single_onset_time_buttonpress,nanmean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('FC Button Press During Experiment','FontSize',18)


single_onset_time_buttonpress = linspace(erp_window_start_time,erp_window_end_time,size(all_subjects_button_press,3));
curr_channel_index=parietooccipital_channels;
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
hold on
this_data = squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2));
plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'-k')
p1 = shadedErrorBar(single_onset_time_buttonpress,nanmean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('PO Button Press During Experiment','FontSize',18)

