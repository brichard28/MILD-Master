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
<<<<<<< HEAD
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40', ...
    'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_46','mild_master_48'); % char();
num_subjects = size(curr_subject_ID,1);
=======
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40'); % char();

>>>>>>> parent of 70d79e5 (yur)

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

<<<<<<< HEAD
erp_window_start_time = -25; % 100 ms before onset of word
erp_window_end_time = 1200; % 750 ms after onset of word
=======
erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 1500; % 750 ms after onset of word
>>>>>>> parent of 70d79e5 (yur)

small_itd_cond = [3,7];
large_itd_cond = [6,8];
small_ild_cond = [1,4];
large_ild_cond = [2,5];

%% WITHOUT BUTTON PRESS SUBTRACTION

% for isubject = 1:size(curr_subject_ID,1)
%     subID = string(curr_subject_ID(isubject,:));
%     disp(subID)
%     % Load Data
%     load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_no_button_press.mat'))
% 
%     %% Plot all channels, remove noisy ones in time domain
% 
%     button_press_delay =0 ;
%     single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
%     single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
%     frontocentral_channels = [1,2,3,4,5,6,8,23,25,26,27,28,29,30,31,32];
%     parietooccipital_channels = 11:20;
%     cz_index = 32;
% 
%     %% Break up data by position and word type
%     % ITD5
%     itd5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     itd5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     itd5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     itd5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     itd5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     itd5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
% 
%     itd15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     itd15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     itd15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     itd15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     itd15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     itd15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_itd_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     ild5_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     ild5_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     ild5_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     ild5_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     ild5_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     ild5_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,small_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     ild15_lead_bash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     ild15_lead_nonbash_target_lag_bash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     ild15_lead_nonbash_target_lag_nonbash_masker(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Target'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     ild15_lead_bash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'bash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
%     ild15_lead_nonbash_masker_lag_bash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'bash'}))),3));
%     ild15_lead_nonbash_masker_lag_nonbash_target(isubject,:,:)  = squeeze(nanmean(data_by_pair_onset_baselined(:,:,logical(ismember(ERP_info(:).Condition,large_ild_cond)'.*ismember(ERP_info(:).Lead_Stream,{'Masker'}).*ismember(ERP_info(:).Lead_Word,{'dash','gash'}).*ismember(ERP_info(:).Lag_Word,{'dash','gash'}))),3));
% 
%     %% Button Press
% 
%     all_subjects_button_press(isubject,:,:) = squeeze(nanmean(data_by_button_press_baselined,3));
% end

<<<<<<< HEAD
% 
% curr_channel_indices = [31, 32, 13]; % Fz, Cz, and Pz
% for curr_channel_index = curr_channel_indices
%     figure;
%     ymin = -3;
%     ymax = 4.5;
%     % 2 x 6 subplots
%     % TOP ROW: Small ITD Bash Lead, Small ITD Bash Lag, Small ITD No Bash, Small ILD Bash Lead, Small ILD Bash Lag, Small ILD No Bash
%     % BOTTOM ROW: Large ITD Bash Lead, Large ITD Bash Lag, Lare ITD No Bash, Large ILD Bash Lead, Large ILD Bash Lag, Large ILD No Bash
% 
%  
%     subplot(2,6,1) % Small ITD Bash Lead
%     hold on
%     this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     ylabel('Voltage (uV)','FontSize',18)
%     title('Bash Leads','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
%     subplot(2,6,2) % Small ITD Bash Lag
%     hold on
%     this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Lags','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
%     subplot(2,6,3) % Small ITD No Bash
%     hold on
%     this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('No Bash','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
% 
%     subplot(2,6,4) % Small ILD Bash Lead
%     hold on
%     this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Leads','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
%     
%     subplot(2,6,5) % Small ILD Bash Lag
%     hold on
%     this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Lags','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
% 
%     subplot(2,6,6) % Small ILD No Bash
%     hold on
%     this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('No Bash','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
%     legend([p1(1).mainLine,p2(1).mainLine],{'Target','Masker'})
% 
% 
%     subplot(2,6,7) % Large ITD Bash Lead
%     hold on
%     this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     ylabel('Voltage (uV)','FontSize',18)
%     title('Bash Leads','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
% 
%     subplot(2,6,8) % Large ITD Bash Lag
%     hold on
%     this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Lags','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
%     subplot(2,6,9) % Large ITD No Bash
%     hold on
%     this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('No Bash','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
%     xlabel('Time re: lead word onset (ms)','FontSize',24)
% 
%     subplot(2,6,10) % Large ILD Bash Lead
%     hold on
%     this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Leads','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
% 
%     subplot(2,6,11) % Large ILD Bash Lag
%     hold on
%     this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('Bash Lags','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
% 
%     subplot(2,6,12) % Large ILD No Bash
%     hold on
%     this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
%     this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
%     p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
%     p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
%     ylim([ymin,ymax])
%     xlim([erp_window_start_time,erp_window_end_time])
%     title('No Bash','FontSize',18)
%     xline(0,'--k','LineWidth',2)
%     xline(250,'--k','LineWidth',2)
%     legend([p1(1).mainLine,p2(1).mainLine],{'Target','Masker'})
% 
% 
%     if curr_channel_index == 31
%         sgtitle('Fz Button Press Included','FontSize',18)
%     elseif curr_channel_index == 32
%         sgtitle('Cz Button Press Included','FontSize',18)
%     elseif curr_channel_index == 13
%         sgtitle('Pz Button Press Included','FontSize',18)
%     end
% end
=======
    %% Plot all channels, remove noisy ones in time domain

    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,3,4,5,6,8,23,25,26,27,28,29,30,31,32];
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


%% FC Bash Lead vs. Bash Lag vs. No Bash
curr_channel_index = frontocentral_channels;
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 5 deg ITDs Button Press Included','FontSize',24)



%ITD15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 15 deg ITDs Button Press Included','FontSize',24)


%ild5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 5 deg ILDs Button Press Included','FontSize',24)



%ild15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 15 deg ILDs Button Press Included','FontSize',24)




%% PO Bash Lead vs. Bash Lag vs. No Bash

curr_channel_index = parietooccipital_channels;
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 5 deg ITDs Button Press Included','FontSize',24)



%ITD15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 15 deg ITDs Button Press Included','FontSize',24)


%ild5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 5 deg ILDs Button Press Included','FontSize',24)



%ild15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 15 deg ILDs Button Press Included','FontSize',24)

>>>>>>> parent of 70d79e5 (yur)






















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

<<<<<<< HEAD

%% PLOT TARGET VS. MASKER
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
    end
    ymin = -3;
    ymax = 4.5;
    % 2 x 6 subplots
    % TOP ROW: Small ITD Bash Lead, Small ITD Bash Lag, Small ITD No Bash, Small ILD Bash Lead, Small ILD Bash Lag, Small ILD No Bash
    % BOTTOM ROW: Large ITD Bash Lead, Large ITD Bash Lag, Lare ITD No Bash, Large ILD Bash Lead, Large ILD Bash Lag, Large ILD No Bash

 
    subplot(2,6,1) % Small ITD Bash Lead
    hold on
    this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Bash Leads','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    subplot(2,6,2) % Small ITD Bash Lag
    hold on
    this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Lags','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    subplot(2,6,3) % Small ITD No Bash
    hold on
    this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('No Bash','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,6,4) % Small ILD Bash Lead
    hold on
    this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Leads','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    
    subplot(2,6,5) % Small ILD Bash Lag
    hold on
    this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Lags','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,6,6) % Small ILD No Bash
    hold on
    this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('No Bash','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    legend([p1(1).mainLine,p2(1).mainLine],{'Target','Masker'})
    xlim([xlim_min,xlim_max])


    subplot(2,6,7) % Large ITD Bash Lead
    hold on
    this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Bash Leads','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,6,8) % Large ITD Bash Lag
    hold on
    this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Lags','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    subplot(2,6,9) % Large ITD No Bash
    hold on
    this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('No Bash','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlabel('Time re: lead word onset (ms)','FontSize',24)
    xlim([xlim_min,xlim_max])

    subplot(2,6,10) % Large ILD Bash Lead
    hold on
    this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Leads','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,6,11) % Large ILD Bash Lag
    hold on
    this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('Bash Lags','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    subplot(2,6,12) % Large ILD No Bash
    hold on
    this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
    this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));
    p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
    p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    title('No Bash','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    legend([p1(1).mainLine,p2(1).mainLine],{'Target','Masker'})
    xlim([xlim_min,xlim_max])



end


%% PLOT BASH LEADS VS. BASH LAGS VS. NO BASH (Collapsed target vs. masker)
=======
%% FC Bash Lead vs. Bash Lag vs. No Bash
curr_channel_index = frontocentral_channels;
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 5 deg ITDs Button Press Subtracted','FontSize',24)



%ITD15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 15 deg ITDs Button Press Subtracted','FontSize',24)


%ild5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 5 deg ILDs Button Press Subtracted','FontSize',24)



%ild15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('FC 15 deg ILDs Button Press Subtracted','FontSize',24)




%% PO Bash Lead vs. Bash Lag vs. No Bash

curr_channel_index = parietooccipital_channels;
%ITD5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 5 deg ITDs Button Press Subtracted','FontSize',24)



%ITD15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ITD5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 15 deg ITDs Button Press Subtracted','FontSize',24)


%ild5 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 5 deg ILDs Button Press Subtracted','FontSize',24)



%ild15 
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure; % ild5
subplot(1,3,1)
hold on
this_target_data = squeeze(nanmean(ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Leads','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lead Target','Bash Lead Masker'})


subplot(1,3,2)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('Bash Lags','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'Bash Lag Target','Bash Lag Masker'})



subplot(1,3,3)
hold on
this_target_data = squeeze(nanmean(ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:),2));
this_masker_data = squeeze(nanmean(ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_target_data,1),std(this_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_masker_data,1),std(this_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('No Bash','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
xlabel('Time re: lead word onset (ms)','FontSize',18)
legend([p1(1).mainLine,p2(1).mainLine],{'No Bash Target','No Bash Masker'})


sgtitle('PO 15 deg ILDs Button Press Subtracted','FontSize',24)
>>>>>>> parent of 70d79e5 (yur)


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
        xlim_min = -25;
        xlim_max = 450;
    elseif sum(ismember(curr_channel_index,parietooccipital_channels)) == length(parietooccipital_channels)
        xlim_min = -25;
        xlim_max = 1000;
    end

    ymin = -3;
    ymax = 4.5;
    % 2 x 6 subplots
    % TOP ROW: Small ITD Bash Lead, Small ITD Bash Lag, Small ITD No Bash, Small ILD Bash Lead, Small ILD Bash Lag, Small ILD No Bash
    % BOTTOM ROW: Large ITD Bash Lead, Large ITD Bash Lag, Lare ITD No Bash, Large ILD Bash Lead, Large ILD Bash Lag, Large ILD No Bash

 
    subplot(2,2,1) % Small ITD
    hold on
    this_bash_leads_data = squeeze(nanmean([itd5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:);itd5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:)],2));
    this_bash_lags_data = squeeze(nanmean([itd5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:);itd5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:)],2));
    this_no_bash_data = squeeze(nanmean([itd5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:);itd5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:)],2));
    p1 = plot(single_onset_time,nanmean(this_bash_leads_data,1),'-k');
    p2 = plot(single_onset_time,nanmean(this_bash_lags_data,1),'--k');
    p3 = plot(single_onset_time,nanmean(this_no_bash_data,1),':k');
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Small ITDs','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])


    subplot(2,2,2) % Small ILD
    hold on
    this_bash_leads_data = squeeze(nanmean([ild5_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:);ild5_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:)],2));
    this_bash_lags_data = squeeze(nanmean([ild5_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:);ild5_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:)],2));
    this_no_bash_data = squeeze(nanmean([ild5_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:);ild5_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:)],2));
   
     p1 = plot(single_onset_time,nanmean(this_bash_leads_data,1),'-k');
    p2 = plot(single_onset_time,nanmean(this_bash_lags_data,1),'--k');
    p3 = plot(single_onset_time,nanmean(this_no_bash_data,1),':k');
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Small ILDs','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    legend({'Bash Leads','Bash Lags','No Bash'})
    xlim([xlim_min,xlim_max])

    subplot(2,2,3) % Large ITD
    hold on
    this_bash_leads_data = squeeze(nanmean([itd15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:);itd15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:)],2));
    this_bash_lags_data = squeeze(nanmean([itd15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:);itd15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:)],2));
    this_no_bash_data = squeeze(nanmean([itd15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:);itd15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:)],2));
   
     p1 = plot(single_onset_time,nanmean(this_bash_leads_data,1),'-k');
    p2 = plot(single_onset_time,nanmean(this_bash_lags_data,1),'--k');
    p3 = plot(single_onset_time,nanmean(this_no_bash_data,1),':k');
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Large ITDs','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

    subplot(2,2,4) % Large ILD
    hold on
    this_bash_leads_data = squeeze(nanmean([ild15_lead_bash_target_lag_nonbash_masker(:,curr_channel_index,:);ild15_lead_bash_masker_lag_nonbash_target(:,curr_channel_index,:)],2));
    this_bash_lags_data = squeeze(nanmean([ild15_lead_nonbash_masker_lag_bash_target(:,curr_channel_index,:);ild15_lead_nonbash_target_lag_bash_masker(:,curr_channel_index,:)],2));
    this_no_bash_data = squeeze(nanmean([ild15_lead_nonbash_masker_lag_nonbash_target(:,curr_channel_index,:);ild15_lead_nonbash_target_lag_nonbash_masker(:,curr_channel_index,:)],2));
    p1 = plot(single_onset_time,nanmean(this_bash_leads_data,1),'-k');
    p2 = plot(single_onset_time,nanmean(this_bash_lags_data,1),'--k');
    p3 = plot(single_onset_time,nanmean(this_no_bash_data,1),':k');
    ylim([ymin,ymax])
    xlim([erp_window_start_time,erp_window_end_time])
    ylabel('Voltage (uV)','FontSize',18)
    title('Large ILDs','FontSize',18)
    xline(0,'--k','LineWidth',2)
    xline(250,'--k','LineWidth',2)
    xlim([xlim_min,xlim_max])

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
% plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'-k')
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
% plot(single_onset_time_buttonpress,squeeze(nanmean(all_subjects_button_press(:,curr_channel_index,:),2)),'-k')
% p1 = shadedErrorBar(single_onset_time_buttonpress,nanmean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('PO Button Press During Experiment','FontSize',18)

