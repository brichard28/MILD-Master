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
    'mild_master_22',...
    'mild_master_23'); % char();
%char('eeg_pilot_1','eeg_pilot_4','eeg_pilot_5','eeg_pilot_6','eeg_pilot_7','eeg_pilot_8'); % char();
itd5_by_lead_target_onset = [];
itd15_by_lead_target_onset = [];
ild5_by_lead_target_onset = [];
ild15_by_lead_target_onset = [];
itd5_by_lag_target_onset = [];
itd15_by_lag_target_onset = [];
ild5_by_lag_target_onset = [];
ild15_by_lag_target_onset = [];

erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 950; % 750 ms after onset of word

    small_itd_cond = [3,7];
    large_itd_cond = [6,8];
    small_ild_cond = [1,4];
    large_ild_cond = [2,5];

for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'))

    %% Plot all channels, remove noisy ones in time domain
%     figure;
%     for ichannel = 1:32
%         subplot(6,6,ichannel)
%         %[this_spectrum,this_f] = pspectrum(squeeze(nanmean(data_by_target_onset_baselined(ichannel,:,:),2)),256);
%         %plot(this_f/pi,pow2db(this_spectrum));
%         %ylim([0,0.1]);
%         plot(linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset_baselined,2)),squeeze(nanmean(data_by_target_onset_baselined(ichannel,:,:),[1,3])))
%         title(num2str(ichannel))
%     end
%     sgtitle(subID,'FontSize',18)
    %pause
    %channels_to_remove = input('Please enter a comma-separated list of channels to remove (ex. [1,2,3]):');
    %data_by_target_onset_baselined(channels_to_remove,:,:) = nan;


    % Remove noisy ERPs
    
    % Plotting parameters
    erp_window_start_time = -50;
    erp_window_end_time = 950;
    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    parietooccipital_channels = 11:20;
    cz_index = 32;
    

    %% ALL WORDS
    % sort lead data into conditions
    itd5_by_lead_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    itd15_by_lead_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lead_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));
    ild15_by_lead_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'bash','dash','gash'})')),3));

    % sort lag data into conditions
    itd5_by_lag_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    itd15_by_lag_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild5_by_lag_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));
    ild15_by_lag_target_onset(isubject,:,:) = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'bash','dash','gash'})')),3));



    %% Broken up by boisition and word type
    % First position
    itd5_by_lead_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    itd5_by_lead_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,small_itd_cond).*ismember(ERP_info_lead_masker(:).Word,{'bash'})')),3));
    itd5_by_lead_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    itd5_by_lead_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,small_itd_cond).*ismember(ERP_info_lead_masker(:).Word,{'dash','gash'})')),3));
    
    itd15_by_lead_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    itd15_by_lead_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,large_itd_cond).*ismember(ERP_info_lead_masker(:).Word,{'bash'})')),3));
    itd15_by_lead_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_itd_cond).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    itd15_by_lead_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,large_itd_cond).*ismember(ERP_info_lead_masker(:).Word,{'dash','gash'})')),3));

    ild5_by_lead_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5_by_lead_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,small_ild_cond).*ismember(ERP_info_lead_masker(:).Word,{'bash'})')),3));
    ild5_by_lead_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,small_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lead_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,small_ild_cond).*ismember(ERP_info_lead_masker(:).Word,{'dash','gash'})')),3));
    
    ild15_by_lead_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild15_by_lead_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,large_ild_cond).*ismember(ERP_info_lead_masker(:).Word,{'bash'})')),3));
    ild15_by_lead_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,large_ild_cond).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild15_by_lead_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lead_masker_onset_baselined(:,:,logical(ismember(ERP_info_lead_masker(:).Condition,large_ild_cond).*ismember(ERP_info_lead_masker(:).Word,{'dash','gash'})')),3));

    % Second position
    itd5_by_lag_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    itd5_by_lag_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,small_itd_cond).*ismember(ERP_info_lag_masker(:).Word,{'bash'})')),3));
    itd5_by_lag_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    itd5_by_lag_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,small_itd_cond).*ismember(ERP_info_lag_masker(:).Word,{'dash','gash'})')),3));

    itd15_by_lag_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    itd15_by_lag_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,large_itd_cond).*ismember(ERP_info_lag_masker(:).Word,{'bash'})')),3));
    itd15_by_lag_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_itd_cond).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    itd15_by_lag_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,large_itd_cond).*ismember(ERP_info_lag_masker(:).Word,{'dash','gash'})')),3));

    ild5_by_lag_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5_by_lag_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,small_ild_cond).*ismember(ERP_info_lag_masker(:).Word,{'bash'})')),3));
    ild5_by_lag_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,small_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lag_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,small_ild_cond).*ismember(ERP_info_lag_masker(:).Word,{'dash','gash'})')),3));
   
    ild15_by_lag_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild15_by_lag_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,large_ild_cond).*ismember(ERP_info_lag_masker(:).Word,{'bash'})')),3));
    ild15_by_lag_word_not_bash_in_target(isubject,:,:)  = squeeze(nanmean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,large_ild_cond).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild15_by_lag_word_not_bash_in_masker(isubject,:,:)  = squeeze(nanmean(data_by_lag_masker_onset_baselined(:,:,logical(ismember(ERP_info_lag_masker(:).Condition,large_ild_cond).*ismember(ERP_info_lag_masker(:).Word,{'dash','gash'})')),3));


    %% Button Press
    
    all_subjects_button_press(isubject,:,:) = squeeze(nanmean(data_by_button_press_baselined,3));
end

curr_channel_index = frontocentral_channels; 

%% Comparing in first position
% Bash in target vs. Bash in masker vs. No bash at all
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % itd5
hold on
this_bash_target_data = squeeze(nanmean(itd5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd5_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd5_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine,p4(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','DashGash in Target Stream','DashGash in Masker Stream'})

subplot(1,4,2) % itd15
hold on
this_bash_target_data = squeeze(nanmean(itd15_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd15_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd15_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd15_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(nanmean(ild5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild5_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild5_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ild15
hold on
this_bash_target_data = squeeze(nanmean(ild15_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild15_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild15_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild15_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

sgtitle('Frontocentral ERP Leading Position Comparisons','FontSize',18)

%% Comparing in second position
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % itd5
hold on
this_bash_target_data = squeeze(nanmean(itd5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd5_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd5_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine,p4(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','DashGash in Target Stream','DashGash in Masker Stream'})

subplot(1,4,2) % itd15
hold on
this_bash_target_data = squeeze(nanmean(itd15_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd15_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd15_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd15_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(nanmean(ild5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild5_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild5_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)


subplot(1,4,4) % ild15
hold on
this_bash_target_data = squeeze(nanmean(ild15_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild15_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild15_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild15_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)

sgtitle('Frontocentral ERP Laggig Position Comparisons','FontSize',18)


%% P300 plot: BASH in first pos. target minus BASH in first pos. non-target
curr_channel_index = parietooccipital_channels;
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % itd5
hold on
this_bash_target_data = squeeze(nanmean(itd5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd5_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd5_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine,p4(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','DashGash in Target Stream','DashGash in Masker Stream'})

subplot(1,4,2) % itd15
hold on
this_bash_target_data = squeeze(nanmean(itd15_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd15_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd15_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd15_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(nanmean(ild5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild5_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild5_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ild15
hold on
this_bash_target_data = squeeze(nanmean(ild15_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild15_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild15_by_lead_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild15_by_lead_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


sgtitle('Parietooccipital ERP Leading Position Comparisons','FontSize',18)

%% P300 plot: BASH in second pos. target minus BASH in second pos. non-target
curr_channel_index = parietooccipital_channels;
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % itd5
hold on
this_bash_target_data = squeeze(nanmean(itd5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd5_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd5_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine,p4(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','DashGash in Target Stream','DashGash in Masker Stream'})

subplot(1,4,2) % itd15
hold on
this_bash_target_data = squeeze(nanmean(itd15_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(itd15_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(itd15_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(itd15_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ITDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(nanmean(ild5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild5_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild5_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)


subplot(1,4,4) % ild15
hold on
this_bash_target_data = squeeze(nanmean(ild15_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(nanmean(ild15_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data_in_target = squeeze(nanmean(ild15_by_lag_word_not_bash_in_target(:,curr_channel_index,:),2));
this_nonbash_data_in_masker = squeeze(nanmean(ild15_by_lag_word_not_bash_in_masker(:,curr_channel_index,:),2));

p1 = shadedErrorBar(single_onset_time,nanmean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,nanmean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_target,1),std(this_nonbash_data_in_target,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
p4 = shadedErrorBar(single_onset_time,nanmean(this_nonbash_data_in_masker,1),std(this_nonbash_data_in_masker,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-y'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('15 deg ILDs','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(350,'--k','LineWidth',2)

sgtitle('Parietooccipital ERP Lagging Position Comparisons','FontSize',18)



% %% Plot All Words

% curr_channel_index = frontocentral_channels;


% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % itd5
% hold on
% this_lead_data = squeeze(nanmean(itd5_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd5_by_lag_target_onset(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % itd15
% hold on
% this_lead_data = squeeze(nanmean(itd15_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd15_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(nanmean(ild5_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild5_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ild15
% hold on
% this_lead_data = squeeze(nanmean(ild15_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild15_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% sgtitle('All Words','FontSize',18)
% 
% 
% %% Plot Bash Only
% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % itd5
% hold on
% this_lead_data = squeeze(nanmean(itd5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % itd15
% hold on
% this_lead_data = squeeze(nanmean(itd15_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd15_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(nanmean(ild5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ild15
% hold on
% this_lead_data = squeeze(nanmean(ild15_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild15_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% sgtitle('BASH Only','FontSize',18)
% 
% %% Plot {Dash, Gash} Only
% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % itd5
% hold on
% this_lead_data = squeeze(nanmean(itd5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % itd15
% hold on
% this_lead_data = squeeze(nanmean(itd15_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd15_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ITDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(nanmean(ild5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ild15
% hold on
% this_lead_data = squeeze(nanmean(ild15_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild15_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ILDs','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% sgtitle('{DASH,GASH} Only','FontSize',18)
% 
% 
% %% BASH in lead position vs. dash,gash in lead position
% 
% ymin = -6;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % itd5
% hold on
% this_lead_data = squeeze(nanmean(itd5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
% p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ITDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Bash','{Dash,Gash}'})
% 
% subplot(1,4,2) % itd15
% hold on
% this_lead_data = squeeze(nanmean(itd15_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd15_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ITDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(nanmean(ild5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,4) % ild15
% hold on
% this_lead_data = squeeze(nanmean(ild15_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild15_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ILDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% sgtitle('BASH in lead position vs. {DASH,GASH} in lead position','FontSize',18)
% 
% 
% %% Bash in lag position vs. dash, gash in lag position
% 
% ymin = -6;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % itd5
% hold on
% this_lead_data = squeeze(nanmean(itd5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
% p2 = shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ITDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Bash','{Dash,Gash}'})
% 
% subplot(1,4,2) % itd15
% hold on
% this_lead_data = squeeze(nanmean(itd15_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(itd15_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ITDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(nanmean(ild5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,4) % ild15
% hold on
% this_lead_data = squeeze(nanmean(ild15_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(nanmean(ild15_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,nanmean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,nanmean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('15 deg ILDs','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% sgtitle('BASH in lag position vs. {DASH,GASH} in lag position','FontSize',18)
% 
% %% Plot all button press
button_press_delay = 0;
single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(all_subjects_button_press,3));
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
title('Button PressFC','FontSize',18)


button_press_delay = 0;
single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(all_subjects_button_press,3));
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
title('Button PressPO','FontSize',18)

