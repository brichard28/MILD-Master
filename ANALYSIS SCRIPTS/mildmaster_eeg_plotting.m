%% mildmaster_eeg_plotting.m

% Author: Benjamin Richardson
% 10/14/2024

curr_subject_ID = char('eeg_pilot_1','eeg_pilot_4','eeg_pilot_5','eeg_pilot_6','eeg_pilot_7','eeg_pilot_8');
%char('eeg_pilot_1','eeg_pilot_4','eeg_pilot_5','eeg_pilot_6','eeg_pilot_7','eeg_pilot_8'); % char();
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
    load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'))

    % Remove noisy ERPs
    
    % Plotting parameters
    erp_window_start_time = -50;
    erp_window_end_time = 950;
    button_press_delay =0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_target_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    cz_index = 32;
    

    %% ALL WORDS
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

    %% BASH ONLY
    % sort lead data into conditions
    itd50_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    itd500_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5mag_by_lead_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));

    % sort lag data into conditions
    itd50_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    itd500_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5mag_by_lag_target_onset_bash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));



    %% {DASH,GASH} ONLY
    % sort lead data into conditions
    itd50_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    itd500_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lead_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));

    % sort lag data into conditions
    itd50_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    itd500_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lag_target_onset_dashgash(isubject,:,:) = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));


    %% Broken up by boisition and word type
    % First position
    itd50_by_lead_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    itd50_by_lead_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).OtherWord,{'bash'})')),3));
    itd50_by_lead_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    itd500_by_lead_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    itd500_by_lead_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).OtherWord,{'bash'})')),3));
    itd500_by_lead_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lead_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5_by_lead_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).OtherWord,{'bash'})')),3));
    ild5_by_lead_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lead_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).Word,{'bash'})')),3));
    ild5mag_by_lead_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).OtherWord,{'bash'})')),3));
    ild5mag_by_lead_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).Word,{'dash','gash'})')),3));

    % Second position
    itd50_by_lag_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    itd50_by_lag_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[2,5]).*ismember(ERP_info_lead_target(:).OtherWord,{'bash'})')),3));
    itd50_by_lag_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[2,5]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    itd500_by_lag_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    itd500_by_lag_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[1,6]).*ismember(ERP_info_lead_target(:).OtherWord,{'bash'})')),3));
    itd500_by_lag_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[1,6]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5_by_lag_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5_by_lag_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[3,8]).*ismember(ERP_info_lead_target(:).OtherWord,{'bash'})')),3));
    ild5_by_lag_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[3,8]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));
    ild5mag_by_lag_bash_in_target(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).Word,{'bash'})')),3));
    ild5mag_by_lag_bash_in_masker(isubject,:,:)  = squeeze(mean(data_by_lead_target_onset_baselined(:,:,logical(ismember(ERP_info_lead_target(:).Condition,[4,7]).*ismember(ERP_info_lead_target(:).OtherWord,{'bash'})')),3));
    ild5mag_by_lag_word_not_bash(isubject,:,:)  = squeeze(mean(data_by_lag_target_onset_baselined(:,:,logical(ismember(ERP_info_lag_target(:).Condition,[4,7]).*ismember(ERP_info_lag_target(:).Word,{'dash','gash'})')),3));


    %% Button Press
    all_subjects_button_press(isubject,:,:) = squeeze(mean(data_by_button_press_baselined,3));
end

curr_channel_index = frontocentral_channels; 

%% Comparing in first position
% Bash in target vs. Bash in masker vs. No bash at all
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ITD50
hold on
this_bash_target_data = squeeze(mean(itd50_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd50_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(itd50_by_lead_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('50 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','No Bash'})

subplot(1,4,2) % ITD500
hold on
this_bash_target_data = squeeze(mean(itd500_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd500_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(itd500_by_lead_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('500 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(mean(ild5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(ild5_by_lead_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_bash_target_data = squeeze(mean(ild5mag_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5mag_by_lead_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(ild5mag_by_lead_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

sgtitle('COMPARING WITHIN FIRST (LEADING) POSITION ONLY','FontSize',18)

%% Comparing in second position
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ITD50
hold on
this_bash_target_data = squeeze(mean(itd50_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd50_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(itd50_by_lag_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('50 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)
legend([p1(1).mainLine,p2(1).mainLine,p3(1).mainLine],{'Bash in Target Stream','Bash in Masker Stream','No Bash'})

subplot(1,4,2) % ITD500
hold on
this_bash_target_data = squeeze(mean(itd500_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd500_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(itd500_by_lag_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('500 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(mean(ild5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(ild5_by_lag_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_bash_target_data = squeeze(mean(ild5mag_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5mag_by_lag_bash_in_masker(:,curr_channel_index,:),2));
this_nonbash_data = squeeze(mean(ild5mag_by_lag_word_not_bash(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data,1),std(this_bash_target_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
p2 = shadedErrorBar(single_onset_time,mean(this_bash_masker_data,1),std(this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
p3 = shadedErrorBar(single_onset_time,mean(this_nonbash_data,1),std(this_nonbash_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

sgtitle('COMPARING WITHIN SECOND (LAGGING) POSITION ONLY','FontSize',18)


%% P300 plot: BASH in first pos. target minus BASH in first pos. non-target
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ITD50
hold on
this_bash_target_data = squeeze(mean(itd50_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd50_by_lead_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('50 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,2) % ITD500
hold on
this_bash_target_data = squeeze(mean(itd500_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd500_by_lead_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('500 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(mean(ild5_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5_by_lead_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_bash_target_data = squeeze(mean(ild5mag_by_lead_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5mag_by_lead_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

sgtitle('BASH in first pos. target minus BASH in first pos. non-target','FontSize',19)

%% P300 plot: BASH in second pos. target minus BASH in second pos. non-target
ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
subplot(1,4,1) % ITD50
hold on
this_bash_target_data = squeeze(mean(itd50_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd50_by_lag_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('50 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,2) % ITD500
hold on
this_bash_target_data = squeeze(mean(itd500_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(itd500_by_lag_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('500 us ITD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

subplot(1,4,3) % ILD5
hold on
hold on
this_bash_target_data = squeeze(mean(ild5_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5_by_lag_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)


subplot(1,4,4) % ILD5Mag
hold on
this_bash_target_data = squeeze(mean(ild5mag_by_lag_bash_in_target(:,curr_channel_index,:),2));
this_bash_masker_data = squeeze(mean(ild5mag_by_lag_bash_in_masker(:,curr_channel_index,:),2));
p1 = shadedErrorBar(single_onset_time,mean(this_bash_target_data - this_bash_masker_data,1),std(this_bash_target_data - this_bash_masker_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});

ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time])
ylabel('Voltage (uV)','FontSize',18)
title('5 deg ILD + MAG','FontSize',18)
xline(0,'--k','LineWidth',2)
xline(250,'--k','LineWidth',2)

sgtitle('BASH in second pos. target minus BASH in second pos. non-target','FontSize',19)
% %% Plot All Words
% ymin = -4;
% ymax = 5;
% num_subjects = size(curr_subject_ID,1);
% figure;
% subplot(1,4,1) % ITD50
% hold on
% this_lead_data = squeeze(mean(itd50_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd50_by_lag_target_onset(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('50 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % ITD500
% hold on
% this_lead_data = squeeze(mean(itd500_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd500_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('500 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(mean(ild5_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ILD5Mag
% hold on
% this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD + MAG','FontSize',18)
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
% subplot(1,4,1) % ITD50
% hold on
% this_lead_data = squeeze(mean(itd50_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd50_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('50 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % ITD500
% hold on
% this_lead_data = squeeze(mean(itd500_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd500_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('500 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(mean(ild5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ILD5Mag
% hold on
% this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD + MAG','FontSize',18)
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
% subplot(1,4,1) % ITD50
% hold on
% this_lead_data = squeeze(mean(itd50_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd50_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
% p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('50 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Attend Lead','Attend Lag'})
% 
% subplot(1,4,2) % ITD500
% hold on
% this_lead_data = squeeze(mean(itd500_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd500_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('500 us ITD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(mean(ild5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD','FontSize',18)
% xline(0,'--r','LineWidth',2)
% xline(250,'--b','LineWidth',2)
% 
% 
% subplot(1,4,4) % ILD5Mag
% hold on
% this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-b'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD + MAG','FontSize',18)
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
% subplot(1,4,1) % ITD50
% hold on
% this_lead_data = squeeze(mean(itd50_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd50_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
% p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('50 us ITD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Bash','{Dash,Gash}'})
% 
% subplot(1,4,2) % ITD500
% hold on
% this_lead_data = squeeze(mean(itd500_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd500_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('500 us ITD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(mean(ild5_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,4) % ILD5Mag
% hold on
% this_lead_data = squeeze(mean(ild5mag_by_lead_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5mag_by_lead_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD + MAG','FontSize',18)
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
% subplot(1,4,1) % ITD50
% hold on
% this_lead_data = squeeze(mean(itd50_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd50_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% p1 = shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'});
% p2 = shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'});
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('50 us ITD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% legend([p1(1).mainLine,p2(1).mainLine],{'Bash','{Dash,Gash}'})
% 
% subplot(1,4,2) % ITD500
% hold on
% this_lead_data = squeeze(mean(itd500_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(itd500_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('500 us ITD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,3) % ILD5
% hold on
% this_lead_data = squeeze(mean(ild5_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% 
% subplot(1,4,4) % ILD5Mag
% hold on
% this_lead_data = squeeze(mean(ild5mag_by_lag_target_onset_bash(:,curr_channel_index,:),2));
% this_lag_data = squeeze(mean(ild5mag_by_lag_target_onset_dashgash(:,curr_channel_index,:),2));
% shadedErrorBar(single_onset_time,mean(this_lead_data,1),std(this_lead_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-g'})
% shadedErrorBar(single_onset_time,mean(this_lag_data,1),std(this_lag_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-m'})
% ylim([ymin,ymax])
% xlim([erp_window_start_time,erp_window_end_time])
% ylabel('Voltage (uV)','FontSize',18)
% title('5 deg ILD + MAG','FontSize',18)
% xline(0,'--k','LineWidth',2)
% xline(250,'--k','LineWidth',2)
% 
% sgtitle('BASH in lag position vs. {DASH,GASH} in lag position','FontSize',18)
% 
% %% Plot all button press
button_press_delay = 0;
single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(all_subjects_button_press,3));

ymin = -4;
ymax = 5;
num_subjects = size(curr_subject_ID,1);
figure;
hold on
this_data = squeeze(mean(all_subjects_button_press(:,curr_channel_index,:),2));
plot(single_onset_time_buttonpress,squeeze(mean(all_subjects_button_press(:,curr_channel_index,:),2)),'-k')
p1 = shadedErrorBar(single_onset_time_buttonpress,mean(this_data,1),std(this_data,[],1)./(sqrt(num_subjects - 1)),'lineProps',{'-r'});
ylim([ymin,ymax])
xlim([erp_window_start_time,erp_window_end_time - 500])
ylabel('Voltage (uV)','FontSize',18)
title('Button Press','FontSize',18)


