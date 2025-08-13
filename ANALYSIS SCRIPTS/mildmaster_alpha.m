addpath(genpath('C:\Users\benri\Documents\eeglab2023.0'));

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

% Set directories

addpath('C:\Users\benri\Documents\eeglab2023.0');
dir = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
dir_mildmaster = 'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx';
prepro_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';

left_hem_attend_left = nan(size(curr_subject_ID,1),4352);
left_hem_attend_right = nan(size(curr_subject_ID,1),4352);
right_hem_attend_left = nan(size(curr_subject_ID,1),4352);
right_hem_attend_right = nan(size(curr_subject_ID,1),4352);

for isubject= 1:size(curr_subject_ID,1)
    

    %% Load in and prepare the data
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID)

    all_eeg_epoch = load([prepro_folder char(strtrim(string(subID))) 'all_epoch_yes_button_press.mat']);
    all_eeg_epoch_data=all_eeg_epoch.EEG_yes_button_press.data(:,:,:); % channels x time x trials raw eeg data
    all_eeg_epoch_time=all_eeg_epoch.EEG_yes_button_press.times;

    numchannels = size(all_eeg_epoch_data,1);
    numtrials= size(all_eeg_epoch_data,3);
    numtimepoints = size(all_eeg_epoch_data,2);
    epoch_start_time=-1;
    epoch_end_time=16;

    % Load word onset times data
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(subID) + "__Word_Times.csv");

    % Read in click times, find the rows in the table for this subject
    BehaviorTable = readtable(dir_mildmaster,'FileType','spreadsheet','Format','auto');
    rows_this_subject = find(BehaviorTable.S == string(curr_subject_ID(isubject,:))); % find the rows in the spreadsheet which belong to this subject
    conditions = BehaviorTable.Condition(rows_this_subject); % conditions by trial for this subject
    condition_names = {'side=r_itd=500_az=0_mag=0_lpf=0',...
        'side=r_itd=50_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az=5_mag=0_lpf=0',...
        'side=l_itd=0_az=5_mag=1_lpf=0',...
        'side=l_itd=50_az=0_mag=0_lpf=0',...
        'side=l_itd=500_az=0_mag=0_lpf=0',...
        'side=r_itd=0_az= 5_mag=1_lpf=0',...
        'side=l_itd=0_az=5_mag=0_lpf=0'};
    this_subject_table = BehaviorTable(rows_this_subject,:);


    p_o_channels=[11:20];
    f_c_channels= [1,2,4,5,6,8,9,23,25,26,27,29,31,32];
    central= [32];
    pz_channel = [13];

    if subID == "mild_master_20" || subID == "mild_master_27" || subID == "mild_master_37"
        conditions(end) = [];
    end

    attend_right_indices = ismember(conditions,[1,2,3,7]);
    attend_left_indices = ismember(conditions,[4,5,6,8]);


    fs= all_eeg_epoch.EEG_yes_button_press.srate;

    %% Baseline raw eeg data to 100 ms before cue onset
    baseline_duration=1;
    baseline_samples=round(baseline_duration*fs);
    [~,index_time0] = min(abs(all_eeg_epoch.EEG_yes_button_press.times - 0));

    all_eeg_epoch_data_baselined = nan(size(all_eeg_epoch_data));

    for ichannel= 1:numchannels
        for itrial=1:numtrials
            all_eeg_epoch_data_baselined(ichannel,:,itrial) = all_eeg_epoch_data(ichannel,:,itrial) - mean(all_eeg_epoch_data(ichannel,1:index_time0,:),'all');
        end
    end

    %% Calculate IPAF
    % UNDER CONSTRUCTION FOR NOW

    %% Extract alpha power between 8-12 Hz
    lower_alphafreq = 8;
    higher_alphafreq = 12;

    lowpass_cutoff = ((lower_alphafreq)/(all_eeg_epoch.EEG_yes_button_press.srate/2)); % normalize frequency
    highpass_cutoff =((higher_alphafreq)/(all_eeg_epoch.EEG_yes_button_press.srate/2));
    b = fir1(128,[lowpass_cutoff highpass_cutoff]);

    all_alpha_power = nan(size(all_eeg_epoch_data_baselined));

    % SUBTRACT MEAN ERP
    all_eeg_epoch_data_baselined = all_eeg_epoch_data_baselined - mean(all_eeg_epoch_data_baselined,3);


    for ichannel= 1:numchannels
        for itrial=1:numtrials
            this_EEG= double(squeeze(all_eeg_epoch_data_baselined(ichannel,:,itrial)));
            % Apply Filter
            this_current_alpha = filtfilt(b,1,this_EEG); % OR FILTFILT?!?!

            % Apply Hilbert
            this_current_alpha_power = abs(hilbert(this_current_alpha));
            % Save to array
            all_alpha_power(ichannel,:,itrial) = this_current_alpha_power;
        end
        % Baseline for this channel
        all_alpha_power(ichannel,:,:) = all_alpha_power(ichannel,:,:) - squeeze(mean(all_alpha_power(ichannel,1:fs,:),'all'));
    end

    left_hem_attend_left(isubject,:) = squeeze(mean(all_alpha_power(p_o_channels([1,2,4,5]),:,attend_left_indices),[1,3])); % left hemisphere attend left
    left_hem_attend_right(isubject,:) =squeeze(mean(all_alpha_power(p_o_channels([1,2,4,5]),:,attend_right_indices),[1,3])); % left hemisphere attend right
    right_hem_attend_left(isubject,:) = squeeze(mean(all_alpha_power(p_o_channels(7:10),:,attend_left_indices),[1,3])); % right hemisphere attend left
    right_hem_attend_right(isubject,:) = squeeze(mean(all_alpha_power(p_o_channels(7:10),:,attend_right_indices),[1,3])); % right hemisphere attend right
    



%     %% Plot average alpha power topography for this subject
%     figure;
%     % normalize for plotting
%     %all_alpha_power = (all_alpha_power - min(all_alpha_power,[],'all'))./(max(all_alpha_power,[],'all')- min(all_alpha_power,[],'all'));
%     topo_end_time = 7;
%     topoplot_indices = round(0:0.5*fs:(topo_end_time - epoch_start_time)*fs);
%     topoplot_indices(1) = 1;
%     topoplot_times = -1000:500:topo_end_time*1000;
%     cmin = 0;
%     cmax = 5;
% 
%     % All Trials
%     iplot = 1;
%     itime = 1;
%     for itopo = topoplot_indices
%         subplot(3,length(topoplot_indices),iplot);
%         this_data = squeeze(mean(all_alpha_power(:,itopo,:),[2,3]));
%         topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
%         title([num2str(topoplot_times(itime)),' ms'])
%         iplot = iplot + 1;
%         itime = itime + 1;
%     end
%     colorbar
% 
%     % Attend Left
%     itime = 1;
%     for itopo = topoplot_indices
%         subplot(3,length(topoplot_indices),iplot);
%         this_data = squeeze(mean(all_alpha_power(:,itopo,attend_left_indices),[2,3]));
%         topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
%         title([num2str(topoplot_times(itime)),' ms'])
%         iplot = iplot + 1;
%         itime = itime + 1;
%     end
%     colorbar
% 
%     % Attend Right
%     itime = 1;
%     for itopo = topoplot_indices
%         subplot(3,length(topoplot_indices),iplot);
%         this_data = squeeze(mean(all_alpha_power(:,itopo,attend_right_indices),[2,3]));
%         topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
%         title([num2str(topoplot_times(itime)),' ms'])
%         iplot = iplot + 1;
%         itime = itime + 1;
%     end
%     colorbar
% 
% 
%     sgtitle(subID)

    %% TO DO
    % worry about IPAF
    % subtract mean ERP
    % split by condition

    %% Plot CWT spectrogram for this subject at Pz

%     Ndft=128;
%     time_window = linspace(epoch_start_time, epoch_end_time, size(all_eeg_epoch_data, 2));
% 
%      all_spectrograms_this_subject = nan(Ndft/2 + 1, 123, 32, size(all_eeg_epoch_data,3));
%      all_cwt_this_subject= cell(length(p_o_channels),numtrials);
%      freq_range= [1 50];
%     for ichannel=1:length(p_o_channels)
%         channel_idx=p_o_channels(ichannel);
%         for itrial=1:numtrials
%             spect_sub= all_eeg_epoch_data_baselined(channel_idx,:,itrial);
%             %[this_spectrogram,f,t]=spectrogram(spect_sub,window,noverlap,Ndft,fs,'yaxis');
%             %all_spectrograms_this_subject(:,:,ichannel, itrial) = abs(this_spectrogram);
%             [cwt_coeffs, frequencies] = cwt(spect_sub, fs, 'FrequencyLimits', freq_range);
%             all_cwt_this_subject{ichannel,itrial}=abs(cwt_coeffs);
%         end
%     end
% 
% 
% 
%     %all_po_cwt_this_subject = all_cwt_this_subject(p_o_channels,:);
%     po_cwt_to_plot = [];
%     for ichannel=1:length(p_o_channels)
%         for itrial = 1:numtrials
%             po_cwt_to_plot(:,:,ichannel,itrial) = all_cwt_this_subject{ichannel,itrial};
%         end
%     end
% 
% 
%     figure;
%     imagesc(time_window,frequencies,mean(po_cwt_to_plot(:,:,3,:),[3,4]))
%     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
%     axis tight
%     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
%     ylabel('Frequency(Hz)','FontSize',18)
%     xlabel('Time (s)','FontSize',18)
%     title(['CWT Spectrogram at Pz ',subID],'FontSize',18)
%     xline(0,'LineWidth',2)
%     xline(2,'LineWidth',2)
%     xline(11.6,'LineWidth',2)
%     xlim([-1,12])
%     yline(8,'LineWidth',2)
%     yline(12,'LineWidth',2)
% 
% 
%     %% Plot CWT spectrogram for this subject at Left Hemisphere Attend Left vs. Right
% 
%      figure;
%      subplot(1,2,1)
%     imagesc(time_window,frequencies,mean(po_cwt_to_plot(:,:,[1,2,4,5],attend_left_indices),[3,4]))
%     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
%     axis tight
%     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
%     ylabel('Frequency(Hz)','FontSize',18)
%     xlabel('Time (s)','FontSize',18)
%     title(['Left Hemisphere ',subID, ' Attend Left'],'FontSize',18)
%     xline(0,'LineWidth',2)
%     xline(2,'LineWidth',2)
%     xline(11.6,'LineWidth',2)
%     xlim([-1,12])
%     yline(8,'LineWidth',2)
%     yline(12,'LineWidth',2)
% 
% 
%     subplot(1,2,2)
%     imagesc(time_window,frequencies,mean(po_cwt_to_plot(:,:,[1,2,4,5],attend_right_indices),[3,4]))
%     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
%     axis tight
%     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
%     ylabel('Frequency(Hz)','FontSize',18)
%     xlabel('Time (s)','FontSize',18)
%     title(['Left Hemisphere ',subID, ' Attend Right'],'FontSize',18)
%     xline(0,'LineWidth',2)
%     xline(2,'LineWidth',2)
%     xline(11.6,'LineWidth',2)
%     xlim([-1,12])
%     yline(8,'LineWidth',2)
%     yline(12,'LineWidth',2)
% 
%     %% Plot CWT spectrogram for this subject Right Hemisphere Attend Left vs. Right
% 
%     figure;
%     subplot(1,2,1)
%     imagesc(time_window,frequencies,mean(po_cwt_to_plot(:,:,7:10,attend_left_indices),[3,4]))
%     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
%     axis tight
%     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
%     ylabel('Frequency(Hz)','FontSize',18)
%     xlabel('Time (s)','FontSize',18)
%     title(['Right Hemisphere ',subID, ' Attend Left'],'FontSize',18)
%     xline(0,'LineWidth',2)
%     xline(2,'LineWidth',2)
%     xline(11.6,'LineWidth',2)
%     xlim([-1,12])
%     yline(8,'LineWidth',2)
%     yline(12,'LineWidth',2)
% 
% 
%     subplot(1,2,2)
%     imagesc(time_window,frequencies,mean(po_cwt_to_plot(:,:,7:10,attend_right_indices),[3,4]))
%     set(gca,'YScale','log','YMinorTick','off','Ydir','normal')
%     axis tight
%     set(gca,'YTick',freq_range(1):1:freq_range(end),'YTickLabel',freq_range(1):1:freq_range(end))
%     ylabel('Frequency(Hz)','FontSize',18)
%     xlabel('Time (s)','FontSize',18)
%     title(['Right Hemisphere ',subID, ' Attend Right'],'FontSize',18)
%     xline(0,'LineWidth',2)
%     xline(2,'LineWidth',2)
%     xline(11.6,'LineWidth',2)
%     xlim([-1,12])
%     yline(8,'LineWidth',2)
%     yline(12,'LineWidth',2)
end


%% Plot alpha time traces across subjects
figure
% Left Hemisphere
subplot(1,2,1)
hold on
% Attend left
shadedErrorBar(all_eeg_epoch.EEG_yes_button_press.times,nanmean(left_hem_attend_left,1),nanstd(left_hem_attend_left,[],1)./sqrt(size(curr_subject_ID,1)),'lineProps','-r')
% Attend right
shadedErrorBar(all_eeg_epoch.EEG_yes_button_press.times,nanmean(left_hem_attend_right,1),nanstd(left_hem_attend_right,[],1)./sqrt(size(curr_subject_ID,1)),'lineProps','-b')
legend({'Attend Left','Attend Right'})
title('Left PO Channels','FontSize',20)
%ylim([-2,1])

% Right Hemisphere
subplot(1,2,2)
hold on
% Attend left
shadedErrorBar(all_eeg_epoch.EEG_yes_button_press.times,nanmean(right_hem_attend_left,1),nanstd(right_hem_attend_left,[],1)./sqrt(size(curr_subject_ID,1)),'lineProps','-r')
% Attend right
shadedErrorBar(all_eeg_epoch.EEG_yes_button_press.times,nanmean(right_hem_attend_right,1),nanstd(right_hem_attend_right,[],1)./sqrt(size(curr_subject_ID,1)),'lineProps','-b')
legend({'Attend Left','Attend Right'})
title('Right PO Channels','FontSize',20)
%ylim([-2,1])

