addpath("C:\Users\maana\Documents\GitHub\Bash_Dash_Gash\Stimcode\Audio\errorbar_files");
addpath('C:\Users\maana\Documents\MATLAB\eeglab2023.0');
addpath(genpath('C:\Users\maana\Documents\MATLAB\eeglab2023.0\functions\'))

curr_subject_ID =  char('fullpilot1','fullpilot2','fullpilot3');

% Set directories
whos_using = 'Bon';
if all(whos_using == 'Ben')
    addpath('/home/ben/Documents/MATLAB/eeglab2023.1');
    dir = '/home/ben/Documents/GitHub/fNIRSandGerbils/';
    dir_mildmaster = '/home/ben/Documents/GitHub/fNIRSandGerbils/data/fNIRSandGerbils.xlsx';
elseif all(whos_using == 'Bon')
    addpath('C:\Users\benri\Documents\eeglab2023.1');
    dir = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
    dir_mildmaster = 'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx';
    prepro_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';
end
for isubject= 1:length(curr_subject_ID)

    %% Load in and prepare the data
    subID = curr_subject_ID(isubject,:); % set subject ID
    
    all_eeg_epoch = load([prepro_folder subID 'all_epoch.mat']);
    all_eeg_epoch_data=all_eeg_epoch.EEG.data(:,:,:); % channels x time x trials raw eeg data
    all_eeg_epoch_time=all_eeg_epoch.EEG.times;

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


    p_o_channels=[[20:31],[57:64]];
    f_c_channels=[[8:19],[43:56]];
    central=[48];

    attend_right_indices = ismember(conditions,[1,2,3,7]);
    attend_left_indices = ismember(conditions,[4,5,6,8]);

    
    fs= all_eeg_epoch.EEG.srate;

    %% Baseline raw eeg data to 100 ms before cue onset
    baseline_duration=0.1;
    baseline_samples=round(baseline_duration*fs);
    [~,index_time0] = min(abs(all_eeg_epoch.EEG.times - 0));

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

    lowpass_cutoff = ((lower_alphafreq)/(all_eeg_epoch.EEG.srate/2)); % normalize frequency
    highpass_cutoff =((higher_alphafreq)/(all_eeg_epoch.EEG.srate/2));
    b = fir1(256,[lowpass_cutoff highpass_cutoff]);

    all_alpha_power = nan(size(all_eeg_epoch_data_baselined));


    for ichannel= 1:numchannels
        for itrial=1:numtrials
            this_EEG= double(squeeze(all_eeg_epoch_data_baselined(ichannel,:,itrial)));
            % Apply Filter
            this_current_alpha = filter(b,1,this_EEG); % OR FILTFILT?!?!

            % Apply Hilbert 
            this_current_alpha_power = abs(hilbert(this_current_alpha));

            % Save to array
            all_alpha_power(ichannel,:,itrial) = this_current_alpha_power;
        end
    end

    %% Plot average alpha power topography for this subject
    figure;
    % normalize for plotting
    %all_alpha_power = (all_alpha_power - min(all_alpha_power,[],'all'))./(max(all_alpha_power,[],'all')- min(all_alpha_power,[],'all'));
    topo_end_time = 7;
    topoplot_indices = round(0:0.5*fs:(topo_end_time - epoch_start_time)*fs);
    topoplot_indices(1) = 1;
    topoplot_times = -1000:500:topo_end_time*1000;
    cmin = 0;
    cmax = 5;

    % All Trials
    iplot = 1;
    itime = 1;
    for itopo = topoplot_indices
        subplot(3,length(topoplot_indices),iplot);
        this_data = squeeze(mean(all_alpha_power(:,itopo,:),[2,3]));
        topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end
    colorbar

    % Attend Left
    itime = 1;
    for itopo = topoplot_indices
        subplot(3,length(topoplot_indices),iplot);
        this_data = squeeze(mean(all_alpha_power(:,itopo,attend_left_indices),[2,3]));
        topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end
    colorbar

    % Attend Right
    itime = 1;
    for itopo = topoplot_indices
        subplot(3,length(topoplot_indices),iplot);
        this_data = squeeze(mean(all_alpha_power(:,itopo,attend_right_indices),[2,3]));
        topoplot(this_data,all_eeg_epoch.EEG.chanlocs,'maplimits',[cmin, cmax]);
        title([num2str(topoplot_times(itime)),' ms'])
        iplot = iplot + 1;
        itime = itime + 1;
    end
    colorbar


    sgtitle(subID)

    %% TO DO
    % worry about IPAF
    % subtract mean ERP
    % split by condition

    %% Plot spectrogram for this subject
    M = fs/10;
    window=hamming(M,'periodic');
    noverlap=floor(0.75*M);
    Ndft=128;
    all_spectrograms_this_subject = nan(Ndft/2 + 1, 619, 32, size(all_eeg_epoch_data_baselined,3));

    for ichannel = 1:32
        for itrial = 1:size(all_eeg_epoch_data_baselined,3) % find the spectrogram for each trial

            this_data = mean(all_eeg_epoch_data_baselined(ichannel,:,itrial),[1,3]);
            [this_spectrogram,f,t] = spectrogram(this_data,window,noverlap,Ndft,fs,'yaxis'); % will return frequency by time
            all_spectrograms_this_subject(:,:,ichannel, itrial) = abs(this_spectrogram); % frequency x time x channel X trial
        end
    end
    
        % Frontocentral channels
        frontocentral_channels = [1,2,4,5,6,8,9,23,25,26,27,29,31,32];

    figure;
    hold on
    imagesc(t,f,squeeze(mean(all_spectrograms_this_subject(:,:,frontocentral_channels,:),[3,4])))
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([20,70])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    %yline(8,'LineWidth',2)
    %yline(12,'LineWidth',2)
    title('Frontocentral Channels','FontSize',18)
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)

    % Parietooccipital channels BROKEN UP BY ATTEND
    parietooccipital_channels = 11:20;
    left_parietooccipital_channels = [11,12,14,15];
    right_parietooccipital_channels = [17,18,19,20];
%         condition_names = {'side=r_itd=500_az=0_mag=0_lpf=0',...
%         'side=r_itd=50_az=0_mag=0_lpf=0',...
%         'side=r_itd=0_az=5_mag=0_lpf=0',...
%         'side=l_itd=0_az=5_mag=1_lpf=0',...
%         'side=l_itd=50_az=0_mag=0_lpf=0',...
%         'side=l_itd=500_az=0_mag=0_lpf=0',...
%         'side=r_itd=0_az= 5_mag=1_lpf=0',...
%         'side=l_itd=0_az=5_mag=0_lpf=0'};
    figure;
    subplot(2,2,1) % left hemisphere, attend left
    imagesc(t,f,squeeze(mean(all_spectrograms_this_subject(:,:,left_parietooccipital_channels,ismember(conditions,[4,5,6,8])),[3,4])))
    set(gca,'YDir','normal')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([20,50])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Left Hem. Attend Left')
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)

    subplot(2,2,2) % right hemisphere, attend left
    imagesc(t,f',squeeze(mean(all_spectrograms_this_subject(:,:,right_parietooccipital_channels,ismember(conditions,[4,5,6,8])),[3,4])))
    set(gca,'YDir','normal')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([20,50])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Right Hem. Attend Left')
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)

    subplot(2,2,3) % left hemisphere, attend right
    imagesc(t,f',squeeze(mean(all_spectrograms_this_subject(:,:,left_parietooccipital_channels,ismember(conditions,[1,2,3,7])),[3,4])))
    set(gca,'YDir','normal')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([20,50])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Left Hem. Attend Right')
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)

    subplot(2,2,4) % right hemisphere, attend right
    imagesc(t,f',squeeze(mean(all_spectrograms_this_subject(:,:,right_parietooccipital_channels,ismember(conditions,[1,2,3,7])),[3,4])))
    set(gca,'YDir','normal')
    ylim([0,30])
    xlim([0,6])
    xticklabels(0:1:16);
    curr_tick_labels = xticklabels;
    xticklabels(string(str2double(string(curr_tick_labels)) - 1));
    caxis([20,50])
    xline(1,'LineWidth',2)
    xline(2.8,'LineWidth',2)
    yline(8,'LineWidth',2)
    yline(12,'LineWidth',2)
    title('Right Hem. Attend Right')
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)


    sgtitle('Parietooccipital Channels','FontSize',18)


end