%% mildmaster_button_press_plotting
% Benjamin Richardson

% Set directories
whos_using = 'Bon';
if all(whos_using == 'Ben')
    addpath('/home/ben/Documents/MATLAB/eeglab2023.1');
    dir = '/home/ben/Documents/GitHub/MILD-Master/';
    dir_mildmaster = '/home/ben/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER Behavior Files/button-press-pilot.xlsx';
elseif all(whos_using == 'Bon')
    addpath('C:\Users\benri\Documents\eeglab2023.1');
    dir = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
    dir_mildmaster = 'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\button-press-pilot.xlsx';
    prepro_folder = 'C:\Users\benri\Documents\GitHub\MILD-Master\prepro_epoched_data\';
end

curr_subject_ID = char('mild_master_3',...
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
    'mild_master_23'); % char();

% Set analysis parameters
erp_window_start_time = -1000; % 100 ms before onset of word
erp_window_end_time = 1000; % 750 ms after onset of word
nsubjects = size(curr_subject_ID,1);
word_length = 0.3;
frontocentral_channels = 32; %[1,2,4,5,6,8,9,23,25,26,27,29,31,32];
parietooccipital_channels = 13; %11:20;
fs = 256;
all_mean_button_press = [];

for isubject = 1:size(curr_subject_ID,1)
    subID = curr_subject_ID(isubject,:); % set subject ID
    disp(subID)

    % load data
    load(append(strtrim(string(curr_subject_ID(isubject,:))),'_button_press_data.mat'))

    
    PressTimesTable = readtable(append('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_',strtrim(string(curr_subject_ID(isubject,:))),'__button_press_times.csv'),'ReadVariableNames',false,'FileType','spreadsheet','Range','A2:ET2');

    eeg_time = 0:1/fs:((length(button_press_data) - 1)/fs);

    button_press_data_this_subject = [];
    for ipress = 1:width(PressTimesTable)
        this_time=strtrim(string(PressTimesTable{:,ipress}));
        this_time=erase(this_time,"'");
        this_time=erase(this_time,"[");
        this_time=erase(this_time,"]");
        this_time=str2double(this_time);

        this_time = this_time + 702;

        this_search_time = this_time/1000;
        [~,this_eeg_index_start] = min(abs(eeg_time - (this_search_time + (erp_window_start_time/1000))));
        [~,this_eeg_index_end] = min(abs(eeg_time - (this_search_time + (erp_window_end_time/1000))));

        if this_eeg_index_end - this_eeg_index_start == 1537
            this_eeg_index_end = this_eeg_index_end + 1;
        end

        button_press_data_this_subject = cat(3,button_press_data_this_subject,button_press_data(:,this_eeg_index_start:this_eeg_index_end));
    end

    all_mean_button_press(isubject,:,:) = squeeze(mean(button_press_data_this_subject,3));



%     load(append('Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'.mat'))
% 
%     num_button_presses_to_use = 5:5:100; %10:10:size(data_by_button_press_baselined,3);
%     
%     single_onset_time_buttonpress = linspace(erp_window_start_time,erp_window_end_time,size(data_by_button_press_baselined,2));
% 
%     
%     colorMapLength = length(num_button_presses_to_use);
%     red = [1, 0, 0];
%     pink = [255, 192, 203]/255;
%     blue = [0,0,1];
%     colors_p = [linspace(red(1),blue(1),colorMapLength)', linspace(red(2),blue(2),colorMapLength)', linspace(red(3),blue(3),colorMapLength)'];
%    
%     
    % plot in frontocentral channels at several steps of number of button
    % presses
%     figure;
%     hold on
%     colormap(colors_p);
%     for numpresses = num_button_presses_to_use
%         indices_to_use = randsample(size(data_by_button_press_baselined,3),numpresses);
%         this_data_to_plot = squeeze(mean(data_by_button_press_baselined(frontocentral_channels,:,indices_to_use),[1,3]));
%         plot(single_onset_time_buttonpress,this_data_to_plot,'Color',colors_p(find(numpresses == num_button_presses_to_use),:))
% 
% 
%     end
%     cbh = colorbar();
%     caxis([0,1])
%     cbh.Ticks = linspace(0,1,length(num_button_presses_to_use));
%     cbh.TickLabels = num2cell(num_button_presses_to_use);
%     xlim([erp_window_start_time,erp_window_end_time])
%     xlabel('Time re: button press onset (ms)','FontSize',18)
%     ylabel('Voltage (uV)','FontSize',18)
%     title(append('Button Presses in Frontocentral Channels Subject ',num2str(isubject)),'FontSize',18)
%     
% 
%     % plot in parietooccipital channels at several steps of number of
%     % button presses
%     figure;
%     hold on
%     colormap(colors_p);
%     for numpresses = num_button_presses_to_use
%         this_data_to_plot = squeeze(mean(data_by_button_press_baselined(parietooccipital_channels,:,1:numpresses),[1,3]));
%         plot(single_onset_time_buttonpress,this_data_to_plot,'Color',colors_p(find(numpresses == num_button_presses_to_use),:))
% 
% 
%     end
%     cbh = colorbar();
%     caxis([0,1])
%     cbh.Ticks = linspace(0,1,length(num_button_presses_to_use));
%     cbh.TickLabels = num2cell(num_button_presses_to_use);
%     xlim([erp_window_start_time,erp_window_end_time])
%     xlabel('Time re: button press onset (ms)','FontSize',18)
%     ylabel('Voltage (uV)','FontSize',18)
%     title(append('Button Presses in Parietooccipital Channels Subject ',num2str(isubject)),'FontSize',18)

end

single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(all_mean_button_press,3));

%figure;plot(single_onset_time,squeeze(mean(all_mean_button_press(:,[1,2,4,5,6,8,9,23,25,26,27,29,31,32],:),2))')
figure;plot(single_onset_time,squeeze(mean(all_mean_button_press(:,[2,3,29],:),2))')

xlabel('Time re: button press (ms)')
title('Frontocentral Channels','FontSize',18)

figure;plot(single_onset_time,squeeze(mean(all_mean_button_press(:,11:20,:),2))')
xlabel('Time re: button press (ms)')
title('Parietooccipital Channels','FontSize',18)