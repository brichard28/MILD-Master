%% mildmaster_eeg_calculations.m

% Calculate P1, N1, and P300 statistics for mild master
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
    'mild_master_32','mild_master_33','mild_master_34','mild_master_36','mild_master_37','mild_master_38','mild_master_39','mild_master_40', ...
    'mild_master_41','mild_master_42','mild_master_43','mild_master_44','mild_master_46','mild_master_48'); % char();
num_subjects = size(curr_subject_ID,1);

mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master\';
erp_window_start_time = -50; % 100 ms before onset of word
erp_window_end_time = 1500; % 750 ms after onset of word
fs = 256;

all_subs_fig_p1n1p2 = figure();
all_subs_fig_lead_p3 = figure();
all_subs_fig_lag_p3 = figure();

all_subs_p1 = struct('S','','Condition','','Lead_Stream','','Lead_Word','','Lag_Word','','Lead_Amplitude','','Lag_Amplitude','','Electrode','');
all_subs_n1 = struct('S','','Condition','','Lead_Stream','','Lead_Word','','Lag_Word','','Lead_Amplitude','','Lag_Amplitude','','Electrode','');
all_subs_p2 = struct('S','','Condition','','Lead_Stream','','Lead_Word','','Lag_Word','','Lead_Amplitude','','Lag_Amplitude','','Electrode','');
all_subs_p3 = struct('S','','Condition','','Lead_Stream','','Lead_Word','','Lag_Word','','Lead_Amplitude','','Lag_Amplitude','','Electrode','');

electrode_names = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5','P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8','CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2','Fz','Cz'};

structrow = 1;
for isubject = 1:size(curr_subject_ID,1)
    subID = string(curr_subject_ID(isubject,:));
    disp(subID)
    % Load Data
    load(append('C:\Users\benri\Downloads\Results_Subject_',strtrim(string(curr_subject_ID(isubject,:))),'_yes_button_press.mat'))

    %% Plot all channels, remove noisy ones in time domain

    % Plotting parameters
    button_press_delay = 0 ;
    single_onset_time = linspace(erp_window_start_time,erp_window_end_time,size(data_by_pair_onset_baselined,2));
    single_onset_time_buttonpress = linspace(erp_window_start_time + button_press_delay,erp_window_end_time,size(data_by_button_press_baselined,2));
    frontocentral_channels = [31, 5, 26, 8, 32, 23, 9, 22]; % Fz, FC1, FC2, C3, Cz, C4, CP1, and CP2
    parietooccipital_channels = [12, 13, 14, 15, 16, 17, 18, 19] ;%  P3, Pz, PO3, O1, Oz, O2, PO4, and P4
    cz_index = 32;
    fz_index = 31;
    pz_index = 13;


    %% Find this subjects' P1, N1, and P300 latency using the mean over ALL lead lag tokens
    % Do this separately for lead, lag P1 and N1
    % We will use Cz and Fz average to calculate P1 and N1, and Pz to
    % calculate P300

    % For P300, only basing on trials that have bash!

    [~,lead_p1_start_index] = min(abs(single_onset_time - (80)));
    [~,lead_p1_end_index] = min(abs(single_onset_time - (150)));
    [~,lead_n1_start_index] = min(abs(single_onset_time - (120)));
    [~,lead_n1_end_index] = min(abs(single_onset_time - (160)));
    [~,lead_p2_start_index] = min(abs(single_onset_time - (160)));
    [~,lead_p2_end_index] = min(abs(single_onset_time - (230)));
    [~,lead_p3_start_index] = min(abs(single_onset_time - (400)));
    [~,lead_p3_end_index] = min(abs(single_onset_time - (600)));


    [~,lag_p1_start_index] = min(abs(single_onset_time - (340)));
    [~,lag_p1_end_index] = min(abs(single_onset_time - (400)));
    [~,lag_n1_start_index] = min(abs(single_onset_time - (380)));
    [~,lag_n1_end_index] = min(abs(single_onset_time - (460)));
    [~,lag_p2_start_index] = min(abs(single_onset_time - (430)));
    [~,lag_p2_end_index] = min(abs(single_onset_time - (520)));
    [~,lag_p3_start_index] = min(abs(single_onset_time - (650)));
    [~,lag_p3_end_index] = min(abs(single_onset_time - (900)));

    % take the average
    this_sub_cz_fz_average = squeeze(mean(data_by_pair_onset_baselined([frontocentral_channels],:,:),[1,3]));
    this_sub_lead_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lead_Word,{'bash'})),[1,3]));
    this_sub_lag_bash_pz_average = squeeze(mean(data_by_pair_onset_baselined([parietooccipital_channels],:,ismember(ERP_info(:).Lag_Word,{'bash'})),[1,3]));

    % lead p1
    [~,this_sub_lead_p1_index] = max(this_sub_cz_fz_average(lead_p1_start_index:lead_p1_end_index));
    this_sub_lead_p1_index = this_sub_lead_p1_index + lead_p1_start_index - 1;
    this_sub_lead_p1_time = single_onset_time(this_sub_lead_p1_index);

    % lead n1
    [~,this_sub_lead_n1_index] = min(this_sub_cz_fz_average(lead_n1_start_index:lead_n1_end_index));
    this_sub_lead_n1_index = this_sub_lead_n1_index + lead_n1_start_index - 1;
    this_sub_lead_n1_time = single_onset_time(this_sub_lead_n1_index);

    % lead p2
    [~,this_sub_lead_p2_index] = max(this_sub_cz_fz_average(lead_p2_start_index:lead_p2_end_index));
    this_sub_lead_p2_index = this_sub_lead_p2_index + lead_p2_start_index - 1;
    this_sub_lead_p2_time = single_onset_time(this_sub_lead_p2_index);


    % lead p3
    [~,this_sub_lead_p3_index] = max(this_sub_lead_bash_pz_average(lead_p3_start_index:lead_p3_end_index));
    this_sub_lead_p3_index = this_sub_lead_p3_index + lead_p3_start_index - 1;
    this_sub_lead_p3_time = single_onset_time(this_sub_lead_p3_index);


   % lag p1
    [~,this_sub_lag_p1_index] = max(this_sub_cz_fz_average(lag_p1_start_index:lag_p1_end_index));
    this_sub_lag_p1_index = this_sub_lag_p1_index + lag_p1_start_index - 1;
    this_sub_lag_p1_time = single_onset_time(this_sub_lag_p1_index);

    % lag n1
    [~,this_sub_lag_n1_index] = min(this_sub_cz_fz_average(lag_n1_start_index:lag_n1_end_index));
    this_sub_lag_n1_index = this_sub_lag_n1_index + lag_n1_start_index - 1;
    this_sub_lag_n1_time = single_onset_time(this_sub_lag_n1_index);

    % lag p2
    [~,this_sub_lag_p2_index] = max(this_sub_cz_fz_average(lag_p2_start_index:lag_p2_end_index));
    this_sub_lag_p2_index = this_sub_lag_p2_index + lag_p2_start_index - 1;
    this_sub_lag_p2_time = single_onset_time(this_sub_lag_p2_index);


    % lag p3
    [~,this_sub_lag_p3_index] = max(this_sub_lag_bash_pz_average(lag_p3_start_index:lag_p3_end_index));
    this_sub_lag_p3_index = this_sub_lag_p3_index + lag_p3_start_index - 1;
    this_sub_lag_p3_time = single_onset_time(this_sub_lag_p3_index);

    figure(all_subs_fig_p1n1p2)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_cz_fz_average,'k')
    scatter(this_sub_lead_p1_time,this_sub_cz_fz_average(this_sub_lead_p1_index),'or','filled');
    scatter(this_sub_lead_n1_time,this_sub_cz_fz_average(this_sub_lead_n1_index),'ob','filled');
    scatter(this_sub_lead_p2_time,this_sub_cz_fz_average(this_sub_lead_p2_index),'oy','filled');
    scatter(this_sub_lag_p1_time,this_sub_cz_fz_average(this_sub_lag_p1_index),'og','filled');
    scatter(this_sub_lag_n1_time,this_sub_cz_fz_average(this_sub_lag_n1_index),'om','filled');
    scatter(this_sub_lag_p2_time,this_sub_cz_fz_average(this_sub_lag_p2_index),'ok','filled');

    figure(all_subs_fig_lead_p3)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_lead_bash_pz_average,'k')
    scatter(this_sub_lead_p3_time,this_sub_lead_bash_pz_average(this_sub_lead_p3_index),'or','filled');
    title('Lead P3 (bash in leading position ONLY)')

    figure(all_subs_fig_lag_p3)
    xlim([single_onset_time(1),single_onset_time(end)])
    hold on
    plot(single_onset_time,this_sub_lag_bash_pz_average,'k')
    scatter(this_sub_lag_p3_time,this_sub_lag_bash_pz_average(this_sub_lag_p3_index),'og','filled');
    title('lag P3 (bash in lagging position ONLY)')


    % Build data arrays to plot later
    % Num subjects x num conditions x num situations x num channels (32)
    % Condition order: small ITDs, large ITDs, small ILDs, large ILDs
    % ORGANIZATION:
    % (1) lead response, bash leads, target leads
    % (2) lead response, bash lags, target leads
    % (3) lead response, no bash, target leads
    % (4) lag response, bash leads, target leads
    % (5) lag response, bash lags, target leads
    % (6) lag response, no bash, target leads
    % (7) lead response, bash leads, masker leads
    % (8) lead response, bash lags, masker leads
    % (9) lead response, no bash, masker leads
    % (10) lag response, bash leads, masker leads
    % (11) lag response, bash lags, masker leads
    % (12) lag response, no bash, masker leads

    small_itd_cond = [3,7];
    large_itd_cond = [6,8];
    small_ild_cond = [1,4];
    large_ild_cond = [2,5];

    all_conditions = [small_itd_cond;large_itd_cond;small_ild_cond;large_ild_cond];
    condition_names = {'small_itd','large_itd','small_ild','large_ild'};

    peak_integration_time = 0.008; % s

    for ichannel = 1:32
    for icondition = 1:4
        these_conditions = all_conditions(icondition,:);
        for lead_stream = {'Target','Masker'}


            if lead_stream == "Target"
                % Bash happened in this stream, followed by non-bash,
                % want where responded = 1
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;
                % Non-bash happened in this stream, followed by bash,
                % want where responded = 0
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;

                % Non-bash happened in this stream, followed by
                % non-bash, want where responded = 0
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;

            elseif lead_stream == "Masker"
                % Bash happened in this stream, followed by non-bash,
                % want where responded = 0
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*ismember(ERP_info(:).Lead_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;
                % Non-bash happened in this stream, followed by bash,
                % want where responded = 1
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;

                % Non-bash happened in this stream, followed by
                % non-bash, want where responded = 0
                this_lead_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p1_index - round(peak_integration_time*fs):this_sub_lead_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_n1_index - round(peak_integration_time*fs):this_sub_lead_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p2_index - round(peak_integration_time*fs):this_sub_lead_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lead_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lead_p3_index - round(peak_integration_time*fs):this_sub_lead_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                this_lag_p1 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p1_index - round(peak_integration_time*fs):this_sub_lag_p1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_n1 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_n1_index - round(peak_integration_time*fs):this_sub_lag_n1_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p2 =  squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p2_index - round(peak_integration_time*fs):this_sub_lag_p2_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));
                this_lag_p3 = squeeze(mean(data_by_pair_onset_baselined(ichannel,this_sub_lag_p3_index - round(peak_integration_time*fs):this_sub_lag_p3_index + round(peak_integration_time*fs),logical(ismember(ERP_info(:).Condition,these_conditions)'.*ismember(ERP_info(:).Lead_Stream,lead_stream).*~ismember(ERP_info(:).Lead_Word,"bash").*~ismember(ERP_info(:).Lag_Word,"bash"))),[2,3],'omitnan'));

                all_subs_p1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p1,'Lag_Amplitude',this_lag_p1,'Electrode',electrode_names(ichannel));
                all_subs_n1(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_n1,'Lag_Amplitude',this_lag_n1,'Electrode',electrode_names(ichannel));
                all_subs_p2(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p2,'Lag_Amplitude',this_lag_p2,'Electrode',electrode_names(ichannel));
                all_subs_p3(structrow) = struct('S',subID,'Condition',condition_names(icondition),'Lead_Stream',lead_stream,'Lead_Word',"non-bash",'Lag_Word',"non-bash",'Lead_Amplitude',this_lead_p3,'Lag_Amplitude',this_lag_p3,'Electrode',electrode_names(ichannel));

                structrow = structrow + 1;


            end
        end
    end
    end

end

writetable(struct2table(all_subs_p1),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p1.csv')
writetable(struct2table(all_subs_n1),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_n1.csv')
writetable(struct2table(all_subs_p2),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p2.csv')
writetable(struct2table(all_subs_p3),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\all_subs_p3.csv')

