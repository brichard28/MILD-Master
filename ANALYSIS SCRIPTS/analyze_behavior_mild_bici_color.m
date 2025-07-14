% Benjamin Richardson
% Created: September 16th, 2024

% Script to analyze behavioral sensitivity (d-prime) for MILD MASTER
mild_master_root = 'C:\Users\benri\Documents\GitHub\MILD-Master\';

BehaviorTable = readtable(append(mild_master_root,'RESULTS DATA\MILD-BICI Behavior Files\mild-bici-color.xlsx'),'Format','auto');

subject_ID = char('mild_bici_color_1_snr=2');% mild bi-ci
num_conditions = 8;

all_hits = zeros(size(subject_ID,1),num_conditions);
all_lead_hits = zeros(size(subject_ID,1),num_conditions);
all_lag_hits = zeros(size(subject_ID,1),num_conditions);
all_hit_rts = zeros(size(subject_ID,1),num_conditions);

all_FAs = zeros(size(subject_ID,1),num_conditions);
all_lead_FAs = zeros(size(subject_ID,1),num_conditions);
all_lag_FAs = zeros(size(subject_ID,1),num_conditions);
all_FA_rts = zeros(size(subject_ID,1),num_conditions);

all_num_target_color = zeros(size(subject_ID,1),num_conditions);
all_num_masker_color = zeros(size(subject_ID,1),num_conditions);

all_num_lead_target_color = zeros(size(subject_ID,1),num_conditions);
all_num_lead_masker_color = zeros(size(subject_ID,1),num_conditions);

all_num_lag_target_color = zeros(size(subject_ID,1),num_conditions);
all_num_lag_masker_color = zeros(size(subject_ID,1),num_conditions);


color_words = {'red','white','green','blue'};

fs = 44100;
cue_dur = 2;
rt_fig = figure;

for isubject = 1:size(subject_ID,1) % For each subject...
    disp(string(subject_ID(isubject,:)))
    clicks_not_counted = 0;
    total_clicks = 0;

    % Load the word times for this subject
    WordTimesTable = readtable(append(mild_master_root,"RESULTS DATA\MILD-BICI Behavior Files\mild-bici-color__s_" + strtrim(string(subject_ID(isubject,:))) + "__Word_Times.csv"));

    run_count_per_condition = -1*ones(1,num_conditions); % array to keep track of which run in each condition we are on

    % Find the rows associated with this subject
    rows_this_subject = find(BehaviorTable.S == strtrim(string(subject_ID(isubject,:))));


%     all_maskers = {'side=r_itd=0_az=10_mag=0_lpf=0',...
%      'side=l_itd=0_az=10_mag=0_lpf=0',...
%      'side=r_itd=0_az=20_mag=0_lpf=0',...
%      'side=l_itd=0_az=20_mag=0_lpf=0',...
%      'side=r_itd=0_az=30_mag=0_lpf=0',...
%      'side=l_itd=0_az=30_mag=0_lpf=0',...
%      'side=r_itd=0_az=40_mag=0_lpf=0',...
%      'side=l_itd=0_az=40_mag=0_lpf=0',...
%      'side=r_itd=0_az=50_mag=0_lpf=0',...
%      'side=l_itd=0_az=50_mag=0_lpf=0',...
%      'side=r_itd=0_az=60_mag=0_lpf=0',...
%      'side=l_itd=0_az=60_mag=0_lpf=0',...
%      'side=r_itd=0_az=70_mag=0_lpf=0',...
%      'side=l_itd=0_az=70_mag=0_lpf=0',...
%      'side=r_itd=0_az=80_mag=0_lpf=0',...
%      'side=l_itd=0_az=80_mag=0_lpf=0',...
%      'side=r_itd=0_az=90_mag=0_lpf=0',...
%      'side=l_itd=0_az=90_mag=0_lpf=0',...
%       'side=r_itd=0_az=10_mag=1_lpf=0',...
%      'side=l_itd=0_az=10_mag=1_lpf=0',...
%      'side=r_itd=0_az=20_mag=1_lpf=0',...
%      'side=l_itd=0_az=20_mag=1_lpf=0',...
%      'side=r_itd=0_az=30_mag=1_lpf=0',...
%      'side=l_itd=0_az=30_mag=1_lpf=0',...
%      'side=r_itd=0_az=40_mag=1_lpf=0',...
%      'side=l_itd=0_az=40_mag=1_lpf=0',...
%      'side=r_itd=0_az=50_mag=1_lpf=0',...
%      'side=l_itd=0_az=50_mag=1_lpf=0',...
%      'side=r_itd=0_az=60_mag=1_lpf=0',...
%      'side=l_itd=0_az=60_mag=1_lpf=0',...
%      'side=r_itd=0_az=70_mag=1_lpf=0',...
%      'side=l_itd=0_az=70_mag=1_lpf=0',...
%      'side=r_itd=0_az=80_mag=1_lpf=0',...
%      'side=l_itd=0_az=80_mag=1_lpf=0',...
%      'side=r_itd=0_az=90_mag=1_lpf=0',...
%      'side=l_itd=0_az=90_mag=1_lpf=0'};

    all_maskers = {'side=r_itd=0_az=30_az_itd=0_mag=0',...
     'side=l_itd=0_az=30_az_itd=0_mag=0',...
     'side=r_itd=0_az=60_az_itd=0_mag=0',...
     'side=l_itd=0_az=60_az_itd=0_mag=0',...
     'side=r_itd=0_az=30_az_itd=0_mag=1',...
     'side=l_itd=0_az=30_az_itd=0_mag=1',...
     'side=r_itd=0_az=60_az_itd=0_mag=1',...
     'side=l_itd=0_az=60_az_itd=0_mag=1'};


    this_subject_color_click_distances = [];

    all_hit_rts_this_subject = struct('rts',[]);
    all_hit_rts_this_subject = repmat(all_hit_rts_this_subject,num_conditions,1);
    all_FA_rts_this_subject = struct('rts',[]);
    all_FA_rts_this_subject = repmat(all_FA_rts_this_subject,num_conditions,1);
    distances_to_nearest_color = [];

    % For each trial....
    for itrial = 1:length(rows_this_subject)


        %         if mod(itrial,10) == 0
        %             disp(itrial)
        %         end

        this_trial_condition = BehaviorTable.Condition(rows_this_subject(itrial)); % find the condition for this trial
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        this_trial_soundfile = BehaviorTable.Soundfile(rows_this_subject(itrial));
        this_trial_soundfile = extractAfter(string(this_trial_soundfile),append(string(subject_ID(isubject,:)),'/'));
        this_trial_soundfile = extractAfter(this_trial_soundfile,'/');


        run_count_per_condition(string(all_maskers) == string(this_trial_masker)) = run_count_per_condition(string(all_maskers) == string(this_trial_masker)) + 1;

        this_trial_run = run_count_per_condition(string(all_maskers) == string(this_trial_masker)); % find how many runs of this condition have happened already
        this_trial_click_times = table2array(BehaviorTable(rows_this_subject(itrial),8:end)); % find the click times for this trial
        this_trial_click_times(isnan(this_trial_click_times)) = []; % remove NaN from these click times
        this_trial_click_times = (this_trial_click_times/1000);
        % remove double clicks
        click_distances = diff(this_trial_click_times);
        click_distances_to_remove = find(click_distances < 0.2);
        this_trial_click_times(click_distances_to_remove + 1) = [];

        this_trial_click_times = this_trial_click_times + 0.702;

        which_soundfile_this_trial = BehaviorTable.Soundfile(rows_this_subject(itrial));
        which_soundfile_this_trial = cell2mat(which_soundfile_this_trial);
        slash_indices = find(which_soundfile_this_trial == '/');
        slash_index = max(slash_indices);
        which_soundfile_this_trial = which_soundfile_this_trial(slash_index + 1:slash_index + 3);
        if contains(which_soundfile_this_trial, '__')
            which_soundfile_this_trial = str2num(which_soundfile_this_trial(1));
        elseif contains(which_soundfile_this_trial,'_')
            which_soundfile_this_trial = str2num(which_soundfile_this_trial(1:2));
        else
            which_soundfile_this_trial = str2num(which_soundfile_this_trial);
        end
        %disp(which_soundfile_this_trial)

        
        this_trial_target_all = WordTimesTable(string(WordTimesTable.Var2) == this_trial_soundfile & string(WordTimesTable.Var3) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times =  (table2array(this_trial_target_all(:,2:2:end)) ./ fs) + cue_dur; % have to add cue back in

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Var2) == this_trial_soundfile & string(WordTimesTable.Var3) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times =  (table2array(this_trial_target_all(:,2:2:end)) ./ fs) + cue_dur;


        this_trial_whether_target_lead = sign(this_trial_masker_times - this_trial_target_times);
        this_trial_whether_target_lead(this_trial_whether_target_lead == -1) = 0; % 0 when masker leads, 1 when target leads


        target_color_lead = this_trial_whether_target_lead(ismember(this_trial_target_words,color_words)) ; % 0 when color lags (second position within pair), 1 when color leads (first position within pair)
        masker_color_lead = 1 - this_trial_whether_target_lead(ismember(this_trial_masker_words,color_words));

        % Store number of color words in the target and masker
        all_num_target_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_target_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words));
        all_num_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words));

        all_num_lead_target_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_target_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words) & this_trial_whether_target_lead == 1);
        all_num_lag_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words) & this_trial_whether_target_lead == 1);

        all_num_lag_target_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_target_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words) & this_trial_whether_target_lead == 0);
        all_num_lead_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words) & this_trial_whether_target_lead == 0);

        % Find just color times in target and masker
        this_trial_target_color_times = this_trial_target_times(ismember(this_trial_target_words, color_words));
        this_trial_masker_color_times = this_trial_masker_times(ismember(this_trial_masker_words, color_words));

        %% Hit and False Alarm Windows

        threshold_window_start = 0.4; %0.2
        threshold_window_end =  1.5; % 1.0
        tVec = 0:1/44100:16;
        hit_windows = zeros(1,length(tVec)); % create an empty array to define hit windows
        lead_hit_windows = zeros(1,length(tVec));
        lag_hit_windows = zeros(1,length(tVec));
        FA_windows = zeros(1,length(tVec)); % create an empty array to define false alarm windows
        lead_FA_windows = zeros(1,length(tVec));
        lag_FA_windows = zeros(1,length(tVec));

        % specify hit windows
        for i = 1:length(this_trial_target_color_times) % for each of the current target color times...
            [~,start_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_start))); % ...the hit window will start threshold_window_start seconds after the word onset
            [~,end_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_end))); % ...the hit window will end threshold_window_end seconds after the word onset

            hit_windows(start_index_hit_window:end_index_hit_window) = 1; % a value of 1 in the vector hit_windows indicate an area where, if a click falls, it will be counted as a hit
            if target_color_lead(i) == 0 % then target LAGs on this trial
                lag_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            elseif target_color_lead(i) == 1 % then target LEADS on this trial
                lead_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            end

        end

        % specify false alarm windows
        for i = 1:length(this_trial_masker_color_times) % for each of the current masker times...
            [~,start_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_start))); % ...the false alarm window will start threshold_window_start seconds after the word onset
            [~,end_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_end))); % ...the false alarm window will end threshold_window_end seconds after the word onset

            if any(hit_windows(start_index_FA_window:end_index_FA_window) == 1)
                continue
            else
                FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                if masker_color_lead(i) == 0 % then target LAGs on this trial
                    lag_FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                elseif masker_color_lead(i) == 1 % then target LEADS on this trial
                    lead_FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                end
            end
        end

        test_vector(itrial,:) = FA_windows + hit_windows;

        % ...Calculate the hit rate, FA rate in this trial
        for iclick = 1:length(this_trial_click_times)
            [~,current_click_index] = min(abs(tVec - this_trial_click_times(iclick))); % ...find the time index of that click...

            if hit_windows(current_click_index) == 1 % ...if that click falls within a hit window...
                all_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;

                % populate lead and lag
                if lead_hit_windows(current_click_index) == 1
                    all_lead_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_hit_windows(current_click_index) == 1
                    all_lag_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end

                % this reaction time
                distances_to_color = [];
                for icolortime = 1:length(this_trial_target_color_times)
                    distances_to_color = [distances_to_color, this_trial_click_times(iclick) - this_trial_target_color_times(icolortime)];
                end
                distances_to_color(distances_to_color < 0) =[];
                all_hit_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts = [all_hit_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts, min(distances_to_color)];

            elseif FA_windows(current_click_index) == 1 %...otherwise if that click falls within a false alarm window...
                all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                % populate lead and lag
                if lead_FA_windows(current_click_index) == 1
                    all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_FA_windows(current_click_index) == 1
                    all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end


                % this reaction time
                distances_to_color = [];
                for icolortime = 1:length(this_trial_masker_color_times)
                    distances_to_color = [distances_to_color, this_trial_click_times(iclick) - this_trial_masker_color_times(icolortime)];
                end
                distances_to_color(distances_to_color < 0) =[];
                all_FA_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts = [all_FA_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts, min(distances_to_color)];

            else% ...if the click is not counted as either
                clicks_not_counted = clicks_not_counted + 1;
            end
            total_clicks = total_clicks + 1;
        end

        %% Calculate time difference between each click and each color time
        distances_click_to_target_color= [];
        for icolortime = 1:length(this_trial_target_color_times)
            distances_click_to_target_color(icolortime,:) = this_trial_click_times - this_trial_target_color_times(icolortime);
        end
        distances_click_to_target_color(distances_click_to_target_color < 0) = nan;

        %% Find the nearest color time to each click (minimum positive value of click_distances in each column)


        [~,nearest_click] = min(abs(distances_click_to_target_color),[],1); % find the nearest click to each target word
        for i = 1:length(this_trial_click_times)
            if isnan(distances_click_to_target_color(:,i)) == ones(1,length(this_trial_target_color_times)) % all of these clicks were before the first word
                nearest_click(i) = nan;
            else

                this_subject_color_click_distances = [this_subject_color_click_distances, distances_click_to_target_color(nearest_click(i),i)];
            end

        end


    end
    save(append(string(subject_ID(isubject,:)),'color_click_distances.mat'),'this_subject_color_click_distances');
    figure(rt_fig)
    subplot(round(size(subject_ID,1)),2,isubject)
    histogram(this_subject_color_click_distances,100)
    %xlim([0,3])

    [counts,edges] = histcounts(this_subject_color_click_distances,100);
    %[peakValue, peakIndex] = findpeaks(counts); % Find peak value and index
    [peakValue, peakIndex] = max(counts);
    peakXValue = edges(peakIndex); % Get the x-axis value corresponding to the peak



    disp(append(string(subject_ID(isubject,:)),' most frequent color reaction time: ',num2str(peakXValue*1000),' ms '))
    disp(append(string(subject_ID(isubject,:)),': ', num2str((clicks_not_counted/total_clicks)*100), '% of clicks not counted'))
end

attend_right_indices = 1:2:length(all_maskers);
attend_left_indices = 2:2:length(all_maskers);

all_hit_rates = all_hits./all_num_target_color;
all_lead_hit_rates = all_lead_hits./all_num_lead_target_color;
all_lag_hit_rates = all_lag_hits./all_num_lag_target_color;
all_FA_rates = all_FAs./all_num_masker_color;
all_lead_FA_rates = all_lead_FAs./all_num_lead_masker_color;
all_lag_FA_rates = all_lag_FAs./all_num_lag_masker_color;

all_hit_rates_collapsed_left_and_right = sum([all_hits(attend_right_indices);all_hits(attend_left_indices)],1)./sum([all_num_target_color(attend_right_indices);all_num_target_color(attend_left_indices)],1);
all_FA_rates_collapsed_left_and_right = sum([all_FAs(attend_right_indices);all_FAs(attend_left_indices)],1)./sum([all_num_masker_color(attend_right_indices);all_num_masker_color(attend_left_indices)],1);

all_hit_rates_collapsed_left_and_right(all_hit_rates_collapsed_left_and_right  == 0) = 0.001;
all_hit_rates_collapsed_left_and_right (all_hit_rates_collapsed_left_and_right  >= 1) = 0.999;
all_FA_rates_collapsed_left_and_right(all_FA_rates_collapsed_left_and_right  == 0) = 0.001;
all_FA_rates_collapsed_left_and_right (all_FA_rates_collapsed_left_and_right  >= 1) = 0.999;

all_d_primes_collapsed_left_and_right = norminv(all_hit_rates_collapsed_left_and_right) - norminv(all_FA_rates_collapsed_left_and_right);
figure;
subplot(1,2,1)
hold on
scatter(1:length(all_maskers)/4,all_hit_rates_collapsed_left_and_right(1:length(all_maskers)/4),'r','filled')
scatter(1:length(all_maskers)/4,all_hit_rates_collapsed_left_and_right(length(all_maskers)/4 + 1:end),'b','filled')
xticks(1:length(all_maskers)/4);
xticklabels(10:10:90)
ylabel('Hit Rate','FontSize',18)
ylim([0,1])
xlabel('Azimuth (deg.)','FontSize',18)

subplot(1,2,2)
hold on
scatter(1:length(all_maskers)/4,all_FA_rates_collapsed_left_and_right(1:length(all_maskers)/4),'r','filled')
scatter(1:length(all_maskers)/4,all_FA_rates_collapsed_left_and_right(length(all_maskers)/4 + 1:end),'b','filled')
xticks(1:length(all_maskers)/4);
xticklabels(10:10:90)
ylabel('FA Rate','FontSize',18)
legend({'Natural ILDs','Magnified ILDs'},'FontSize',18)
ylim([0,1])
xlabel('Azimuth (deg.)','FontSize',18)

figure;
hold on
scatter(1:length(all_maskers)/4,all_d_primes_collapsed_left_and_right(1:length(all_maskers)/4),'r','filled')
scatter(1:length(all_maskers)/4,all_d_primes_collapsed_left_and_right(length(all_maskers)/4 + 1:end),'b','filled')
xticks(1:length(all_maskers)/4);
xticklabels(10:10:90)
ylabel("d'",'FontSize',18)
ylim([0,3.5])
xlabel('Azimuth (deg.)','FontSize',18)


% %% NEW ORDER = itd50 noise, itd500 noise, ildnat noise, ild10 noise, itd50 speech, itd500 speech, ildnat speech, ild10 speech
% 
% 
% all_hits_collapsed_left_and_right = [];
% all_hits_collapsed_left_and_right(1,:) = sum(all_hits(:,small_itd_cond),2); % itd50 / ild10
% all_hits_collapsed_left_and_right(2,:) = sum(all_hits(:,large_itd_cond),2); % itd100 / ild10mag
% all_hits_collapsed_left_and_right(3,:) = sum(all_hits(:,small_ild_cond),2); % itd200 /ild60
% all_hits_collapsed_left_and_right(4,:) = sum(all_hits(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_lead_hits_collapsed_left_and_right = [];
% all_lead_hits_collapsed_left_and_right(1,:) = sum(all_lead_hits(:,small_itd_cond),2); % itd50 / ild10
% all_lead_hits_collapsed_left_and_right(2,:) = sum(all_lead_hits(:,large_itd_cond),2); % itd100 / ild10mag
% all_lead_hits_collapsed_left_and_right(3,:) = sum(all_lead_hits(:,small_ild_cond),2); % itd200 /ild60
% all_lead_hits_collapsed_left_and_right(4,:) = sum(all_lead_hits(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_lag_hits_collapsed_left_and_right = [];
% all_lag_hits_collapsed_left_and_right(1,:) = sum(all_lag_hits(:,small_itd_cond),2); % itd50 / ild10
% all_lag_hits_collapsed_left_and_right(2,:) = sum(all_lag_hits(:,large_itd_cond),2); % itd100 / ild10mag
% all_lag_hits_collapsed_left_and_right(3,:) = sum(all_lag_hits(:,small_ild_cond),2); % itd200 /ild60
% all_lag_hits_collapsed_left_and_right(4,:) = sum(all_lag_hits(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_FAs_collapsed_left_and_right = [];
% all_FAs_collapsed_left_and_right(1,:) = sum(all_FAs(:,small_itd_cond),2); % itd5
% all_FAs_collapsed_left_and_right(2,:) = sum(all_FAs(:,large_itd_cond),2); % itd15
% all_FAs_collapsed_left_and_right(3,:) = sum(all_FAs(:,small_ild_cond),2); % ild5
% all_FAs_collapsed_left_and_right(4,:) = sum(all_FAs(:,large_ild_cond),2); % ild15
% 
% all_lead_FAs_collapsed_left_and_right = [];
% all_lead_FAs_collapsed_left_and_right(1,:) = sum(all_lead_FAs(:,small_itd_cond),2); % itd50 / ild10
% all_lead_FAs_collapsed_left_and_right(2,:) = sum(all_lead_FAs(:,large_itd_cond),2); % itd100 / ild10mag
% all_lead_FAs_collapsed_left_and_right(3,:) = sum(all_lead_FAs(:,small_ild_cond),2); % itd200 /ild60
% all_lead_FAs_collapsed_left_and_right(4,:) = sum(all_lead_FAs(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_lag_FAs_collapsed_left_and_right = [];
% 
% all_lag_FAs_collapsed_left_and_right(1,:) = sum(all_lag_FAs(:,small_itd_cond),2); % itd50 / ild10
% all_lag_FAs_collapsed_left_and_right(2,:) = sum(all_lag_FAs(:,large_itd_cond),2); % itd100 / ild10mag
% all_lag_FAs_collapsed_left_and_right(3,:) = sum(all_lag_FAs(:,small_ild_cond),2); % itd200 /ild60
% all_lag_FAs_collapsed_left_and_right(4,:) = sum(all_lag_FAs(:,large_ild_cond),2); % itd400 / ild60mag
% 
% 
% all_num_target_color_collapsed_left_and_right = [];
% all_num_target_color_collapsed_left_and_right(1,:) = sum(all_num_target_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_target_color_collapsed_left_and_right(2,:) = sum(all_num_target_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_target_color_collapsed_left_and_right(3,:) = sum(all_num_target_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_target_color_collapsed_left_and_right(4,:) = sum(all_num_target_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_num_lead_target_color_collapsed_left_and_right = [];
% all_num_lead_target_color_collapsed_left_and_right(1,:) = sum(all_num_lead_target_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_lead_target_color_collapsed_left_and_right(2,:) = sum(all_num_lead_target_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_lead_target_color_collapsed_left_and_right(3,:) = sum(all_num_lead_target_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_lead_target_color_collapsed_left_and_right(4,:) = sum(all_num_lead_target_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_num_lag_target_color_collapsed_left_and_right = [];
% all_num_lag_target_color_collapsed_left_and_right(1,:) = sum(all_num_lag_target_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_lag_target_color_collapsed_left_and_right(2,:) = sum(all_num_lag_target_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_lag_target_color_collapsed_left_and_right(3,:) = sum(all_num_lag_target_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_lag_target_color_collapsed_left_and_right(4,:) = sum(all_num_lag_target_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_num_masker_color_collapsed_left_and_right = [];
% all_num_masker_color_collapsed_left_and_right(1,:) = sum(all_num_masker_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_masker_color_collapsed_left_and_right(2,:) = sum(all_num_masker_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_masker_color_collapsed_left_and_right(3,:) = sum(all_num_masker_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_masker_color_collapsed_left_and_right(4,:) = sum(all_num_masker_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_num_lead_masker_color_collapsed_left_and_right = [];
% all_num_lead_masker_color_collapsed_left_and_right(1,:) = sum(all_num_lead_masker_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_lead_masker_color_collapsed_left_and_right(2,:) = sum(all_num_lead_masker_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_lead_masker_color_collapsed_left_and_right(3,:) = sum(all_num_lead_masker_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_lead_masker_color_collapsed_left_and_right(4,:) = sum(all_num_lead_masker_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% all_num_lag_masker_color_collapsed_left_and_right = [];
% all_num_lag_masker_color_collapsed_left_and_right(1,:) = sum(all_num_lag_masker_color(:,small_itd_cond),2); % itd50 / ild10
% all_num_lag_masker_color_collapsed_left_and_right(2,:) = sum(all_num_lag_masker_color(:,large_itd_cond),2); % itd100 / ild10mag
% all_num_lag_masker_color_collapsed_left_and_right(3,:) = sum(all_num_lag_masker_color(:,small_ild_cond),2); % itd200 /ild60
% all_num_lag_masker_color_collapsed_left_and_right(4,:) = sum(all_num_lag_masker_color(:,large_ild_cond),2); % itd400 / ild60mag
% 
% 
% %% Calculate Rates
% % Hit Rates

% 
% all_hit_rates_collapsed = all_hits_collapsed_left_and_right./all_num_target_color_collapsed_left_and_right;
% all_lead_hit_rates_collapsed = all_lead_hits_collapsed_left_and_right./all_num_lead_target_color_collapsed_left_and_right;
% all_lag_hit_rates_collapsed = all_lag_hits_collapsed_left_and_right./all_num_lag_target_color_collapsed_left_and_right;
% 
% all_hit_rates(all_hit_rates == 0) = 0.001;
% all_hit_rates(all_hit_rates >= 1) = 0.999;
% all_lead_hit_rates(all_lead_hit_rates == 0) = 0.001;
% all_lead_hit_rates(all_lead_hit_rates >= 1) = 0.999;
% all_lag_hit_rates(all_lag_hit_rates == 0) = 0.001;
% all_lag_hit_rates(all_lag_hit_rates >= 1) = 0.999;
% 
% all_hit_rates_collapsed(all_hit_rates_collapsed == 0) = 0.001;
% all_hit_rates_collapsed(all_hit_rates_collapsed >= 1) = 0.999;
% all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed == 0) = 0.001;
% all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed >= 1) = 0.999;
% all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed == 0) = 0.001;
% all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed >= 1) = 0.999;
% 
% % FA Rates
% all_FA_rates = all_FAs./all_num_masker_color;
% all_lead_FA_rates = all_lead_FAs./all_num_lead_masker_color;
% all_lag_FA_rates = all_lag_FAs./all_num_lag_masker_color;
% 
% all_FA_rates_collapsed = all_FAs_collapsed_left_and_right./all_num_masker_color_collapsed_left_and_right;
% all_lead_FA_rates_collapsed = all_lead_FAs_collapsed_left_and_right./all_num_lead_masker_color_collapsed_left_and_right;
% all_lag_FA_rates_collapsed = all_lag_FAs_collapsed_left_and_right./all_num_lag_masker_color_collapsed_left_and_right;
% 
% all_FA_rates(all_FA_rates == 0) = 0.001;
% all_FA_rates(all_FA_rates >= 1) = 0.999;
% all_lead_FA_rates(all_lead_FA_rates == 0) = 0.001;
% all_lead_FA_rates(all_lead_FA_rates >= 1) = 0.999;
% all_lag_FA_rates(all_lag_FA_rates == 0) = 0.001;
% all_lag_FA_rates(all_lag_FA_rates >= 1) = 0.999;
% 
% all_FA_rates_collapsed(all_FA_rates_collapsed == 0) = 0.001;
% all_FA_rates_collapsed(all_FA_rates_collapsed >= 1) = 0.999;
% all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed == 0) = 0.001;
% all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed >= 1) = 0.999;
% all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed == 0) = 0.001;
% all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed >= 1) = 0.999;
% 
% 
% %% D-prime calculation
% d_primes_all = norminv(all_hit_rates) - norminv(all_FA_rates);
% d_primes_collapsed = [];
% d_primes_collapsed(1,:) = mean(d_primes_all(:,small_itd_cond),2); % itd50 / ild10
% d_primes_collapsed(2,:) = mean(d_primes_all(:,large_itd_cond),2); % itd100 / ild10mag
% d_primes_collapsed(3,:) = mean(d_primes_all(:,small_ild_cond),2); % itd200 /ild60
% d_primes_collapsed(4,:) = mean(d_primes_all(:,large_ild_cond),2); % itd400 / ild60mag
% 
% lead_d_primes_all = norminv(all_lead_hit_rates) - norminv(all_lead_FA_rates);
% lead_d_primes_collapsed = [];
% lead_d_primes_collapsed(1,:) = mean(lead_d_primes_all(:,small_itd_cond),2); % itlead_d50 / illead_d10
% lead_d_primes_collapsed(2,:) = mean(lead_d_primes_all(:,large_itd_cond),2); % itlead_d100 / illead_d10mag
% lead_d_primes_collapsed(3,:) = mean(lead_d_primes_all(:,small_ild_cond),2); % itlead_d200 /illead_d60
% lead_d_primes_collapsed(4,:) = mean(lead_d_primes_all(:,large_ild_cond),2); % itlead_d400 / illead_d60mag
% 
% lag_d_primes_all = norminv(all_lag_hit_rates) - norminv(all_lag_FA_rates);
% lag_d_primes_collapsed = [];
% lag_d_primes_collapsed(1,:) = mean(lag_d_primes_all(:,small_itd_cond),2); % itlag_d50 / illag_d10
% lag_d_primes_collapsed(2,:) = mean(lag_d_primes_all(:,large_itd_cond),2); % itlag_d100 / illag_d10mag
% lag_d_primes_collapsed(3,:) = mean(lag_d_primes_all(:,small_ild_cond),2); % itlag_d200 /illag_d60
% lag_d_primes_collapsed(4,:) = mean(lag_d_primes_all(:,large_ild_cond),2); % itlag_d400 / illag_d60mag
% 
% %% PLOT ALL
% figure;
% colors_to_use = ["g","m","g","m"];
% 
% subplot(1,3,1)
% % plot(1:2,d_primes_collapsed(1:2,d_primes_collapsed(1,:) > d_primes_collapsed(2,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(1:2,d_primes_collapsed(1:2,d_primes_collapsed(1,:) < d_primes_collapsed(2,:)),'-','Color',colors_to_use(2));
% % 
% % plot(3:4,d_primes_collapsed(3:4,d_primes_collapsed(3,:) > d_primes_collapsed(4,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(3:4,d_primes_collapsed(3:4,d_primes_collapsed(3,:) < d_primes_collapsed(4,:)),'-','Color',colors_to_use(2));
% 
% 
% 
% for iplot = 1:4
%     errorbar(iplot,mean(d_primes_collapsed(iplot,:),2),std(d_primes_collapsed(iplot,:),[],2)/sqrt(size(subject_ID,1) - 1),'o','Color',colors_to_use(iplot),'LineWidth',3)
%     hold on
% end
% ylabel("d'",'FontSize',18)
% %xticklabels({'50','100','200','400'})
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% 
% subplot(1,3,2)
% % plot(1:2,all_hit_rates_collapsed(1:2,all_hit_rates_collapsed(1,:) > all_hit_rates_collapsed(2,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(1:2,all_hit_rates_collapsed(1:2,all_hit_rates_collapsed(1,:) < all_hit_rates_collapsed(2,:)),'-','Color',colors_to_use(2));
% % 
% % plot(3:4,all_hit_rates_collapsed(3:4,all_hit_rates_collapsed(3,:) > all_hit_rates_collapsed(4,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(3:4,all_hit_rates_collapsed(3:4,all_hit_rates_collapsed(3,:) < all_hit_rates_collapsed(4,:)),'-','Color',colors_to_use(2));
% % 
% 
% 
% for iplot = 1:4
%     errorbar(iplot,mean(all_hit_rates_collapsed(iplot,:),2),std(all_hit_rates_collapsed(iplot,:),[],2)/sqrt(size(subject_ID,1) - 1),'o','Color',colors_to_use(iplot),'LineWidth',3)
%     hold on
% end
% ylabel("Hit Rate",'FontSize',18)
% %xticklabels({'50','100','200','400'})
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% ylim([0,1])
% 
% subplot(1,3,3)
% % plot(1:2,all_FA_rates_collapsed(1:2,all_FA_rates_collapsed(1,:) > all_FA_rates_collapsed(2,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(1:2,all_FA_rates_collapsed(1:2,all_FA_rates_collapsed(1,:) < all_FA_rates_collapsed(2,:)),'-','Color',colors_to_use(2));
% % 
% % plot(3:4,all_FA_rates_collapsed(3:4,all_FA_rates_collapsed(3,:) > all_FA_rates_collapsed(4,:)),'-','Color',colors_to_use(1));
% % hold on
% % plot(3:4,all_FA_rates_collapsed(3:4,all_FA_rates_collapsed(3,:) < all_FA_rates_collapsed(4,:)),'-','Color',colors_to_use(2));
% 
% 
% 
% for iplot = 1:4
%     errorbar(iplot,mean(all_FA_rates_collapsed(iplot,:),2),std(all_FA_rates_collapsed(iplot,:),[],2)/sqrt(size(subject_ID,1) - 1),'o','Color',colors_to_use(iplot),'LineWidth',3)
%     hold on
% end
% ylabel("False Alarm Rate",'FontSize',18)
% %xticklabels({'50','100','200','400'})
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% ylim([0,1])
% 
% %% PLOT BROKEN UP BY LEAD LAG
% % Lead vs. lag d-prime
% x_offset_lead = -0.1;
% x_offset_lag = 0.1;
% figure;
% subplot(1,3,1)
% hold on
% p1 = errorbar((1:2:4) + x_offset_lead,mean(lead_d_primes_collapsed(1:2:4,:),2),std(lead_d_primes_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% p2 = errorbar((2:2:4) + x_offset_lead,mean(lead_d_primes_collapsed(2:2:4,:),2),std(lead_d_primes_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% xticks(1:4)
% ylim([0 7])
% ylabel("d'",'FontSize',18)
% 
% hold on
% p3 = errorbar((1:2:4) + x_offset_lag,mean(lag_d_primes_collapsed(1:2:4,:),2),std(lag_d_primes_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% p4 = errorbar((2:2:4) + x_offset_lag,mean(lag_d_primes_collapsed(2:2:4,:),2),std(lag_d_primes_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% xticks(1:4)
% ylim([0,7])
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% legend([p1(1), p3(1)],{'Lead','Lag'})
% 
% % lead vs. lag hit rate
% subplot(1,3,2)
% hold on
% p1 = errorbar((1:2:4) + x_offset_lead,mean(all_lead_hit_rates_collapsed(1:2:4,:),2),std(all_lead_hit_rates_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% p2 = errorbar((2:2:4) + x_offset_lead,mean(all_lead_hit_rates_collapsed(2:2:4,:),2),std(all_lead_hit_rates_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% ylim([0 1])
% ylabel("Hit Rate",'FontSize',18)
% xticks(1:4)
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% legend([p1(1), p3(1)],{'Lead','Lag'})
% 
% hold on
% p3 = errorbar((1:2:4) + x_offset_lag,mean(all_lag_hit_rates_collapsed(1:2:4,:),2),std(all_lag_hit_rates_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% p4 = errorbar((2:2:4) + x_offset_lag,mean(all_lag_hit_rates_collapsed(2:2:4,:),2),std(all_lag_hit_rates_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% ylim([0 1])
% xticks(1:4)
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% legend([p1(1), p3(1)],{'Lead','Lag'})
% 
% % lead vs. lag FA rate
% subplot(1,3,3)
% hold on
% p1 = errorbar((1:2:4) + x_offset_lead,mean(all_lead_FA_rates_collapsed(1:2:4,:),2),std(all_lead_FA_rates_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% p2 = errorbar((2:2:4) + x_offset_lead,mean(all_lead_FA_rates_collapsed(2:2:4,:),2),std(all_lead_FA_rates_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','c','MarkerFaceColor','c');
% ylim([0 1])
% ylabel("FA Rate",'FontSize',18)
% xticks(1:4)
% 
% hold on
% p3 = errorbar((1:2:4) + x_offset_lag,mean(all_lag_FA_rates_collapsed(1:2:4,:),2),std(all_lag_FA_rates_collapsed(1:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% p4 = errorbar((2:2:4) + x_offset_lag,mean(all_lag_FA_rates_collapsed(2:2:4,:),2),std(all_lag_FA_rates_collapsed(2:2:4,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','m','MarkerFaceColor','m');
% ylim([0 1])
% xticks(1:4)
% xticklabels({'60 deg ILDs','60 deg ILDs Mag','90 deg ILDs','90 deg ILDs Mag'})
% legend([p1(1), p3(1)],{'Lead','Lag'})


%% Save data
%save('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Behavior_Results.mat','d_primes_collapsed','d_primes_collapsed','all_hit_rates_collapsed','all_FA_rates_collapsed')

hit_rate_table = array2table(all_hit_rates_collapsed);
writetable(rows2vars(hit_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Hit_Rates.csv')

FA_rate_table = array2table(all_FA_rates_collapsed(5:end,:));
writetable(rows2vars(FA_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_FA_Rates.csv')

d_prime_table = array2table(d_primes_collapsed(5:end,:));
writetable(rows2vars(d_prime_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_d_primes.csv')