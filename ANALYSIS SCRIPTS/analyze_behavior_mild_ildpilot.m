% Benjamin Richardson
% Created: September 16th, 2024

% Script to analyze behavioral sensitivity (d-prime) for MILD MASTER

BehaviorTable = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx','Format','auto');

subject_ID = char('ildpilot3','ildpilot4','ildpilot5','ildpilot6','ildpilot7'); % ild pilot % 'ildpilot1','ildpilot2'
num_conditions = 16;

all_hits = zeros(size(subject_ID,1),num_conditions);
all_lead_hits = zeros(size(subject_ID,1),num_conditions);
all_lag_hits = zeros(size(subject_ID,1),num_conditions);
all_hit_rts = struct();
all_FAs = zeros(size(subject_ID,1),num_conditions);
all_lead_FAs = zeros(size(subject_ID,1),num_conditions);
all_lag_FAs = zeros(size(subject_ID,1),num_conditions);
all_FA_rts = struct();

all_num_target_bash = zeros(size(subject_ID,1),num_conditions);
all_num_masker_bash = zeros(size(subject_ID,1),num_conditions);

all_num_lead_target_bash = zeros(size(subject_ID,1),num_conditions);
all_num_lead_masker_bash = zeros(size(subject_ID,1),num_conditions);

all_num_lag_target_bash = zeros(size(subject_ID,1),num_conditions);
all_num_lag_masker_bash = zeros(size(subject_ID,1),num_conditions);

all_maskers = {'side=r_itd=0_az=5_mag=0',...
'side=l_itd=0_az=5_mag=0',...
'side=r_itd=0_az=5_mag=1',...
'side=l_itd=0_az=5_mag=1',...
'side=r_itd=0_az=10_mag=0',...
'side=l_itd=0_az=10_mag=0',...
'side=r_itd=0_az=10_mag=1',...
'side=l_itd=0_az=10_mag=1',...
'side=r_itd=0_az=20_mag=0',...
'side=l_itd=0_az=20_mag=0',...
'side=r_itd=0_az=20_mag=1',...
'side=l_itd=0_az=20_mag=1',...
'side=r_itd=0_az=30_mag=0',...
'side=l_itd=0_az=30_mag=0',...
'side=r_itd=0_az=30_mag=1',...
'side=l_itd=0_az=30_mag=1'};

fs = 44100;
cue_dur = 1.5;

for isubject = 1:size(subject_ID,1) % For each subject...
    disp(string(subject_ID(isubject,:)))
    clicks_not_counted = 0;
    total_clicks = 0;

    % Load the word times for this subject
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(string(subject_ID(isubject,:))) + "__Word_Times.csv");

    run_count_per_condition = -1*ones(1,num_conditions); % array to keep track of which run in each condition we are on

    % Find the rows associated with this subject
    rows_this_subject = find(BehaviorTable.S == strtrim(string(subject_ID(isubject,:))));


    % For each trial....
    for itrial = 1:length(rows_this_subject)

        this_trial_condition = BehaviorTable.Condition(rows_this_subject(itrial)); % find the condition for this trial
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        run_count_per_condition(string(all_maskers) == string(this_trial_masker)) = run_count_per_condition(string(all_maskers) == string(this_trial_masker)) + 1;

        this_trial_run = run_count_per_condition(string(all_maskers) == string(this_trial_masker)); % find how many runs of this condition have happened already
        this_trial_click_times = table2array(BehaviorTable(rows_this_subject(itrial),8:end)); % find the click times for this trial
        this_trial_click_times(isnan(this_trial_click_times)) = []; % remove NaN from these click times
        this_trial_click_times = this_trial_click_times/1000;
        % remove double clicks
        click_distances = diff(this_trial_click_times);
        click_distances_to_remove = find(click_distances < 0.2);
        this_trial_click_times(click_distances_to_remove + 1) = [];

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


        this_trial_target_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = (table2array(this_trial_target_all(:,2:2:end)) ./ fs) + cue_dur; % have to add cue back in

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Var1) == append(string(which_soundfile_this_trial),"_",this_trial_masker) & string(WordTimesTable.Var3) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = (table2array(this_trial_masker_all(:,2:2:end)) ./ fs) + cue_dur;

%         if this_trial_target_words(1) == "bash"
%             this_trial_target_words(1) = [];
%             this_trial_masker_words(1) = [];
%             this_trial_target_times(1) = [];
%             this_trial_masker_times(1) = [];
%         elseif this_trial_masker_words(1) == "bash"
%             this_trial_target_words(1) = [];
%             this_trial_masker_words(1) = [];
%             this_trial_target_times(1) = [];
%             this_trial_masker_times(1) = [];
%         end

        this_trial_whether_target_lead = sign(this_trial_masker_times - this_trial_target_times);
        this_trial_whether_target_lead(this_trial_whether_target_lead == -1) = 0; % 0 when masker leads, 1 when target leads

        
        
        % Find just color times in target and masker
        this_trial_target_bash_times = this_trial_target_times(this_trial_target_words == "bash");
        this_trial_masker_bash_times = this_trial_masker_times(this_trial_masker_words == "bash");

        target_bash_lead = this_trial_whether_target_lead(this_trial_target_words == "bash") ; % 0 when bash lags (second position within pair), 1 when bash leads (first position within pair)
        masker_bash_lead = 1 - this_trial_whether_target_lead(this_trial_masker_words == "bash");

        % Store number of color words in the target and masker
        all_num_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_target_words == "bash");
        all_num_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_masker_words == "bash");
        
        all_num_lead_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_target_words == "bash" & this_trial_whether_target_lead == 1);
        all_num_lag_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_masker_words == "bash" & this_trial_whether_target_lead == 1);
        
        all_num_lag_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_target_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_target_words == "bash" & this_trial_whether_target_lead == 0);
        all_num_lead_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_masker_bash(isubject,string(all_maskers) == string(this_trial_masker)) + sum(this_trial_masker_words == "bash" & this_trial_whether_target_lead == 0);
        
       %% Hit and False Alarm Windows

       threshold_window_start = 0.2; %0.2
       threshold_window_end =  1.3; % 1.0
       tVec = 0:1/44100:16;
       hit_windows = zeros(1,length(tVec)); % create an empty array to define hit windows
       lead_hit_windows = zeros(1,length(tVec));
       lag_hit_windows = zeros(1,length(tVec));
       FA_windows = zeros(1,length(tVec)); % create an empty array to define false alarm windows
       lead_FA_windows = zeros(1,length(tVec));
       lag_FA_windows = zeros(1,length(tVec));

        % specify hit windows
        for i = 1:length(this_trial_target_bash_times) % for each of the current target color times...
            [~,start_index_hit_window] = min(abs(tVec - (this_trial_target_bash_times(i)+threshold_window_start))); % ...the hit window will start threshold_window_start seconds after the word onset
            [~,end_index_hit_window] = min(abs(tVec - (this_trial_target_bash_times(i)+threshold_window_end))); % ...the hit window will end threshold_window_end seconds after the word onset

            hit_windows(start_index_hit_window:end_index_hit_window) = 1; % a value of 1 in the vector hit_windows indicate an area where, if a click falls, it will be counted as a hit
            if target_bash_lead(i) == 0 % then target LAGs on this trial
                lag_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            elseif target_bash_lead(i) == 1 % then target LEADS on this trial
                lead_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            end

        end

        % specify false alarm windows
        for i = 1:length(this_trial_masker_bash_times) % for each of the current masker times...
            [~,start_index_FA_window] = min(abs(tVec - (this_trial_masker_bash_times(i)+threshold_window_start))); % ...the false alarm window will start threshold_window_start seconds after the word onset
            [~,end_index_FA_window] = min(abs(tVec - (this_trial_masker_bash_times(i)+threshold_window_end))); % ...the false alarm window will end threshold_window_end seconds after the word onset

            if any(hit_windows(start_index_FA_window:end_index_FA_window) == 1)
                continue
            else
                FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                if masker_bash_lead(i) == 0 % then target LAGs on this trial
                    lag_FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                elseif masker_bash_lead(i) == 1 % then target LEADS on this trial
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
                %all_hit_rts
            elseif FA_windows(current_click_index) == 1 %...otherwise if that click falls within a false alarm window...
                all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                % populate lead and lag
                if lead_FA_windows(current_click_index) == 1
                    all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_FA_windows(current_click_index) == 1
                    all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end
                %all_FA_rts
            else% ...if the click is not counted as either
                clicks_not_counted = clicks_not_counted + 1;
            end
            total_clicks = total_clicks + 1;
        end

        % associate it with the correct condition

        
    end

disp(append(string(subject_ID(isubject,:)),': ', num2str((clicks_not_counted/total_clicks)*100), '% of clicks not counted'))
end

%% COLLAPSE OVER CONDITION
all_hits_collapsed_left_and_right = [];
all_hits_collapsed_left_and_right(1,:) = sum(all_hits(:,[1,2]),2); % itd50 / ild10
all_hits_collapsed_left_and_right(2,:) = sum(all_hits(:,[3,4]),2); % itd100 / ild10mag
all_hits_collapsed_left_and_right(3,:) = sum(all_hits(:,[5,6]),2); % itd200 /ild60
all_hits_collapsed_left_and_right(4,:) = sum(all_hits(:,[7,8]),2); % itd400 / ild60mag
all_hits_collapsed_left_and_right(5,:) = sum(all_hits(:,[9,10]),2); % itd400 / ild60mag
all_hits_collapsed_left_and_right(6,:) = sum(all_hits(:,[11,12]),2); % itd50 / ild10
all_hits_collapsed_left_and_right(7,:) = sum(all_hits(:,[13,14]),2); % itd100 / ild10mag
all_hits_collapsed_left_and_right(8,:) = sum(all_hits(:,[15,16]),2); % itd200 /ild60

all_lead_hits_collapsed_left_and_right = [];
all_lead_hits_collapsed_left_and_right(1,:) = sum(all_lead_hits(:,[1,2]),2); % itd50 / ild10
all_lead_hits_collapsed_left_and_right(2,:) = sum(all_lead_hits(:,[3,4]),2); % itd100 / ild10mag
all_lead_hits_collapsed_left_and_right(3,:) = sum(all_lead_hits(:,[5,6]),2); % itd200 /ild60
all_lead_hits_collapsed_left_and_right(4,:) = sum(all_lead_hits(:,[7,8]),2); % itd400 / ild60mag
all_lead_hits_collapsed_left_and_right(5,:) = sum(all_lead_hits(:,[9,10]),2); % itd400 / ild60mag
all_lead_hits_collapsed_left_and_right(6,:) = sum(all_lead_hits(:,[11,12]),2); % itd50 / ild10
all_lead_hits_collapsed_left_and_right(7,:) = sum(all_lead_hits(:,[13,14]),2); % itd100 / ild10mag
all_lead_hits_collapsed_left_and_right(8,:) = sum(all_lead_hits(:,[15,16]),2); % itd200 /ild60

all_lag_hits_collapsed_left_and_right = [];
all_lag_hits_collapsed_left_and_right(1,:) = sum(all_lag_hits(:,[1,2]),2); % itd50 / ild10
all_lag_hits_collapsed_left_and_right(2,:) = sum(all_lag_hits(:,[3,4]),2); % itd100 / ild10mag
all_lag_hits_collapsed_left_and_right(3,:) = sum(all_lag_hits(:,[5,6]),2); % itd200 /ild60
all_lag_hits_collapsed_left_and_right(4,:) = sum(all_lag_hits(:,[7,8]),2); % itd400 / ild60mag
all_lag_hits_collapsed_left_and_right(5,:) = sum(all_lag_hits(:,[9,10]),2); % itd400 / ild60mag
all_lag_hits_collapsed_left_and_right(6,:) = sum(all_lag_hits(:,[11,12]),2); % itd50 / ild10
all_lag_hits_collapsed_left_and_right(7,:) = sum(all_lag_hits(:,[13,14]),2); % itd100 / ild10mag
all_lag_hits_collapsed_left_and_right(8,:) = sum(all_lag_hits(:,[15,16]),2); % itd200 /ild60

all_FAs_collapsed_left_and_right = [];
all_FAs_collapsed_left_and_right(1,:) = sum(all_FAs(:,[1,2]),2); % itd50 / ild10
all_FAs_collapsed_left_and_right(2,:) = sum(all_FAs(:,[3,4]),2); % itd100 / ild10mag
all_FAs_collapsed_left_and_right(3,:) = sum(all_FAs(:,[5,6]),2); % itd200 /ild60
all_FAs_collapsed_left_and_right(4,:) = sum(all_FAs(:,[7,8]),2); % itd400 / ild60mag
all_FAs_collapsed_left_and_right(5,:) = sum(all_FAs(:,[9,10]),2); % itd400 / ild60mag
all_FAs_collapsed_left_and_right(6,:) = sum(all_FAs(:,[11,12]),2); % itd50 / ild10
all_FAs_collapsed_left_and_right(7,:) = sum(all_FAs(:,[13,14]),2); % itd100 / ild10mag
all_FAs_collapsed_left_and_right(8,:) = sum(all_FAs(:,[15,16]),2); % itd200 /ild60

all_lead_FAs_collapsed_left_and_right = [];
all_lead_FAs_collapsed_left_and_right(1,:) = sum(all_lead_FAs(:,[1,2]),2); % itd50 / ild10
all_lead_FAs_collapsed_left_and_right(2,:) = sum(all_lead_FAs(:,[3,4]),2); % itd100 / ild10mag
all_lead_FAs_collapsed_left_and_right(3,:) = sum(all_lead_FAs(:,[5,6]),2); % itd200 /ild60
all_lead_FAs_collapsed_left_and_right(4,:) = sum(all_lead_FAs(:,[7,8]),2); % itd400 / ild60mag
all_lead_FAs_collapsed_left_and_right(5,:) = sum(all_lead_FAs(:,[9,10]),2); % itd400 / ild60mag
all_lead_FAs_collapsed_left_and_right(6,:) = sum(all_lead_FAs(:,[11,12]),2); % itd50 / ild10
all_lead_FAs_collapsed_left_and_right(7,:) = sum(all_lead_FAs(:,[13,14]),2); % itd100 / ild10mag
all_lead_FAs_collapsed_left_and_right(8,:) = sum(all_lead_FAs(:,[15,16]),2); % itd200 /ild60
all_lag_FAs_collapsed_left_and_right = [];

all_lag_FAs_collapsed_left_and_right(1,:) = sum(all_lag_FAs(:,[1,2]),2); % itd50 / ild10
all_lag_FAs_collapsed_left_and_right(2,:) = sum(all_lag_FAs(:,[3,4]),2); % itd100 / ild10mag
all_lag_FAs_collapsed_left_and_right(3,:) = sum(all_lag_FAs(:,[5,6]),2); % itd200 /ild60
all_lag_FAs_collapsed_left_and_right(4,:) = sum(all_lag_FAs(:,[7,8]),2); % itd400 / ild60mag
all_lag_FAs_collapsed_left_and_right(5,:) = sum(all_lag_FAs(:,[9,10]),2); % itd400 / ild60mag
all_lag_FAs_collapsed_left_and_right(6,:) = sum(all_lag_FAs(:,[11,12]),2); % itd50 / ild10
all_lag_FAs_collapsed_left_and_right(7,:) = sum(all_lag_FAs(:,[13,14]),2); % itd100 / ild10mag
all_lag_FAs_collapsed_left_and_right(8,:) = sum(all_lag_FAs(:,[15,16]),2); % itd200 /ild60

all_num_target_bash_collapsed_left_and_right = [];
all_num_target_bash_collapsed_left_and_right(1,:) = sum(all_num_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_target_bash_collapsed_left_and_right(2,:) = sum(all_num_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_target_bash_collapsed_left_and_right(3,:) = sum(all_num_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_target_bash_collapsed_left_and_right(4,:) = sum(all_num_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_target_bash_collapsed_left_and_right(5,:) = sum(all_num_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_target_bash_collapsed_left_and_right(6,:) = sum(all_num_target_bash(:,[11,12]),2); % itd200 /ild60
all_num_target_bash_collapsed_left_and_right(7,:) = sum(all_num_target_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_target_bash_collapsed_left_and_right(8,:) = sum(all_num_target_bash(:,[15,16]),2); % itd400 / ild60mag

all_num_lead_target_bash_collapsed_left_and_right = [];
all_num_lead_target_bash_collapsed_left_and_right(1,:) = sum(all_num_lead_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_lead_target_bash_collapsed_left_and_right(2,:) = sum(all_num_lead_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lead_target_bash_collapsed_left_and_right(3,:) = sum(all_num_lead_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_lead_target_bash_collapsed_left_and_right(4,:) = sum(all_num_lead_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lead_target_bash_collapsed_left_and_right(5,:) = sum(all_num_lead_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lead_target_bash_collapsed_left_and_right(6,:) = sum(all_num_lead_target_bash(:,[11,12]),2); % itd200 /ild60
all_num_lead_target_bash_collapsed_left_and_right(7,:) = sum(all_num_lead_target_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_lead_target_bash_collapsed_left_and_right(8,:) = sum(all_num_lead_target_bash(:,[15,16]),2); % itd400 / ild60mag

all_num_lag_target_bash_collapsed_left_and_right = [];
all_num_lag_target_bash_collapsed_left_and_right(1,:) = sum(all_num_lag_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_lag_target_bash_collapsed_left_and_right(2,:) = sum(all_num_lag_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lag_target_bash_collapsed_left_and_right(3,:) = sum(all_num_lag_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_lag_target_bash_collapsed_left_and_right(4,:) = sum(all_num_lag_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lag_target_bash_collapsed_left_and_right(5,:) = sum(all_num_lag_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lag_target_bash_collapsed_left_and_right(6,:) = sum(all_num_lag_target_bash(:,[11,12]),2); % itd200 /ild60
all_num_lag_target_bash_collapsed_left_and_right(7,:) = sum(all_num_lag_target_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_lag_target_bash_collapsed_left_and_right(8,:) = sum(all_num_lag_target_bash(:,[15,16]),2); % itd400 / ild60mag

all_num_masker_bash_collapsed_left_and_right = [];
all_num_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_masker_bash(:,[11,12]),2); % itd200 /ild60
all_num_masker_bash_collapsed_left_and_right(7,:) = sum(all_num_masker_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_masker_bash_collapsed_left_and_right(8,:) = sum(all_num_masker_bash(:,[15,16]),2); % itd400 / ild60mag

all_num_lead_masker_bash_collapsed_left_and_right = [];
all_num_lead_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_lead_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_lead_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_lead_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lead_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_lead_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_lead_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_lead_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lead_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_lead_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lead_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_lead_masker_bash(:,[11,12]),2); % itd200 /ild60
all_num_lead_masker_bash_collapsed_left_and_right(7,:) = sum(all_num_lead_masker_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_lead_masker_bash_collapsed_left_and_right(8,:) = sum(all_num_lead_masker_bash(:,[15,16]),2); % itd400 / ild60mag

all_num_lag_masker_bash_collapsed_left_and_right = [];
all_num_lag_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_lag_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_lag_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_lag_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lag_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_lag_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_lag_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_lag_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lag_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_lag_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lag_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_lag_masker_bash(:,[11,12]),2); % itd200 /ild60
all_num_lag_masker_bash_collapsed_left_and_right(7,:) = sum(all_num_lag_masker_bash(:,[13,14]),2); % itd400 / ild60mag
all_num_lag_masker_bash_collapsed_left_and_right(8,:) = sum(all_num_lag_masker_bash(:,[15,16]),2); % itd400 / ild60mag

%% Calculate Rates
% Hit Rates
all_hit_rates = all_hits./all_num_target_bash;
all_lead_hit_rates = all_lead_hits./all_num_lead_target_bash;
all_lag_hit_rates = all_lag_hits./all_num_lag_target_bash;

all_hit_rates_collapsed = all_hits_collapsed_left_and_right./all_num_target_bash_collapsed_left_and_right;
all_lead_hit_rates_collapsed = all_lead_hits_collapsed_left_and_right./all_num_lead_target_bash_collapsed_left_and_right;
all_lag_hit_rates_collapsed = all_lag_hits_collapsed_left_and_right./all_num_lag_target_bash_collapsed_left_and_right;

all_hit_rates(all_hit_rates == 0) = 0.001;
all_hit_rates(all_hit_rates >= 1) = 0.999;
all_lead_hit_rates(all_lead_hit_rates == 0) = 0.001;
all_lead_hit_rates(all_lead_hit_rates >= 1) = 0.999;
all_lag_hit_rates(all_lag_hit_rates == 0) = 0.001;
all_lag_hit_rates(all_lag_hit_rates >= 1) = 0.999;

all_hit_rates_collapsed(all_hit_rates_collapsed == 0) = 0.001;
all_hit_rates_collapsed(all_hit_rates_collapsed >= 1) = 0.999;
all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed == 0) = 0.001;
all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed >= 1) = 0.999;
all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed == 0) = 0.001;
all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed >= 1) = 0.999;

% FA Rates
all_FA_rates = all_FAs./all_num_masker_bash;
all_lead_FA_rates = all_lead_FAs./all_num_lead_masker_bash;
all_lag_FA_rates = all_lag_FAs./all_num_lag_masker_bash;

all_FA_rates_collapsed = all_FAs_collapsed_left_and_right./all_num_masker_bash_collapsed_left_and_right;
all_lead_FA_rates_collapsed = all_lead_FAs_collapsed_left_and_right./all_num_lead_masker_bash_collapsed_left_and_right;
all_lag_FA_rates_collapsed = all_lag_FAs_collapsed_left_and_right./all_num_lag_masker_bash_collapsed_left_and_right;

all_FA_rates(all_FA_rates == 0) = 0.001;
all_FA_rates(all_FA_rates >= 1) = 0.999;
all_lead_FA_rates(all_lead_FA_rates == 0) = 0.001;
all_lead_FA_rates(all_lead_FA_rates >= 1) = 0.999;
all_lag_FA_rates(all_lag_FA_rates == 0) = 0.001;
all_lag_FA_rates(all_lag_FA_rates >= 1) = 0.999;

all_FA_rates_collapsed(all_FA_rates_collapsed == 0) = 0.001;
all_FA_rates_collapsed(all_FA_rates_collapsed >= 1) = 0.999;
all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed == 0) = 0.001;
all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed >= 1) = 0.999;
all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed == 0) = 0.001;
all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed >= 1) = 0.999;


%% D-prime calculation
d_primes_all = norminv(all_hit_rates) - norminv(all_FA_rates);
d_primes_collapsed = [];
d_primes_collapsed(1,:) = mean(d_primes_all(:,[1,2]),2); % itd50 / ild10
d_primes_collapsed(2,:) = mean(d_primes_all(:,[3,4]),2); % itd100 / ild10mag
d_primes_collapsed(3,:) = mean(d_primes_all(:,[5,6]),2); % itd200 /ild60
d_primes_collapsed(4,:) = mean(d_primes_all(:,[7,8]),2); % itd400 / ild60mag
d_primes_collapsed(5,:) = mean(d_primes_all(:,[9,10]),2); % itd50 / ild10
d_primes_collapsed(6,:) = mean(d_primes_all(:,[11,12]),2); % itd100 / ild10mag
d_primes_collapsed(7,:) = mean(d_primes_all(:,[13,14]),2); % itd200 /ild60
d_primes_collapsed(8,:) = mean(d_primes_all(:,[15,16]),2); % itd400 / ild60mag

lead_d_primes_all = norminv(all_lead_hit_rates) - norminv(all_lead_FA_rates);
lead_d_primes_collapsed = [];
lead_d_primes_collapsed(1,:) = mean(lead_d_primes_all(:,[1,2]),2); % itlead_d50 / illead_d10
lead_d_primes_collapsed(2,:) = mean(lead_d_primes_all(:,[3,4]),2); % itlead_d100 / illead_d10mag
lead_d_primes_collapsed(3,:) = mean(lead_d_primes_all(:,[5,6]),2); % itlead_d200 /illead_d60
lead_d_primes_collapsed(4,:) = mean(lead_d_primes_all(:,[7,8]),2); % itlead_d400 / illead_d60mag
lead_d_primes_collapsed(5,:) = mean(lead_d_primes_all(:,[9,10]),2); % itlead_d50 / illead_d10
lead_d_primes_collapsed(6,:) = mean(lead_d_primes_all(:,[11,12]),2); % itlead_d100 / illead_d10mag
lead_d_primes_collapsed(7,:) = mean(lead_d_primes_all(:,[13,14]),2); % itlead_d200 /illead_d60
lead_d_primes_collapsed(8,:) = mean(lead_d_primes_all(:,[15,16]),2); % itlead_d400 / illead_d60mag

lag_d_primes_all = norminv(all_lag_hit_rates) - norminv(all_lag_FA_rates);
lag_d_primes_collapsed = [];
lag_d_primes_collapsed(1,:) = mean(lag_d_primes_all(:,[1,2]),2); % itlag_d50 / illag_d10
lag_d_primes_collapsed(2,:) = mean(lag_d_primes_all(:,[3,4]),2); % itlag_d100 / illag_d10mag
lag_d_primes_collapsed(3,:) = mean(lag_d_primes_all(:,[5,6]),2); % itlag_d200 /illag_d60
lag_d_primes_collapsed(4,:) = mean(lag_d_primes_all(:,[7,8]),2); % itlag_d400 / illag_d60mag
lag_d_primes_collapsed(5,:) = mean(lag_d_primes_all(:,[9,10]),2); % itlag_d50 / illag_d10
lag_d_primes_collapsed(6,:) = mean(lag_d_primes_all(:,[11,12]),2); % itlag_d100 / illag_d10mag
lag_d_primes_collapsed(7,:) = mean(lag_d_primes_all(:,[13,14]),2); % itlag_d200 /illag_d60
lag_d_primes_collapsed(8,:) = mean(lag_d_primes_all(:,[15,16]),2); % itlag_d400 / illag_d60mag

%% PLOT
% Overall d', hit rate, FA rate
figure;
subplot(1,3,1)
hold on
plot(1:8,d_primes_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(d_primes_collapsed(1:2:8,:),2),std(d_primes_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(d_primes_collapsed(2:2:8,:),2),std(d_primes_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylabel("d'")
ylim([0 7])
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'Natural','Magnified'})

subplot(1,3,2)
hold on
plot(1:8,all_hit_rates_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(all_hit_rates_collapsed(1:2:8,:),2),std(all_hit_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(all_hit_rates_collapsed(2:2:8,:),2),std(all_hit_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
ylabel('Hit Rate')
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'','Natural','Magnified'})

subplot(1,3,3)
hold on
plot(1:8,all_FA_rates_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(all_FA_rates_collapsed(1:2:8,:),2),std(all_FA_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(all_FA_rates_collapsed(2:2:8,:),2),std(all_FA_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
ylabel('FA Rate')
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'Natural','Magnified'})


% Lead vs. lag d-prime
figure;
subplot(2,3,1)
hold on
plot(1:8,lead_d_primes_collapsed,'LineStyle','-','Color','k')
p1 = errorbar(1:2:8,mean(lead_d_primes_collapsed(1:2:8,:),2),std(lead_d_primes_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
p2 = errorbar(2:2:8,mean(lead_d_primes_collapsed(2:2:8,:),2),std(lead_d_primes_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
xticks(1:8)
ylim([0 7])
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lead','fontsize',18)
%suptitle("d'")

subplot(2,3,4)
hold on
plot(1:8,lag_d_primes_collapsed,'LineStyle','-','Color','k')
p1 = errorbar(1:2:8,mean(lag_d_primes_collapsed(1:2:8,:),2),std(lag_d_primes_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r');
p2 = errorbar(2:2:8,mean(lag_d_primes_collapsed(2:2:8,:),2),std(lag_d_primes_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b');
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lag','fontsize',18)
ylim([0,7])

% lead vs. lag hit rate
subplot(2,3,2)
hold on
plot(1:8,all_lead_hit_rates_collapsed,'LineStyle','-','Color','k')
p1 = errorbar(1:2:8,mean(all_lead_hit_rates_collapsed(1:2:8,:),2),std(all_lead_hit_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
p2 = errorbar(2:2:8,mean(all_lead_hit_rates_collapsed(2:2:8,:),2),std(all_lead_hit_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
title('Hit Rate')
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lead','fontsize',18)

subplot(2,3,5)
hold on
plot(1:8,all_lag_hit_rates_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(all_lag_hit_rates_collapsed(1:2:8,:),2),std(all_lag_hit_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(all_lag_hit_rates_collapsed(2:2:8,:),2),std(all_lag_hit_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
title('Hit Rate')
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'','Natural','Magnified'})
title('lag','fontsize',18)
%suptitle("Hit Rates")

% lead vs. lag FA rate
subplot(2,3,3)
hold on
plot(1:8,all_lead_FA_rates_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(all_lead_FA_rates_collapsed(1:2:8,:),2),std(all_lead_FA_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(all_lead_FA_rates_collapsed(2:2:8,:),2),std(all_lead_FA_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
title('FA Rate')
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'','Natural','Magnified'})
title('lead','fontsize',18)

subplot(2,3,6)
hold on
plot(1:8,all_lag_FA_rates_collapsed,'LineStyle','-','Color','k')
errorbar(1:2:8,mean(all_lag_FA_rates_collapsed(1:2:8,:),2),std(all_lag_FA_rates_collapsed(1:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r')
errorbar(2:2:8,mean(all_lag_FA_rates_collapsed(2:2:8,:),2),std(all_lag_FA_rates_collapsed(2:2:8,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b')
ylim([0 1])
xticks(1:8)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag','30deg','30degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend({'','Natural','Magnified'})
title('lag','fontsize',18)
%suptitle("FA Rates")

%% Save data
save('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Behavior_Results.mat','d_primes_collapsed','d_primes_collapsed','all_hit_rates_collapsed','all_FA_rates_collapsed')

hit_rate_table = array2table(all_hit_rates_collapsed);
writetable(rows2vars(hit_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Hit_Rates.csv')

FA_rate_table = array2table(all_FA_rates_collapsed(5:end,:));
writetable(rows2vars(FA_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_FA_Rates.csv')

d_prime_table = array2table(d_primes_collapsed(5:end,:));
writetable(rows2vars(d_prime_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_d_primes.csv')