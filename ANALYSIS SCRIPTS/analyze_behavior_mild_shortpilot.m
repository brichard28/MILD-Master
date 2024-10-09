% Benjamin Richardson
% Created: September 16th, 2024

% Script to analyze behavioral sensitivity (d-prime) for MILD MASTER

BehaviorTable = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master.xlsx','Format','auto');

subject_ID = char('shorttest2','shorttest3','shorttest5','shorttest6'); % ild pilot % ,'shorttest4'
num_conditions = 12;
min_click_times = [];
all_click_times = [];

all_hits = zeros(size(subject_ID,1),num_conditions);
all_lead_hits = zeros(size(subject_ID,1),num_conditions);
all_lag_hits = zeros(size(subject_ID,1),num_conditions);
all_hit_rts = zeros(size(subject_ID,1),num_conditions);
all_FAs = zeros(size(subject_ID,1),num_conditions);
all_lead_FAs = zeros(size(subject_ID,1),num_conditions);
all_lag_FAs = zeros(size(subject_ID,1),num_conditions);
all_FA_rts = zeros(size(subject_ID,1),num_conditions);

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
'side=l_itd=0_az=20_mag=1'};

fs = 44100;
cue_dur = 1.3;

for isubject = 1:size(subject_ID,1) % For each subject...
    disp(string(subject_ID(isubject,:)))
    clicks_not_counted = 0;
    total_clicks = 0;

    % Load the word times for this subject
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\mild-master__s_" + strtrim(string(subject_ID(isubject,:))) + "__Word_Times.csv");

    run_count_per_condition = -1*ones(1,num_conditions); % array to keep track of which run in each condition we are on

    % Find the rows associated with this subject
    rows_this_subject = find(BehaviorTable.S == strtrim(string(subject_ID(isubject,:))));

    all_hit_rts_this_subject = struct('rts',[]);
    all_hit_rts_this_subject = repmat(all_hit_rts_this_subject,num_conditions,1);
    all_FA_rts_this_subject = struct('rts',[]);
    all_FA_rts_this_subject = repmat(all_FA_rts_this_subject,num_conditions,1);
    distances_to_nearest_bash = [];

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

       threshold_window_start = 0.1; %0.2
       threshold_window_end =  0.8; % 1.15
       tVec = 0:1/44100:10;
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
        min_click_times = [min_click_times, min(this_trial_click_times)];
        all_click_times = [all_click_times, this_trial_click_times];
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
                distances_to_bash = [];
                for ibashtime = 1:length(this_trial_target_bash_times)
                    distances_to_bash = [distances_to_bash, this_trial_click_times(iclick) - this_trial_target_bash_times(ibashtime)];
                end
                distances_to_bash(distances_to_bash < 0) =[];
                all_hit_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts = [all_hit_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts, min(distances_to_bash)];
            
            elseif FA_windows(current_click_index) == 1 %...otherwise if that click falls within a false alarm window...
                all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                % populate lead and lag
                if lead_FA_windows(current_click_index) == 1
                    all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_FA_windows(current_click_index) == 1
                    all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end
                

                % this reaction time
                distances_to_bash = [];
                for ibashtime = 1:length(this_trial_masker_bash_times)
                    distances_to_bash = [distances_to_bash, this_trial_click_times(iclick) - this_trial_masker_bash_times(ibashtime)];
                end
                distances_to_bash(distances_to_bash < 0) =[];
                all_FA_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts = [all_FA_rts_this_subject(string(all_maskers) == string(this_trial_masker)).rts, min(distances_to_bash)];
            
            else% ...if the click is not counted as either
                clicks_not_counted = clicks_not_counted + 1;
            end
            total_clicks = total_clicks + 1;
        end

        % Calculate distances to nearest click
        all_bash_click_distances= [];
        for ibashtime = 1:length(this_trial_target_bash_times)
            all_bash_click_distances(ibashtime,:) = this_trial_click_times - this_trial_target_bash_times(ibashtime);
        end
        all_bash_click_distances(all_bash_click_distances < 0) = nan;

        %% Find the nearest color time to each click (minimum positive value of click_distances in each column)
        [~,nearest_click] = min(abs(all_bash_click_distances),[],1); % find the nearest click to each target word
        for i = 1:length(this_trial_click_times)
            if isnan(all_bash_click_distances(:,i)) == ones(1,length(this_trial_target_bash_times)) % all of these clicks were before the first word
                nearest_click(i) = nan;
            else
                distances_to_nearest_bash = [distances_to_nearest_bash, all_bash_click_distances(nearest_click(i),i)];
            end

        end

        
    end
    all_distances_to_nearest_bash(isubject).distances = distances_to_nearest_bash;

    for icondition = 1:num_conditions
    all_hit_rts(isubject,icondition) = nanmean([all_hit_rts_this_subject(icondition).rts]);
    all_FA_rts(isubject,icondition) = nanmean([all_FA_rts_this_subject(icondition).rts]);
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

all_lead_hits_collapsed_left_and_right = [];
all_lead_hits_collapsed_left_and_right(1,:) = sum(all_lead_hits(:,[1,2]),2); % itd50 / ild10
all_lead_hits_collapsed_left_and_right(2,:) = sum(all_lead_hits(:,[3,4]),2); % itd100 / ild10mag
all_lead_hits_collapsed_left_and_right(3,:) = sum(all_lead_hits(:,[5,6]),2); % itd200 /ild60
all_lead_hits_collapsed_left_and_right(4,:) = sum(all_lead_hits(:,[7,8]),2); % itd400 / ild60mag
all_lead_hits_collapsed_left_and_right(5,:) = sum(all_lead_hits(:,[9,10]),2); % itd400 / ild60mag
all_lead_hits_collapsed_left_and_right(6,:) = sum(all_lead_hits(:,[11,12]),2); % itd50 / ild10

all_lag_hits_collapsed_left_and_right = [];
all_lag_hits_collapsed_left_and_right(1,:) = sum(all_lag_hits(:,[1,2]),2); % itd50 / ild10
all_lag_hits_collapsed_left_and_right(2,:) = sum(all_lag_hits(:,[3,4]),2); % itd100 / ild10mag
all_lag_hits_collapsed_left_and_right(3,:) = sum(all_lag_hits(:,[5,6]),2); % itd200 /ild60
all_lag_hits_collapsed_left_and_right(4,:) = sum(all_lag_hits(:,[7,8]),2); % itd400 / ild60mag
all_lag_hits_collapsed_left_and_right(5,:) = sum(all_lag_hits(:,[9,10]),2); % itd400 / ild60mag
all_lag_hits_collapsed_left_and_right(6,:) = sum(all_lag_hits(:,[11,12]),2); % itd50 / ild10

all_FAs_collapsed_left_and_right = [];
all_FAs_collapsed_left_and_right(1,:) = sum(all_FAs(:,[1,2]),2); % itd50 / ild10
all_FAs_collapsed_left_and_right(2,:) = sum(all_FAs(:,[3,4]),2); % itd100 / ild10mag
all_FAs_collapsed_left_and_right(3,:) = sum(all_FAs(:,[5,6]),2); % itd200 /ild60
all_FAs_collapsed_left_and_right(4,:) = sum(all_FAs(:,[7,8]),2); % itd400 / ild60mag
all_FAs_collapsed_left_and_right(5,:) = sum(all_FAs(:,[9,10]),2); % itd400 / ild60mag
all_FAs_collapsed_left_and_right(6,:) = sum(all_FAs(:,[11,12]),2); % itd50 / ild10

all_lead_FAs_collapsed_left_and_right = [];
all_lead_FAs_collapsed_left_and_right(1,:) = sum(all_lead_FAs(:,[1,2]),2); % itd50 / ild10
all_lead_FAs_collapsed_left_and_right(2,:) = sum(all_lead_FAs(:,[3,4]),2); % itd100 / ild10mag
all_lead_FAs_collapsed_left_and_right(3,:) = sum(all_lead_FAs(:,[5,6]),2); % itd200 /ild60
all_lead_FAs_collapsed_left_and_right(4,:) = sum(all_lead_FAs(:,[7,8]),2); % itd400 / ild60mag
all_lead_FAs_collapsed_left_and_right(5,:) = sum(all_lead_FAs(:,[9,10]),2); % itd400 / ild60mag
all_lead_FAs_collapsed_left_and_right(6,:) = sum(all_lead_FAs(:,[11,12]),2); % itd50 / ild10
all_lag_FAs_collapsed_left_and_right = [];

all_lag_FAs_collapsed_left_and_right(1,:) = sum(all_lag_FAs(:,[1,2]),2); % itd50 / ild10
all_lag_FAs_collapsed_left_and_right(2,:) = sum(all_lag_FAs(:,[3,4]),2); % itd100 / ild10mag
all_lag_FAs_collapsed_left_and_right(3,:) = sum(all_lag_FAs(:,[5,6]),2); % itd200 /ild60
all_lag_FAs_collapsed_left_and_right(4,:) = sum(all_lag_FAs(:,[7,8]),2); % itd400 / ild60mag
all_lag_FAs_collapsed_left_and_right(5,:) = sum(all_lag_FAs(:,[9,10]),2); % itd400 / ild60mag
all_lag_FAs_collapsed_left_and_right(6,:) = sum(all_lag_FAs(:,[11,12]),2); % itd50 / ild10

all_num_target_bash_collapsed_left_and_right = [];
all_num_target_bash_collapsed_left_and_right(1,:) = sum(all_num_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_target_bash_collapsed_left_and_right(2,:) = sum(all_num_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_target_bash_collapsed_left_and_right(3,:) = sum(all_num_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_target_bash_collapsed_left_and_right(4,:) = sum(all_num_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_target_bash_collapsed_left_and_right(5,:) = sum(all_num_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_target_bash_collapsed_left_and_right(6,:) = sum(all_num_target_bash(:,[11,12]),2); % itd200 /ild60

all_num_lead_target_bash_collapsed_left_and_right = [];
all_num_lead_target_bash_collapsed_left_and_right(1,:) = sum(all_num_lead_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_lead_target_bash_collapsed_left_and_right(2,:) = sum(all_num_lead_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lead_target_bash_collapsed_left_and_right(3,:) = sum(all_num_lead_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_lead_target_bash_collapsed_left_and_right(4,:) = sum(all_num_lead_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lead_target_bash_collapsed_left_and_right(5,:) = sum(all_num_lead_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lead_target_bash_collapsed_left_and_right(6,:) = sum(all_num_lead_target_bash(:,[11,12]),2); % itd200 /ild60

all_num_lag_target_bash_collapsed_left_and_right = [];
all_num_lag_target_bash_collapsed_left_and_right(1,:) = sum(all_num_lag_target_bash(:,[1,2]),2); % itd50 / ild10
all_num_lag_target_bash_collapsed_left_and_right(2,:) = sum(all_num_lag_target_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lag_target_bash_collapsed_left_and_right(3,:) = sum(all_num_lag_target_bash(:,[5,6]),2); % itd200 /ild60
all_num_lag_target_bash_collapsed_left_and_right(4,:) = sum(all_num_lag_target_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lag_target_bash_collapsed_left_and_right(5,:) = sum(all_num_lag_target_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lag_target_bash_collapsed_left_and_right(6,:) = sum(all_num_lag_target_bash(:,[11,12]),2); % itd200 /ild60

all_num_masker_bash_collapsed_left_and_right = [];
all_num_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_masker_bash(:,[11,12]),2); % itd200 /ild60

all_num_lead_masker_bash_collapsed_left_and_right = [];
all_num_lead_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_lead_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_lead_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_lead_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lead_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_lead_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_lead_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_lead_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lead_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_lead_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lead_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_lead_masker_bash(:,[11,12]),2); % itd200 /ild60

all_num_lag_masker_bash_collapsed_left_and_right = [];
all_num_lag_masker_bash_collapsed_left_and_right(1,:) = sum(all_num_lag_masker_bash(:,[1,2]),2); % itd50 / ild10
all_num_lag_masker_bash_collapsed_left_and_right(2,:) = sum(all_num_lag_masker_bash(:,[3,4]),2); % itd100 / ild10mag
all_num_lag_masker_bash_collapsed_left_and_right(3,:) = sum(all_num_lag_masker_bash(:,[5,6]),2); % itd200 /ild60
all_num_lag_masker_bash_collapsed_left_and_right(4,:) = sum(all_num_lag_masker_bash(:,[7,8]),2); % itd400 / ild60mag
all_num_lag_masker_bash_collapsed_left_and_right(5,:) = sum(all_num_lag_masker_bash(:,[9,10]),2); % itd400 / ild60mag
all_num_lag_masker_bash_collapsed_left_and_right(6,:) = sum(all_num_lag_masker_bash(:,[11,12]),2); % itd200 /ild60

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

lead_d_primes_all = norminv(all_lead_hit_rates) - norminv(all_lead_FA_rates);
lead_d_primes_collapsed = [];
lead_d_primes_collapsed(1,:) = mean(lead_d_primes_all(:,[1,2]),2); % itlead_d50 / illead_d10
lead_d_primes_collapsed(2,:) = mean(lead_d_primes_all(:,[3,4]),2); % itlead_d100 / illead_d10mag
lead_d_primes_collapsed(3,:) = mean(lead_d_primes_all(:,[5,6]),2); % itlead_d200 /illead_d60
lead_d_primes_collapsed(4,:) = mean(lead_d_primes_all(:,[7,8]),2); % itlead_d400 / illead_d60mag
lead_d_primes_collapsed(5,:) = mean(lead_d_primes_all(:,[9,10]),2); % itlead_d50 / illead_d10
lead_d_primes_collapsed(6,:) = mean(lead_d_primes_all(:,[11,12]),2); % itlead_d100 / illead_d10mag

lag_d_primes_all = norminv(all_lag_hit_rates) - norminv(all_lag_FA_rates);
lag_d_primes_collapsed = [];
lag_d_primes_collapsed(1,:) = mean(lag_d_primes_all(:,[1,2]),2); % itlag_d50 / illag_d10
lag_d_primes_collapsed(2,:) = mean(lag_d_primes_all(:,[3,4]),2); % itlag_d100 / illag_d10mag
lag_d_primes_collapsed(3,:) = mean(lag_d_primes_all(:,[5,6]),2); % itlag_d200 /illag_d60
lag_d_primes_collapsed(4,:) = mean(lag_d_primes_all(:,[7,8]),2); % itlag_d400 / illag_d60mag
lag_d_primes_collapsed(5,:) = mean(lag_d_primes_all(:,[9,10]),2); % itlag_d50 / illag_d10
lag_d_primes_collapsed(6,:) = mean(lag_d_primes_all(:,[11,12]),2); % itlag_d100 / illag_d10mag

%% PLOT
% Overall d', hit rate, FA rate
figure;
subplot(1,3,1)
hold on
for i = [1,3,5]
plot(i:i+1,d_primes_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(d_primes_collapsed(1:2:6,:),2),std(d_primes_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(d_primes_collapsed(2:2:6,:),2),std(d_primes_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylabel("d'",'FontSize',18)
ylim([0 7])
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1), p2(1)],{'Natural','Magnified'})

subplot(1,3,2)
hold on
for i = [1,3,5]
plot(i:i+1,all_hit_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_hit_rates_collapsed(1:2:6,:),2),std(all_hit_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_hit_rates_collapsed(2:2:6,:),2),std(all_hit_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel('Hit Rate','FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})

subplot(1,3,3)
hold on
for i = [1,3,5]
plot(i:i+1,all_FA_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_FA_rates_collapsed(1:2:6,:),2),std(all_FA_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_FA_rates_collapsed(2:2:6,:),2),std(all_FA_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel('FA Rate','FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})


% Lead vs. lag d-prime
figure;
subplot(2,3,1)
hold on
for i = [1,3,5]
plot(i:i+1,lead_d_primes_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(lead_d_primes_collapsed(1:2:6,:),2),std(lead_d_primes_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(lead_d_primes_collapsed(2:2:6,:),2),std(lead_d_primes_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
xticks(1:6)
ylim([0 7])
ylabel("d'",'FontSize',18)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
%suptitle("d'")

subplot(2,3,4)
hold on
for i = [1,3,5]
plot(i:i+1,lag_d_primes_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(lag_d_primes_collapsed(1:2:6,:),2),std(lag_d_primes_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(lag_d_primes_collapsed(2:2:6,:),2),std(lag_d_primes_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
xticks(1:6)
ylabel("d'",'FontSize',18)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
ylim([0,7])

% lead vs. lag hit rate
subplot(2,3,2)
hold on
for i = [1,3,5]
plot(i:i+1,all_lead_hit_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_lead_hit_rates_collapsed(1:2:6,:),2),std(all_lead_hit_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_lead_hit_rates_collapsed(2:2:6,:),2),std(all_lead_hit_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel("Hit Rate",'FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lead','fontsize',18)

subplot(2,3,5)
hold on
for i = [1,3,5]
plot(i:i+1,all_lag_hit_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_lag_hit_rates_collapsed(1:2:6,:),2),std(all_lag_hit_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_lag_hit_rates_collapsed(2:2:6,:),2),std(all_lag_hit_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel("Hit Rate",'FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lag','fontsize',18)
%suptitle("Hit Rates")

% lead vs. lag FA rate
subplot(2,3,3)
hold on
for i = [1,3,5]
plot(i:i+1,all_lead_FA_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_lead_FA_rates_collapsed(1:2:6,:),2),std(all_lead_FA_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_lead_FA_rates_collapsed(2:2:6,:),2),std(all_lead_FA_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel("FA Rate",'FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lead','fontsize',18)

subplot(2,3,6)
hold on
for i = [1,3,5]
plot(i:i+1,all_lag_FA_rates_collapsed(i:i+1,:),'LineStyle','-','Color','k')
end
p1 = errorbar(1:2:6,mean(all_lag_FA_rates_collapsed(1:2:6,:),2),std(all_lag_FA_rates_collapsed(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,mean(all_lag_FA_rates_collapsed(2:2:6,:),2),std(all_lag_FA_rates_collapsed(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylim([0 1])
ylabel("FA Rate",'FontSize',18)
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)','FontSize',18)
legend([p1(1),p2(1)],{'Natural','Magnified'})
title('lag','fontsize',18)
%suptitle("FA Rates")

%% Click Times Histograms
figure;
for isubject = 1:size(subject_ID,1)
    hold on;
    histogram(all_distances_to_nearest_bash(isubject).distances,50)
end
xlabel('Time re: closest bash onset','FontSize',18)
ylabel('Number of Clicks','Fontsize',18)


%% Reaction times
all_hit_rts_collapsed_left_and_right = [];
all_hit_rts_collapsed_left_and_right(1,:) = nanmean(all_hit_rts(:,[1,2]),2); % itd50 / ild10
all_hit_rts_collapsed_left_and_right(2,:) = nanmean(all_hit_rts(:,[3,4]),2); % itd100 / ild10mag
all_hit_rts_collapsed_left_and_right(3,:) = nanmean(all_hit_rts(:,[5,6]),2); % itd200 /ild60
all_hit_rts_collapsed_left_and_right(4,:) = nanmean(all_hit_rts(:,[7,8]),2); % itd400 / ild60mag
all_hit_rts_collapsed_left_and_right(5,:) = nanmean(all_hit_rts(:,[9,10]),2); % itd400 / ild60mag
all_hit_rts_collapsed_left_and_right(6,:) = nanmean(all_hit_rts(:,[11,12]),2); % itd50 / ild10
all_FA_rts_collapsed_left_and_right = [];
all_FA_rts_collapsed_left_and_right(1,:) = nanmean(all_FA_rts(:,[1,2]),2); % itd50 / ild10
all_FA_rts_collapsed_left_and_right(2,:) = nanmean(all_FA_rts(:,[3,4]),2); % itd100 / ild10mag
all_FA_rts_collapsed_left_and_right(3,:) = nanmean(all_FA_rts(:,[5,6]),2); % itd200 /ild60
all_FA_rts_collapsed_left_and_right(4,:) = nanmean(all_FA_rts(:,[7,8]),2); % itd400 / ild60mag
all_FA_rts_collapsed_left_and_right(5,:) = nanmean(all_FA_rts(:,[9,10]),2); % itd400 / ild60mag
all_FA_rts_collapsed_left_and_right(6,:) = nanmean(all_FA_rts(:,[11,12]),2); % itd50 / ild10

% Hit RTs
figure;
hold on
plot(1:6,all_hit_rts_collapsed_left_and_right,'LineStyle','-','Color','k')
p1 = errorbar(1:2:6,nanmean(all_hit_rts_collapsed_left_and_right(1:2:6,:),2),std(all_hit_rts_collapsed_left_and_right(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,nanmean(all_hit_rts_collapsed_left_and_right(2:2:6,:),2),std(all_hit_rts_collapsed_left_and_right(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylabel("Hit Reaction Time (s)")
ylim([0 1.5])
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend([p1(1), p2(1)],{'Natural','Magnified'})

% FA RTs
figure;
hold on
plot(1:6,all_FA_rts_collapsed_left_and_right,'LineStyle','-','Color','k')
p1 = errorbar(1:2:6,nanmean(all_FA_rts_collapsed_left_and_right(1:2:6,:),2),std(all_FA_rts_collapsed_left_and_right(1:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','r','MarkerFaceColor','r');
p2 = errorbar(2:2:6,nanmean(all_FA_rts_collapsed_left_and_right(2:2:6,:),2),std(all_FA_rts_collapsed_left_and_right(2:2:6,:),[],2)/(sqrt(size(subject_ID,1))-1),'o','Color','b','MarkerFaceColor','b');
ylabel("FA Reaction Time (s)")
ylim([0 1.5])
xticks(1:6)
xticklabels({'5deg','5degMag','10deg','10degMag','20deg','20degMag'})
%xlabel('ITD (us)')
xlabel('Natural ILD (deg)')
legend([p1(1), p2(1)],{'Natural','Magnified'})
%% Save data
save('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Behavior_Results.mat','d_primes_collapsed','d_primes_collapsed','all_hit_rates_collapsed','all_FA_rates_collapsed')

hit_rate_table = array2table(all_hit_rates_collapsed);
writetable(rows2vars(hit_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_Hit_Rates.csv')

FA_rate_table = array2table(all_FA_rates_collapsed(5:end,:));
writetable(rows2vars(FA_rate_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_FA_Rates.csv')

d_prime_table = array2table(d_primes_collapsed(5:end,:));
writetable(rows2vars(d_prime_table),'C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER_d_primes.csv')