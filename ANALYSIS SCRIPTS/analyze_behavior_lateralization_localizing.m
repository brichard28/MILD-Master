%% analyze_behavior_lateralization_matching.m

BehaviorTable = readtable('C:\Users\benri\Documents\GitHub\MILD-Master\RESULTS DATA\MILD-MASTER Behavior Files\lateralization-localizing.xlsx','Format','auto');
%subject_ID = char('latpilot1','latpilot2','latmatch3','latpilot4','latpilot5','latpilot6'); %
subject_ID = char('lateralization_1_ITD');
for isubject = 1:size(subject_ID,1) % For each subject...
    rows_this_subject = find(BehaviorTable.S == strtrim(string(subject_ID(isubject,:))));

    num_trials_per_condition = 6;
    ref_itds = [-800,-500,-400,-200,-100,-50,50,100,200,400,500,800];
    ref_ilds = [-90,-60,-30,-15,-10,-5,5,10,15,30,60,90];
    num_conditions = length(ref_itds);
    run_count_per_condition = zeros(1,num_conditions);

    responses_this_subject = nan(length(ref_itds),num_trials_per_condition);
    for itrial = 1:length(rows_this_subject)

        this_trial_ref_type = BehaviorTable.cue_type(rows_this_subject(itrial));
        this_trial_ref_value = BehaviorTable.location(rows_this_subject(itrial)); % find the condition for this trial

        if this_trial_ref_type == "itd"

            run_count_per_condition(ref_itds == this_trial_ref_value) = run_count_per_condition(ref_itds == this_trial_ref_value) + 1;
    
            this_trial_response = BehaviorTable.Resp(rows_this_subject(itrial));
            responses_this_subject(find(ref_itds == this_trial_ref_value),run_count_per_condition(ref_itds == this_trial_ref_value)) = this_trial_response;
        elseif this_trial_ref_type == "ild"
            run_count_per_condition(ref_ilds == this_trial_ref_value) = run_count_per_condition(ref_ilds == this_trial_ref_value) + 1;
    
            this_trial_response = BehaviorTable.Resp(rows_this_subject(itrial));
            responses_this_subject(find(ref_ilds == this_trial_ref_value),run_count_per_condition(ref_ilds == this_trial_ref_value)) = this_trial_response;

        end
    end

    all_responses(isubject,:) = mean(responses_this_subject,2);
end

figure;
hold on
plot(ref_itds,all_responses);%,'-','Color','black')
%errorbar(ref_itds,mean(all_responses,1),std(all_responses,[],1)/(sqrt(size(all_responses,1)) - 1),'-o','Color','red','MarkerFaceColor','red')
xlabel('ITD (us)','FontSize',18)
xticks(ref_itds)
xticklabels(string(ref_itds))
ylim([-90 90])
ylabel('Matched Degree az','FontSize',18)
ax = gca;
ax.XTickLabelRotation = 60;
title('SSD Subject with ITDs','FontSize',18)