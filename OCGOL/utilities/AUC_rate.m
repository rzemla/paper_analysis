function [session_vars] = AUC_rate(session_vars,options)

%get AUC/min calculation for all sig calcium events in run epochs 

%% Define parameters

sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;
%frame period (s) - first session, A trial - should be the same for all
dt = session_vars{1}.Imaging_split{selectTrial(1)}.dt;

%% Get run time in each subset of trials

for ss=sessionSelect
    %get time (min) for - run
    A.run.time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(1)}.run_ones ==1))*dt)/60;
    B.run.time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(2)}.run_ones ==1))*dt)/60;
    %norun
    A.norun.time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(1)}.run_ones ==0))*dt)/60;
    B.norun.time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(2)}.run_ones ==0))*dt)/60;
end

%get AUC sum for each ROI
for ss=sessionSelect
    %for each ROI, get AUC sum - run
    AUC_sum{ss}.run.A = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC,'UniformOutput',true);
    AUC_sum{ss}.run.B = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC,'UniformOutput',true);
    %no run
    AUC_sum{ss}.norun.A = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.AUC,'UniformOutput',true);
    AUC_sum{ss}.norun.B = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.AUC,'UniformOutput',true);

    %for each ROI get AUC/min - run
    AUC_min{ss}.run.A = AUC_sum{ss}.run.A./A.run.time(ss);
    AUC_min{ss}.run.B = AUC_sum{ss}.run.B./B.run.time(ss);

    %no run
    AUC_min{ss}.norun.A = AUC_sum{ss}.norun.A./A.norun.time(ss);
    AUC_min{ss}.norun.B = AUC_sum{ss}.norun.B./B.norun.time(ss);
end

%get run event frequency (nb_events/min)
for ss=sessionSelect
    %run
    event_freq{ss}.run.A = session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.nb_events./A.run.time(ss);
    event_freq{ss}.run.B = session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.nb_events./B.run.time(ss);
    %no run
    event_freq{ss}.norun.A = session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.nb_events./A.norun.time(ss);
    event_freq{ss}.norun.B = session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.nb_events./B.norun.time(ss);    
end

%% Assign to session_vars.struct for export
for ss=sessionSelect
    %run event freq assignment
    session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.event_freq = event_freq{ss}.run.A;
    session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.event_freq = event_freq{ss}.run.B;
    %run AUC/min assignment
    session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min = AUC_min{ss}.run.A;
    session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min = AUC_min{ss}.run.B;

    %no run event freq assignment
    session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.event_freq = event_freq{ss}.norun.A;
    session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.event_freq = event_freq{ss}.norun.B;

    %no run AUC/min assignment
    session_vars{ss}.Events_split{selectTrial(1)}.NoRun.properties.AUC_min = AUC_min{ss}.norun.A;
    session_vars{ss}.Events_split{selectTrial(2)}.NoRun.properties.AUC_min = AUC_min{ss}.norun.B;
end

