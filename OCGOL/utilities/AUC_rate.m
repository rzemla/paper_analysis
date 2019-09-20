function [session_vars] = AUC_rate(session_vars,options)

%get AUC/min calculation for all sig calcium events in run epochs 

%% Define parameters

sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;
%frame period (s) - first session, A trial - should be the same for all
dt = session_vars{1}.Imaging_split{selectTrial(1)}.dt;

%% Get run time in each subset of trials

for ss=sessionSelect
    %get time (min) for
    A_run_time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(1)}.run_ones ==1))*dt)/60;
    B_run_time(ss) = (length(find(session_vars{ss}.Behavior_split{selectTrial(2)}.run_ones ==1))*dt)/60;
end

%get AUC sum for each ROI
for ss=sessionSelect
    %for each ROI, get AUC sum
    AUC_sum.A{ss} = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC,'UniformOutput',true);
    AUC_sum.B{ss} = cellfun(@sum,session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC,'UniformOutput',true);

    %for each ROI get AUC/min
    AUC_min.A{ss} = AUC_sum.A{ss}./A_run_time(ss);
    AUC_min.B{ss} = AUC_sum.B{ss}./B_run_time(ss);
end

%get run event frequency (nb_events/min)
for ss=sessionSelect
    event_freq.A{ss} = session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.nb_events./A_run_time(ss);
    event_freq.B{ss} = session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.nb_events./B_run_time(ss);
end

%% Assign to session_vars.struct for export
for ss=sessionSelect
    session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.event_freq = event_freq.A{ss};
    session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.event_freq = event_freq.B{ss};
    
    session_vars{ss}.Events_split{selectTrial(1)}.Run.properties.AUC_min = AUC_min.A{ss};
    session_vars{ss}.Events_split{selectTrial(2)}.Run.properties.AUC_min = AUC_min.B{ss};
end

