function [Behavior_split,Imaging_split,Events_split,Behavior_split_lap,Events_split_lap] = split_trials_OCGOL(Behavior, Imaging, Events)
% takes behavioral and calcium data as inputs from all laps and splits them
% into two structures in each behavioral, imaging, and calcium cell 
%INPUTS (structs):
%Behavior - cell containing the time, position, lick, cumulative position
%and lap restricted data
%Imaging - calcium trace and time data and restricted calcium trace and time data
%Events - onset of significant calcium events

%trialOrder = vector of trial types for each lap starting with first odor
%trial --> either 2 or 3
%2 = far reward, odor A (amyl acetate)
%3 = near reward, odor B (alpha-pinene)

%V3 - creates dF/F mask for calculating the calcium properties by trial
%type and generates other split inputs for event_properties script

%% Import behavioral variables

%behavior variables
%time resampled to matching imaging rate
time = Behavior.resampled.time;

% 1-~200 position for each lap
position = Behavior.resampled.position;

%normalized position
position_norm = Behavior.resampled.normalizedposition;

%position of animal on each frame
lapFrames = Behavior.resampled.lapNb;

%times when animal is running (all laps) 
run_ones = Behavior.run_ones;

%start and stop times (s) for each lap
Behavior.restricted.lap;

%Imaging variables
%Imaging time 
imaging_time = Imaging.time_restricted;
%ROI traced
imaging_trace = Imaging.trace_restricted;

%Events variables
%onsets and offset indices of all significant events
onset_offset = Events.onset_offset;
%1's demarking the entire event from onset to offset 
onset_ones = Events.onset_ones;
%1 where ever there is a significant Ca2+ onset
onset_binary =  Events.onset_binary;

%run events
run_onset_offset = Events.Run.run_onset_offset;
run_onset_binary = Events.Run.run_onset_binary;
run_onset_ones = Events.Run.run_onset_ones;

%no run events
norun_onset_offset = Events.NoRun.norun_onset_offset;
norun_onset_binary = Events.NoRun.norun_onset_binary;
norun_onset_ones = Events.NoRun.norun_onset_ones;

%lap trial labels - double check to make sure trials are not flipped
lap_order = Behavior.performance.trialOrder;

%number of events for each ROI
for rr =1:size(onset_offset,2)
    nbEvents(rr) = size(onset_offset{rr},1);
end

%number of events for each RUN during run epochs
for rr =1:size(run_onset_offset,2)
    nbEvents_run(rr) = size(run_onset_offset{rr},1);
end

%number of events for each NO RUN during run epochs
for rr =1:size(norun_onset_offset,2)
    nbEvents_norun(rr) = size(norun_onset_offset{rr},1);
end

%use Behavior{1}.resampled.lapNb to define frame start and stop indices
startLap = min(Behavior.resampled.lapNb);
endLap = max(Behavior.resampled.lapNb);

%break each variables into separate laps that are later used to define
%combined event variables according to which laps are selected
%for each lap
for ii=startLap:endLap
    %behavioral data
    %associated imaging indices for each lap
    lapFramesIdx{ii} = find(Behavior.resampled.lapNb == ii);
    lapTime{ii} = time(lapFramesIdx{ii});
    lapPosition{ii} = position(lapFramesIdx{ii});
    lapPositionNorm{ii} = position_norm(lapFramesIdx{ii});
    lapNb{ii} = ii*ones(1,size(lapFramesIdx{ii},1));
    
    %from resampled behavioral time behavioral time
    lapStartStopTime(ii,:) = [lapTime{ii}(1), lapTime{ii}(end)];
    
    %run epochs (run_ones split by lap)
    lapRunEpochs{ii} = run_ones(lapFramesIdx{ii});
    
    %select run frames for time position and lapNb (only run frames
    %included)
    lapTimeRun{ii} = lapTime{ii}(lapRunEpochs{ii} ==1);
    lapPositionRun{ii} = lapPosition{ii}(lapRunEpochs{ii} ==1);
    lapPositionNormRun{ii} = lapPositionNorm{ii}(lapRunEpochs{ii} ==1);
    lapFramesIdxRun{ii} = lapFramesIdx{ii}(lapRunEpochs{ii} ==1);
    lapNbRun{ii} = ii*ones(1,size(lapFramesIdxRun{ii},1));
    
    %imaging data
    lapImagingTime{ii} = imaging_time(lapFramesIdx{ii});
    lapTraces{ii} = imaging_trace(lapFramesIdx{ii},:);
    
    %event data
    onset_binary_lap{ii} = onset_binary(lapFramesIdx{ii},:);
    onset_ones_lap{ii} = onset_ones(lapFramesIdx{ii},:);
    
    %event run data
    run_onset_binary_lap{ii} = run_onset_binary(lapFramesIdx{ii},:);
    run_onset_ones_lap{ii} = run_onset_ones(lapFramesIdx{ii},:);
    
end

%% Create cell that assigns each event for each ROI to a given lap based on event onset frame

%do also for frames during run activity %run_onset_offset_laps{ii} 

%first extract all event onset/offset frames and convert to vector
eventsR = cell2mat(reshape(onset_offset,[],1));
%for run events (note: across all frames - not just run frames)
eventsR_run = cell2mat(reshape(run_onset_offset,[],1));
%for no run events (note: across all frames - not just run frames)
eventsR_norun = cell2mat(reshape(norun_onset_offset,[],1));

%take onsets and find associated lap
eventLapR = lapFrames(eventsR(:,1));
%for run events
eventLapR_run = lapFrames(eventsR_run(:,1));
%for no run events
eventLapR_norun = lapFrames(eventsR_norun(:,1));

%reshape back into cell that corresponds to onset_offset events cell
endEventIdx = cumsum(nbEvents);
startEventIdx = endEventIdx + 1;
startEventIdx = [1,startEventIdx(1:end-1)];
%for run events
endEventIdx_run = cumsum(nbEvents_run);
startEventIdx_run = endEventIdx_run + 1;
startEventIdx_run = [1,startEventIdx_run(1:end-1)];
%for No run events
endEventIdx_norun = cumsum(nbEvents_norun);
startEventIdx_norun = endEventIdx_norun + 1;
startEventIdx_norun = [1,startEventIdx_norun(1:end-1)];

%reshape back into cell that matched onset_offset cell
%assigns a lap # to each event
for rr =1:size(onset_offset,2)
    onset_offset_lapIdx{rr} = eventLapR(startEventIdx(rr):endEventIdx(rr));
end

%for run events - laps corresponding to each run event for each ROI
for rr =1:size(run_onset_offset,2)
    onset_offset_lapIdx_run{rr} = eventLapR_run(startEventIdx_run(rr):endEventIdx_run(rr));
end

%for No run events - laps corresponding to each run event for each ROI
for rr =1:size(norun_onset_offset,2)
    onset_offset_lapIdx_norun{rr} = eventLapR_norun(startEventIdx_norun(rr):endEventIdx_norun(rr));
end

%% Plot calcium trace and position for a sample ROI across all laps

%A trial laps (2 = correct, A trials)
lapSelect{1} = find(lap_order == 2);
%B trial laps (3 = correct, B trials)
lapSelect{2} = find(lap_order == 3);
%all laps
lapSelect{3} = (1:size(lap_order,1))';

%% Holder cells
%holder cell for the combined traces/positions
traces_split = cell(1,size(lapSelect,2));
traces_time_split = cell(1,size(lapSelect,2));

positions_split = cell(1,size(lapSelect,2));
positions_norm_split = cell(1,size(lapSelect,2));
time_split = cell(1,size(lapSelect,2));

%split for place cell analysis
onset_binaries_split = cell(1,size(lapSelect,2));
onset_ones_split = cell(1,size(lapSelect,2));

%run events
run_onset_ones_split = cell(1,size(lapSelect,2));
run_onset_binary_split = cell(1,size(lapSelect,2));
lapNb_split = cell(1,size(lapSelect,2));
run_ones_split = cell(1,size(lapSelect,2));

%% 
%for each set of laps
for ii=1:size(lapSelect,2)
    
    %split lap start/stop times
    lapStartStopTime_split{ii} = lapStartStopTime(lapSelect{ii},:);
    
    for jj=1:size(lapSelect{ii},1)
        
        %behavior
        positions_split{ii} = [positions_split{ii}; lapPosition{lapSelect{ii}(jj)}];
        positions_norm_split{ii} = [positions_norm_split{ii}; lapPositionNorm{lapSelect{ii}(jj)}];
        
        time_split{ii} = [time_split{ii}; lapTime{lapSelect{ii}(jj)}];
        lapNb_split{ii} = [lapNb_split{ii}, lapNb{lapSelect{ii}(jj)}];
        
        %imaging
        traces_split{ii} = [traces_split{ii}; lapTraces{lapSelect{ii}(jj)}];
        traces_time_split{ii} = [traces_time_split{ii}; lapImagingTime{lapSelect{ii}(jj)}];       
        
        %events
        onset_binaries_split{ii} = [onset_binaries_split{ii}; onset_binary_lap{lapSelect{ii}(jj)}];
        onset_ones_split{ii} = [onset_ones_split{ii}; onset_ones_lap{lapSelect{ii}(jj)}];
        run_onset_binary_split{ii} = [run_onset_binary_split{ii}; run_onset_binary_lap{lapSelect{ii}(jj)}];
        run_onset_ones_split{ii} = [run_onset_ones_split{ii}; run_onset_ones_lap{lapSelect{ii}(jj)}];
        run_ones_split{ii}  = [run_ones_split{ii}; lapRunEpochs{lapSelect{ii}(jj)}];
    end
    
    %extract run time, position and labNb for each lap
    run_position_split{ii} = positions_split{ii}(run_ones_split{ii}==1);
    run_position_norm_split{ii} = positions_norm_split{ii}(run_ones_split{ii}==1);
    run_time_split{ii} = time_split{ii}(run_ones_split{ii}==1);
    run_lapNb_split{ii} = lapNb_split{ii}(run_ones_split{ii}==1);
    
end

%% Put lap start/stop times into cells for input to place cell code

for ii=1:size(lapStartStopTime_split,2)
    %for each lap in trial split array
    for ll=1:size(lapStartStopTime_split{ii},1)
        lapStartStopTime_split_cell{ii}{ll} = lapStartStopTime_split{ii}(ll,:);
    end
end

%% Select onset_offset frames for each event for each ROI to corresponding laps 

%select only event indices associated with lap of interest
for ii=1:size(lapSelect,2)
    %all events that associated with given trial laps (all)
    onset_offset_select{ii} = ismember(eventLapR,lapSelect{ii});
    
    %all events that associated with given trial laps (run only)
    onset_offset_select_run{ii} = ismember(eventLapR_run,lapSelect{ii});
        
    %all events that associated with given trial laps (no run only)
    onset_offset_select_norun{ii} = ismember(eventLapR_norun,lapSelect{ii});
end


%take out selected events based on the filtered indices above
for ii=1:size(lapSelect,2)
    
    %for all events, select the laps associated with specifc trial type
    for rr =1:size(onset_offset,2)
        onset_offset_selectIdx{ii}{rr} = onset_offset_select{ii}(startEventIdx(rr):endEventIdx(rr));
        onset_offset_split{ii}{rr} = onset_offset{rr}(onset_offset_selectIdx{ii}{rr},:);
    end
    
    %for run events, select the laps associated with specifc trial type
    for rr =1:size(run_onset_offset,2)
        onset_offset_selectIdx_run{ii}{rr} = onset_offset_select_run{ii}(startEventIdx_run(rr):endEventIdx_run(rr));
        onset_offset_run_split{ii}{rr} = run_onset_offset{rr}(onset_offset_selectIdx_run{ii}{rr},:);
    end
    
    %for No run events, select the laps associated with specifc trial type
    for rr =1:size(norun_onset_offset,2)
        onset_offset_selectIdx_norun{ii}{rr} = onset_offset_select_norun{ii}(startEventIdx_norun(rr):endEventIdx_norun(rr));
        onset_offset_norun_split{ii}{rr} = norun_onset_offset{rr}(onset_offset_selectIdx_norun{ii}{rr},:);
    end
    
end

%% Assemble structures for each sets of laps/trials for input to place cell analysis

for ii=1:size(lapSelect,2)
    %Position, time, and lap # of resampled+restricted behavior
    Behavior_split{ii}.resampled.position = positions_split{ii};
    Behavior_split{ii}.resampled.position_norm = positions_norm_split{ii};
    Behavior_split{ii}.resampled.time = time_split{ii};
    Behavior_split{ii}.resampled.lapNb = lapNb_split{ii}';
    Behavior_split{ii}.restricted.lap = lapStartStopTime_split_cell{ii};
    %copy the same info outside restricted for spiral event code
    Behavior_split{ii}.lap = lapStartStopTime_split_cell{ii};
    
    %Run behavior
    Behavior_split{ii}.resampled.run_time = run_time_split{ii};
    Behavior_split{ii}.resampled.run_position = run_position_split{ii};
    Behavior_split{ii}.resampled.run_position_norm = run_position_norm_split{ii};
    Behavior_split{ii}.resampled.run_lapNb = run_lapNb_split{ii}';
    
    %run activity across all restricted laps
    Behavior_split{ii}.run_ones = run_ones_split{ii};
    
    %Imaging of restricted behavior (add non-restricted in the future)
    Imaging_split{ii}.trace_restricted = traces_split{ii};
    Imaging_split{ii}.time_restricted  = traces_time_split{ii};
    %imaging period (s)
    Imaging_split{ii}.dt = Imaging.dt;
    
    %Events
    Events_split{ii}.Run.run_onset_ones = run_onset_ones_split{ii};
    Events_split{ii}.Run.run_onset_binary = run_onset_binary_split{ii};
    Events_split{ii}.Run.run_onset_offset = onset_offset_run_split{ii};
    
    Events_split{ii}.NoRun.norun_onset_offset = onset_offset_norun_split{ii};
    
    Events_split{ii}.onset_offset = onset_offset_split{ii};
    Events_split{ii}.options = Events.options;
    
    %add time vector 
    Events_split{ii}.dff_time = Events.dff_time;
    
    
end

%% split by lap (all laps/trials) for spiral map input

%behavior
Behavior_split_lap.position = lapPosition;
Behavior_split_lap.Run.position = lapPositionRun;
Behavior_split_lap.time = lapTime;
Behavior_split_lap.Run.time = lapTimeRun;
Behavior_split_lap.run_ones = lapRunEpochs;
Behavior_split_lap.lap = Behavior.lap;

%events
Events_split_lap.Run.run_onset_binary = run_onset_binary_lap;
Events_split_lap.onset_binary = onset_binary_lap;

%% Break run epochs into start/stop frames in A or B trials

%do this for each split
for ii=1:size(Behavior_split,2)
    
    % Run on/off idx matrix
    %on/off indices within each trial set
    %pad with 0 at each end
    run_on_idx{ii} = find([0;diff([0;Behavior_split{ii}.run_ones;0])] == 1)-1;
    run_off_idx{ii} = find([0;diff([0;Behavior_split{ii}.run_ones;0])] == -1)-2;
    
    %create a run on/off matrix of indices
    %work around the onset indices
    run_on_off_idx{ii} = [run_on_idx{ii}, run_off_idx{ii}];
    
    % No run on/off idx matrix
    
    norun_on_idx{ii} = find([0;diff([0;~Behavior_split{ii}.run_ones;0])] == 1)-1;
    norun_off_idx{ii} = find([0;diff([0;~Behavior_split{ii}.run_ones;0])] == -1)-2;
    
    %create a run on/off matrix of indices
    %work around the onset indices
    norun_on_off_idx{ii} = [norun_on_idx{ii}, norun_off_idx{ii}];
    
end

%QC check of on/off
figure;
subplot(2,1,1)
hold on
title('Run on intervals (correct A and B)');
plot(Behavior_split{ii}.run_ones, 'k')
stem(run_on_idx{ii}, 2*ones(size(run_on_idx{ii},1)),'g');
stem(run_off_idx{ii}, 2*ones(size(run_off_idx{ii},1)),'r');
hold off

subplot(2,1,2)
hold on
title('Run off intervals (correct A and B)');
plot(~Behavior_split{ii}.run_ones, 'k')
stem(norun_on_idx{ii}, 2*ones(size(norun_on_idx{ii},1)),'g');
stem(norun_off_idx{ii}, 2*ones(size(norun_off_idx{ii},1)),'r');
hold off

%reconstruction QC - RESUME HERE
%create blank
for ii=1:size(Behavior_split,2)
    %create blank for each type
    run_ones_recon_run_on{ii} = false(size(Behavior_split{ii}.run_ones,1),1);
    
    for jj=1:size(run_on_off_idx{ii},1)
        run_ones_recon_run_on{ii}(run_on_off_idx{ii}(jj,1):run_on_off_idx{ii}(jj,2)) = 1;
    end
    
    for jj=1:size(norun_on_off_idx{ii},1)
        run_ones_recon_run_on{ii}(norun_on_off_idx{ii}(jj,1):norun_on_off_idx{ii}(jj,2)) = 0;
    end
end

for ii=1:size(Behavior_split,2)
check_run_intervals(ii) = isequal(Behavior_split{ii}.run_ones,run_ones_recon_run_on{ii});
end

%display that check passes
if sum(check_run_intervals) == size(Behavior_split,2)
    disp('Run intervals for each sub-split agree.');
end


%% Plot a sample ROI from each trial type
ROI = randi(size(onset_offset,2));
ROI = 401;

%plot all laps
figure
%A laps
subplot(2,1,1)
hold on
title('A trials')

%plot whole A trace as continuous frames
plot(traces_split{1}(:,ROI),'k')

%shade traces that are in run epochs green
for ii=1:size(run_on_off_idx{1},1)
    plot(run_on_off_idx{1}(ii,1):run_on_off_idx{1}(ii,2),...
        traces_split{1}(run_on_off_idx{1}(ii,1):run_on_off_idx{1}(ii,2),ROI),'Color', [0 100 0]/255);
end
%shade traces that are in not in run epochs red
for ii=1:size(norun_on_off_idx{1},1)
    plot(norun_on_off_idx{1}(ii,1):norun_on_off_idx{1}(ii,2),...
        traces_split{1}(norun_on_off_idx{1}(ii,1):norun_on_off_idx{1}(ii,2),ROI),'r');
end

subplot(2,1,2)
hold on
title('B trials')

%plot whole A trace as continuous frames
plot(traces_split{2}(:,ROI),'k')

%shade traces that are in run epochs green
for ii=1:size(run_on_off_idx{2},1)
    plot(run_on_off_idx{2}(ii,1):run_on_off_idx{2}(ii,2),...
        traces_split{2}(run_on_off_idx{2}(ii,1):run_on_off_idx{2}(ii,2),ROI),'Color', [0 100 0]/255);
end
%shade traces that are in not in run epochs red
for ii=1:size(norun_on_off_idx{2},1)
    plot(norun_on_off_idx{2}(ii,1):norun_on_off_idx{2}(ii,2),...
        traces_split{2}(norun_on_off_idx{2}(ii,1):norun_on_off_idx{2}(ii,2),ROI),'r');
end


%% plot the split laps for visual confirmation

figure;
subplot(2,1,1)
hold on
title('A trials');
plot(traces_split{1}(:,ROI),'k')
plot((positions_split{1}/200) -1.2, 'r')
ylim([-2 3]);
xlim([0 size(traces_split{1},1)]);
hold off

%B laps
subplot(2,1,2)
hold on
title('B trials');
plot(traces_split{2}(:,ROI),'k')
plot((positions_split{2}/200) -1.2, 'r')
ylim([-2 3]);
xlim([0 size(traces_split{2},1)]);
hold off

end

