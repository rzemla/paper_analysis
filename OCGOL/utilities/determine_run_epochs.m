function [Events, Behavior] = determine_run_epochs(Events,Behavior,Imaging, updated_dff, options)

%% Set parameters and variables

%whether to run the analysis on only full laps
%options.restricted=Events.options.restricted;

%minimum duration of running epoch
mindur = options.mindur;

%minimum speed of running epock
minspeed = options.minspeed;

%time below which neighboring run epochs are merged 
mergdur = options.merge;

%width of moving mean filter for Cossart (also filter for speed)
mov_wind = options.moving_window;

%Time and behavior
%if only doing analysis on restricted laps, set the following calcium frame
%time, behavior time cumulative position and lap start relative position
%all behavior variable should match in sample length
%C_df time should match against C_df time legnth
if options.restrict==true
    
    %original dF/F signal
    C_df=Imaging.trace_restricted;
    
    %restricted updated dF/F signal
    updated_C_df = updated_dff.F_dff_exp;
    
    %imaging time (restricted)
    Cdf_time=Imaging.time_restricted;
    %behavior time (restricted)
    time=Behavior.restricted.time;
    %cumulative position (restricted)
    cum_position=Behavior.restricted.cumulativeposition;
    %normalizesd position
    norm_position = Behavior.restricted.normalizedposition; 
    %relative position to start of lap (cm)
    position=Behavior.restricted.position;
    %index of lap at each behavior index
    lapNb = Behavior.restricted.lapNb;
    
elseif options.restrict==false
    %original dF/F signal
    C_df=Imaging.trace;
    
    %non-restricted updated dF/F signal
    updated_dff = F_dff_exp;
    %imaging time
    Cdf_time=Imaging.time;
    %behavior time
    time=Behavior.time;
    %cumulative position
    cum_position=Behavior.cumulativeposition;
    %position since lap start
    position=Behavior.position;
    %lap nb - get unrestricted version here
    lapNb = Behavior.restricted.lapNb;
    
end

%% Load resampled behavior variables

%resample behavioral time
res_time = Behavior.resampled.time;

%resampled behavioral position
res_position = Behavior.resampled.position;

%cumulative positon for speed calculation
res_cum_position = Behavior.resampled.cumulativeposition;

%normalized position
res_norm_position = Behavior.resampled.normalizedposition;

%resampled lap indicator
res_lapNb = Behavior.resampled.lapNb;

%% Import detected events

%onset and offset frame of each event for each ROI
onset_offset=Events.onset_offset;

%mark all events for each ROI in time with 1's (for the duration of the
%event
onset_ones=Events.onset_ones;

%mark only the onset of the events in time with 1's (just onset)
onset_binary=Events.onset_binary;

%% Speed of animal and 

% Mean measured framerate of imaging (Hz - frames/second)
avg_fr = 1/Imaging.dt;

%Speed
speed=[0;diff(res_cum_position)]*avg_fr;

%average speed using mean filter
speed = movmean(speed,[mov_wind]);


%% Find running epochs

switch options.method
    case 'peak' %(Danielson 2016)
        disp('minimum peak method')
        
        %set minimum peak speed
        pks_thr = options.minpeak; 
        
        %Find periods of forward motion (speed >0)
        %list the indices where the speed is greater than 0
        %will need to this to reconstruct 1's matrix
        run_idx = find(speed > minspeed);
        
        %extract the corresponding frames times when the speed is greater than 0
        run_time = Cdf_time(run_idx);
        
        %extract the corresponding speeds when the speed is greater than 0
        run_speed = speed(run_idx);
        
        %find time difference between each frame time when the animal was
        %running
        dist_epochs=diff(run_time);
        
        %% MERGING CRITERION IS SELECTED FOR HERE
        %Find epochs separated by more than the merging threshold
        
        %find all run epochs separated by more than the minimum time
        %(indices)
        %get a list of indices where the epochs between running is greater
        %than 0.5s (or whatever minimum)
        %not meeting the criteria
        epochs_end_idx=find(dist_epochs>=mergdur);
        
        %frame indices where any animal was running are taken into consideration(all frames are running
        %frames)
        %start of first running epoch
        %start of next running epochs which are the frames immediately
        %after the end 
        %+ 1 frames after each epoch where the not run epoch was too long
        %(and did not get merged)
        epochs_start_idx=[run_idx(1);epochs_end_idx+1];
        
        %frames where epochs are greater than minmerge apart ;last run frame
        epochs_end_idx=[epochs_end_idx ;length(run_time)];
        
        %take start and end epoch frames and merge into one matrix where
        %each row corresponds to a run epoch with correponding start and
        %stop frame AMONG the RUN frames
        run_epochs_idx=[epochs_start_idx epochs_end_idx];
        
        %extract the corresponding frame times in the original input
        run_epochs_time=run_time(run_epochs_idx);
        
        
        %% MIMIMUM DURATION OF EPOCH IS SELECTED FOR HERE
        %Minimum duration for running epoch
        %duration of each run epoch in imaging time domain
        run_epochs_dur=run_epochs_time(:,2)-run_epochs_time(:,1);
        
        %for each run epoch
        for i=1:size(run_epochs_dur,1)
            
            %if duration is less than mindur for that epoch - NaN the start
            %and end frames of that epoch and the associated imaging times
            if run_epochs_dur(i)<mindur ==1
                
                run_epochs_time(i,:)=NaN;
                run_epochs_idx(i,:)=NaN;
            end
        end
        
        %extract only corresponding frame times that meet the minimum
        %duration criteria - use epoch end times
        run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
        
        %extract only corresponding frame indices that meet the minimum
        %duration criteria
        run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);
        
        %% PEAK SPEED CRITERION IS MATCHED HERE
        %Find if peaks speed value in running epochs, count how many:
        %how many epochs that are above 5 cm/s
        
        
        %for each run epoch that met above criteria
        for i=1:size(run_epochs_idx,1)
            
            %get the corresponding speed between those frame indices when
            %speed was greater than 0 into separate cells (velocity across
            %each epoch)
            speed_run_epochs{i}=run_speed(run_epochs_idx(i,1):run_epochs_idx(i,2));
            
            %at how many frames did the speed exceed the minimum threshold
            pks_speed_cnt(i)=sum(speed_run_epochs{i}>pks_thr);
            
            %Remove run epochs if no peak speed
            %NaN the corresponding epoch frame indices and frame times at
            %which no single sample passes the peak speed
            if pks_speed_cnt(i)==0
                run_epochs_time(i,:)=NaN;
                run_epochs_idx(i,:)=NaN;
            end
        end
        
        %extract only the associated times and corresponding indices that
        %meet the minimum peak speed criterion
        run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
        
        run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);
        
        %norun_epochs_idx =
        
    case 'speed' %Cossart approach
        disp('average speed')
        run_thr=options.minspeed;
        
        %Find periods when speed > threshold
        
        %find frames where speed exceeds the threshold (unlike Danielson
        %where frames with speeds >0 are taken
        run_idx=find(speed>run_thr);
        
        %get corresponding frame times
        run_time=Cdf_time(run_idx);
        
        %get corresponding speed in those times
        run_speed=speed(run_idx);
        
        %find the time between corresponding epochs
        dist_epochs=diff(run_time);
        
        %Find epochs separated by more than the merging threshold
        %will become the endpoint indices
        epochs_end_idx=find(dist_epochs>=mergdur);
        
        %first run frame and +1 frame to epochs end idx form start frames
        %of epochs
        epochs_start_idx=[run_idx(1);epochs_end_idx+1];
        
        %where epochs separated by more than minimum duration for the end
        %and last frame of running frames
        epochs_end_idx=[epochs_end_idx ;length(run_idx)];
        
        %combine start and stop run frames together
        run_epochs_idx=[epochs_start_idx epochs_end_idx];
        
        %extract corresponding imaging times for these epochs
        run_epochs_time=run_time(run_epochs_idx);
        
        
        %MINIMUM DURATION CRITERION - RZ
        %Minimum duration for running epoch
        
        %calculate duration of each epoch (based on frame time)
        run_epochs_dur=run_epochs_time(:,2)-run_epochs_time(:,1);
        
        %for each run epoch
        for i=1:size(run_epochs_dur,1)
            
            %NaN those which are less than the minimum duration
            if run_epochs_dur(i)<mindur ==1
                run_epochs_time(i,:)=NaN;
                run_epochs_idx(i,:)=NaN;
            end
        end
        
        %remove corresoning frame idxs and times that do not meet minimum
        run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
        
        run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);
        
    otherwise
        disp('options.method should be peak or speed')
end

%% Parse behavioral time and postion based on run epochs

%Run idx - frames in the 'run domain' that correspond to whole frame domain
%start and stop frames of filtered run epochs filtered that correspond to
%the actual input frames - RZ
run_int_idx=run_idx(run_epochs_idx);

%make binary

%make vector filled with 0's that corresponds to length of imaging session
%this will be used for the marking ONLY the onsets of run epochs
run_binary=zeros(size(Cdf_time,1), 1);

%make copy to run_ones - 1's correponds to run epochs
%1's make the all the frames when the animal is in a run epoch
run_ones=run_binary;

%take the frames of run epoch onsets and set to 1
run_binary(run_int_idx(:,1),1)=1;

%make ones

%for each run epoch
for i=1:size(run_int_idx,1)
    
    %for each frame between the start and end of each run epoch fill the 
    %frames in between with 1's
    run_ones(run_int_idx(i,1):run_int_idx(i,2),:)=1;
end

%on resample behavior data matched to sampling rate of imaging recording
%get the behavioral time of the animal during run epochs
runtime=res_time(run_ones==1);

%get the behavioral position of the animal during run epochs
run_position=res_position(run_ones==1);

%get resampled lap index
run_lapNb = res_lapNb(run_ones==1);

%get the event onsets ONLY during run epochs
%on2=onset_binary(run_ones==1,:);


%% Generate no run index matrix here

%start and stop indices of run intervals across all imaging time
%lengh of imaging ses.
time_len =  size(Cdf_time,1);

%check if  first run epoch starts on 1st position of imaging
if run_int_idx(1,1) == 1
    for ii = 1:(size(run_int_idx,1)-1)
        norun_int_idx(ii,:) = [run_int_idx(ii,2)+1, run_int_idx(ii+1,1)-1];
    end
    
    %if the last run index not equal to end of imaging time length
    if run_int_idx(end,2) ~= time_len
        %take last index of run + 1 until end and add to matrix
        norun_int_idx = [norun_int_idx; [(run_int_idx(end,2)+1), time_len]];
    end
    
    %if not equal to 1
elseif run_int_idx(1,1) ~= 1
    %fill in first non run epochs
    norun_int_idx(1,:) = [1, run_int_idx(1,1)-1];
    %fill the rest of no run intervals (starting at 2nd position)
    for ii = 1:(size(run_int_idx,1)-1)
        norun_int_idx(ii+1,:) = [run_int_idx(ii,2)+1, run_int_idx(ii+1,1)-1];
    end
    
        %if the last run index not equal to end of imaging time length
        %(same as above)
    if run_int_idx(end,2) ~= time_len
        %take last index of run + 1 until end and add to matrix
        norun_int_idx = [norun_int_idx; [(run_int_idx(end,2)+1), time_len]];
    end
    
end

%QC check that intervals in matrices are properly constructed
%reconstruct run_ones using the defined run and no run intervals and check
%that it is the same
%generate blank
run_ones_recon = false(size(run_ones,1),1);

%insert the run intervals
for ii=1:size(run_int_idx,1)
run_ones_recon(run_int_idx(ii,1):run_int_idx(ii,2)) = 1;
end

%insert the no run intervals
for ii=1:size(norun_int_idx,1)
run_ones_recon(norun_int_idx(ii,1):norun_int_idx(ii,2)) = 0;
end

if isequal(run_ones, run_ones_recon)
    disp('Reconstructed run_ones interval agrees with original.');
else
    disp('Reconstructed run_ones interval does NOT agree with original.');
end
%print



%% Parse events based on running epochs

% restrict onset to running epochs

%for each ROI with associated event onset and offsets (frames)
for i=1:size(onset_offset,2)
    
    %for each run epoch
    for ii=1:size(run_int_idx,1)
        
        %if at least 1 event occured
        if isempty(onset_offset{i})==0
            
            %take all the frames of onset for that ROI and store in temp
            %variable onset
            %can place this outside this loop since it only depends on i
            %onset events during all epochs (not just run at this point)
            onset=onset_offset{i}(:,1);
            
            %for that ROI, give logical output for each event if that event
            %is within the range of run epoch frames for that run epoch
            
            %a given event should be only within 1 run epoch 
            %each column is the corresponding run epoch
            %yields a column vector on each interation where each row
            %corresponds to whether that event onset occurred in that epoch
            %- RZ
            on_R{i}(:,ii)= onset >= run_int_idx(ii,1) & onset <= run_int_idx(ii,2);
            
            %sum along all columns (run epochs) to see in which events
            %occurred in a run epoch - get a column vector which indicates
            %whether that event occurred in a run epoch or not (does not
            %indicate which one in occurred in) - RZ
            %this will act as a vector telling which events should be kept
            %and which discarded for run events output - RZ
            onset_R_keep{i}=sum(on_R{i},2);
        
            %if that ROI has no onsets, the the keep vector to NaN
        elseif isempty(onset_offset{i})==1
            onset_R_keep{i}=nan;
        end
    end
end

%for each ROI and associated onset and offset frames
for i=1:size(onset_offset,2)
    
    %take the binary vector of which events are within run epochs and
    %place in temporary variable (column vector)
    keep=onset_R_keep{i};
    
    %keep onset and offset frames of each event that are within a run epoch
    %(set keep to a logical by equating to 1)
    run_onset_offset{i}=onset_offset{i}(keep==1,:);
    
    %keep onset and offset frames of each event that are outside a run
    %epoch
    norun_onset_offset{i}=onset_offset{i}(keep==0,:);
end


%% Make binary outputs

%Make binary

%frames x ROI - make empty 0's matrix (for run events)
run_onset_binary=zeros(size(onset_binary,1), size(onset_binary,2));

%copy 0's matrix for norun events
norun_onset_binary=run_onset_binary;

%for each ROI
for i=1:size(onset_offset,2)
    
    %if at least one event for that ROI
    if isempty(onset_offset{i})==0
        
        %take the run onset frames and for that column (frames for that ROI)
        %set to 1 at those onsets 
        run_onset_binary(run_onset_offset{i}(:,1),i)=1;
        
        %take the norun onset frames and for that column (frames for that ROI)
        %set to 1 at those onsets 
        norun_onset_binary(norun_onset_offset{i}(:,1),i)=1;
    end
end

%Make ones

%frames x ROI (zeros from run onset until end of event) for run events
run_onset_ones=zeros(size(onset_ones,1), size(onset_ones,2));

%copy of 0's matrix for norun events
norun_onset_ones=run_onset_ones;

%for each ROI
for i=1:size(onset_offset,2)
    
    %for each run event and its associated onset and offset frames
    for irun=1:size(run_onset_offset{i},1)
        
        %for all frames for the ROI, set to 1 all values within that range
        run_onset_ones(run_onset_offset{i}(irun,1):run_onset_offset{i}(irun,2),i)=1;
    end
    
        %for each norun event and its associated onset and offset frames
    for inorun=1:size(norun_onset_offset{i},1)
        
        %for all frames for the ROI, set to 1 all values within that range
        norun_onset_ones(norun_onset_offset{i}(inorun,1):norun_onset_offset{i}(inorun,2),i)=1;
    end
end


%% Save variables to Behavior structure

%for all frames, which ones occur during run epochs
Behavior.run_ones=run_ones;

%split into cell of indices for plotting purposes (avoid discontinuities)
%Behavior.norunbinary=norun;

%speed of the animal at each frame - spikes because of lap difference - but
%excluded during run epoch parsing
Behavior.speed = speed;

%corresponding behavioral time matched to imaging frame rate (resampled) during run
%epochs 
Behavior.resampled.run_time = runtime;
%Behavior.resampled.no_run_time=noruntime;

%position during run epochs
Behavior.resampled.run_position = run_position;

%running epoch intervals in imaging time space
Behavior.resampled.run_idx = run_int_idx;

%no run epochs intervals in imaging time space
Behavior.resampled.norun_idx = norun_int_idx;

%lap idx during running
Behavior.resampled.run_lapNb = run_lapNb;

%% Save variables to Events structure 

%Run events

%onset and offset frames of events in input frames space during running
Events.Run.run_onset_offset=run_onset_offset;

%onset events in input framexROI space
Events.Run.run_onset_binary=run_onset_binary;

%span of onset - offset events in input framexROIspace
Events.Run.run_onset_ones=run_onset_ones;


%No run events -  same as above except for norun events

Events.NoRun.norun_onset_offset=norun_onset_offset;

Events.NoRun.norun_onset_binary=norun_onset_binary;

Events.NoRun.norun_onset_ones=norun_onset_ones;

%carry over the options
Events.options.run_epochs=options;

%% Figure
if options.dispfig==true
    ROI=options.c2plot;

    figure;
    subplot(2,1,1);
    hold on;
    title('Marked events during respective epochs')
    yyaxis left
    ylim([-1.5 2.5])
    ylabel('dF/F')
    %all events
    onset_idx =find(onset_binary(:,ROI) == 1);
    %plot all events
    stem(Cdf_time(onset_idx), 2*onset_binary(onset_idx,ROI),'k--');
    %find run events
    run_onset_idx =find(run_onset_binary(:,ROI) == 1);
    %mark run events by filling in stem circle
    scatter(Cdf_time(run_onset_idx), 2*run_onset_binary(run_onset_idx,ROI),'*','MarkerEdgeColor', [0 100 0]/255);
    %find no run events
    norun_onset_idx =find(norun_onset_binary(:,ROI) == 1);
    %mark run events by filling in stem circle
    scatter(Cdf_time(norun_onset_idx), 2*norun_onset_binary(norun_onset_idx,ROI),'*','MarkerEdgeColor', 'r');
    
    %plot run interval parts of the trace as green
    for ii=1:size(run_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(run_int_idx(ii,1):run_int_idx(ii,2)),C_df(run_int_idx(ii,1):run_int_idx(ii,2),ROI),'LineStyle', '-',...
            'Marker','none','Color', [0 100 0]/255)
        %p1.LineStyle = '-';
    end
    %plot no run interval parts of the trace as red
    for ii=1:size(norun_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(norun_int_idx(ii,1):norun_int_idx(ii,2)),C_df(norun_int_idx(ii,1):norun_int_idx(ii,2),ROI), 'Color','r',...
            'LineStyle','-','Marker','none')
    end
    
    
    %plot normalized position below _shift
    %plot(Cdf_time, res_norm_position - 1.2, 'k');
    for ii=1:size(run_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(run_int_idx(ii,1):run_int_idx(ii,2)),res_norm_position(run_int_idx(ii,1):run_int_idx(ii,2))-1.2,'LineStyle', '-',...
            'Marker','none','Color', [0 100 0]/255)
    end
    for ii=1:size(norun_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(norun_int_idx(ii,1):norun_int_idx(ii,2)),res_norm_position(norun_int_idx(ii,1):norun_int_idx(ii,2))-1.2, 'Color','r',...
            'LineStyle','-','Marker','none')
    end
    
    %mark run and no run events on position plot
    
    
    yyaxis right
    ylim([-1.5 2.5])
    yticks([-1.2 -0.2])
    yticklabels({'0','1'});
    ylabel('Normalized position');
    hold off
    
    %%%%
    %no run epochs
    s2 = subplot(2,1,2);
    hold on;
    yyaxis left
    title('Speed vs. Position');
    plot(Cdf_time, speed,'LineStyle','-','Marker','none','Color','k')
    xlabel('Time [s]');
    ylabel('Speed [cm/s]');


    yyaxis right
    ylim([-0.2 1.2])
    yticks([0 1])
    yticklabels({'0','1'});
    ylabel('Normalized position');
        %plot normalized position below _shift
    %plot(Cdf_time, res_norm_position , 'k');
    for ii=1:size(run_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(run_int_idx(ii,1):run_int_idx(ii,2)),res_norm_position(run_int_idx(ii,1):run_int_idx(ii,2)),'LineStyle', '-',...
            'Marker','none','Color', [0 100 0]/255)
    end
    for ii=1:size(norun_int_idx,1)
        %plot trace in green that is run
        plot(Cdf_time(norun_int_idx(ii,1):norun_int_idx(ii,2)),res_norm_position(norun_int_idx(ii,1):norun_int_idx(ii,2)), 'Color','r',...
            'LineStyle','-','Marker','none')
    end

    
end

end



