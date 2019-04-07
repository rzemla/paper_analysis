function [Events, Behavior] = resample_and_determine_runV2(Events,Behavior,Imaging, options)

%% Set parameters and variables

%whether to run the analysis on only full laps
%options.restricted=Events.options.restricted;

%minimum duration of running epoch
mindur = options.mindur;

%time below which neighboring run epochs are merged 
mergdur = options.merge;

%width of moving mean filter for Cossart
mov_wind = options.moving_window;

%Select whether to lap restricted (complete laps) imaging data or all
%recorded data
if options.restrict==true
    C_df=Imaging.trace_restricted;
elseif options.restrict==false
    C_df=Imaging.trace;

end

%Time and behavior
%if only doing analysis on restricted laps, set the following calcium frame
%time, behavior time cumulative position and lap start relative position
%all behavior variable should match in sample length
%C_df time should match against C_df time legnth
if options.restrict==true
    %imaging time (restricted)
    Cdf_time=Imaging.time_restricted;
    %behavior time (restricted)
    time=Behavior.restricted.time;
    %cumulative position (restricted)
    cum_postion=Behavior.restricted.cumulativeposition;
    %relative position to start of lap (cm)
    position=Behavior.restricted.position;
    %index of lap at each behavior index
    lapNb = Behavior.restricted.lapNb;
    
elseif options.restrict==false
    Cdf_time=Imaging.time;
    time=Behavior.time;
    cum_postion=Behavior.cumulativeposition;
    position=Behavior.position;
    lapNb = Behavior.restricted.lapNb;
end

%Import the associatied Events - all should matched restricted calcium
%all events are restricted

%onset and offset frame of each event for each ROI
onset_offset=Events.onset_offset;

%mark all events for each ROI in time with 1's (for the duration of the
%event
onset_ones=Events.onset_ones;

%mark only the onset of the events in time with 1's (just onset)
onset_binary=Events.onset_binary;

%put behavioral time and position into 1 matrix together
time_position=[time position];

%% Compute speed - resample behavioral data to match imaging frequency
%resample - change resolution while keeping the size the same! - RZ

%Resample behavior at same frequency than imaging
%downsample to match imaging frequency (10 kHz --> ~30 Hz

%take all the imaging times and see how many counts fall within each
%behavioral 
%bin will be the same size as Cdf_time - represents ths bin indices

%see which bins of behavioral time (fed raw - higher sample) does the lower
%sampled imaging time fall into (Cdf_time)
%N is bin counts - how many events fall between each edges of time
%bin take each imaging time value and get corresponding bin that it
%corresponds to in the time (behavioral) vector
%histc return N that is the same size as time (edges) - should be 1 less
%which histcount fixed
%idea is that the imaging time will only fall within 1 bin of the much
%higher sampled signal (behavior)
%tells which in which behavioral bin (timepoint) each imaging timepoint
%belongs
%time (behavioral) - edges of bins
%use discretize since histc is no longer recommended
%[N,bin]=histc(Cdf_time,time);

%how many counts found in each bin
%does not work like histc
%[N2,edges,bin2] = histcounts(Cdf_time,time);

%this yields the same result as histc
%edges = time
%bin2 - indices of bins in time (behavior - time range) that contain a
%given imaging time point in
[bin,edges] = discretize(Cdf_time,time);


%If missing behavior time - should never run, but valid check
%should have an assignment to a behavioral timepoint
%for all the imaging timepoints
for idx=1:length(bin)
    %make sure that it has a corresponding behavioral recording timpoints
    %(should have if both recorded and no behavioral data lost- RZ)
    if bin(idx)==0
        disp('!!! Missing behavior recording data !!!')
        %if missing, assign the earlier corresponds bin to present one
        bin(idx)=bin(idx-1);
    end
end

%create index variable which is 1 more than the the corresponding higher sampled bin in
%bin
%NOT USED - RZ
index=bin+1;

%QC
%compare the imaging timestamps to the resampled behavior timestamps

%take the actual imaging time and correponding resampled behavior time and
%take the difference

%assign the resamples values into separate vectors/matrices that correspond
%to the length of acquired imaging frames
%does the same regardless assign from changing the indexing - fractional
%differences
if abs(Cdf_time-time(bin))<abs(Cdf_time-time(bin+1))

    index=bin; %%shifts index back
    
    res_time_position=time_position(bin,:);
    res_cum_postion=cum_postion(bin);
    res_position=position(bin);
    
    %resampled lap idx
    res_lapNb = lapNb(bin);
    
else %this will run all the time and 
    %take the corresponding resampled bin in the behavioral voltage and assign to
    %imaging timeframe
    
    %combined time and position (start offset)
    res_time_position=time_position(bin,:);
    
    %cumulative position (without offset)
    res_cum_position=cum_postion(bin);
    
    %start offset cumulative position (start offset)
    res_position=position(bin);
    
    %resampled lap idx
    res_lapNb = lapNb(bin);
end

%resample time
res_time=res_time_position(:,1);

%Smooth position
%res_cum_pos_sm=smooth(res_cum_postion,3);

%Calculate this value in the imaging code and import here
%average framerate will differ when recording session is split apart for
%OCGOL

% Mean measured framerate of imaging (Hz - frames/second)
%RZ - changed from mean to median given how the laps are split
%avg_fr=1/(median(diff(Cdf_time)));
avg_fr = 1/Imaging.dt;

%Speed - transition points have high velocties - calculate this on
%non-split version
speed=[0;diff(res_cum_position)]*avg_fr;

%average using moving mean - overwrites previous using using of moving mean
%window of length defined in mov_wind

speed = movmean(speed,[mov_wind]);



%% Find running epochs

switch options.method
    case 'peak' %(Danielson 2016)
        disp('minimum peak method')
        
        %set minimum peak speed
        pks_thr=options.minpeak; 
        
        %Find periods of forward motion (speed >0)
        %list the indices where the speed is greater than 0
        %will need to this to reconstruct 1's matrix
        run_idx=find(speed>0);
        
        %extract the corresponding frames times when the speed is greater than 0
        run_time=Cdf_time(run_idx);
        
        %extract the corresponding speeds when the speed is greater than 0
        run_speed=speed(run_idx);
        
        %find time difference between each frame time when the animal was
        %running
        dist_epochs=diff(run_time);
        
        %Find epochs separated by more than the merging threshold
        %MERGING CRITERION IS SELECTED FOR HERE - RZ
        
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
        
        
        %MIMIMUM DURATION OF EPOCH IS SELECTED FOR HERE
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
        
        %PEAK SPEED CRITERION IS MATCHED HERE - RZ
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
        
        
    case 'speed'
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


%% Find run time and position - extract the original frames that correspond to 
%run criteria filtered frames - RZ

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
on2=onset_binary(run_ones==1,:);


%% Running events

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

%% Make structure

%for all frames, which ones occur during run epochs
Behavior.run_ones=run_ones;
%Behavior.norunbinary=norun;

%speed of the animal at each frame - spikes because of lap difference - but
%excluded during run epoch parsing
Behavior.speed=speed;


%resampled behavior

%runtime (during run epochs) ~= run_time (speed greater > 0 --> more frames)

%corresponding behavioral time matched to imaging frame rate (resampled) during run
%epochs 
Behavior.resampled.run_time=runtime;
%Behavior.resampled.no_run_time=noruntime;

%resampled behavioral time to match imaging rate
Behavior.resampled.time=res_time;

%resampled behavioral lap index to match imaging rate
Behavior.resampled.lapNb=res_lapNb;

%cumulative position (absolute since beginning to recording session)
Behavior.resampled.cumulativeposition=res_cum_position;

%relative position to start of each lap (not normalized to each lap) -all
%restricted frames
Behavior.resampled.position=res_position;

%position during run epochs
Behavior.resampled.run_position=run_position;

%running epoch intervals in input frame space
Behavior.resampled.run_idx=run_int_idx;

%lap idx during running
Behavior.resampled.run_lapNb=run_lapNb;

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
    c2plot=options.c2plot;
    
    %choose random ROI
    %c2plot= randi(size(C_df,2));
    
    figure;  
    %run epochs
    subplot(2,1,1);
    hold on;

    onset_idx =find(onset_binary(:,c2plot) == 1);
    
    stem(Cdf_time(onset_idx), 2*onset_binary(onset_idx,c2plot),'bo');

    run_onset_idx =find(run_onset_binary(:,c2plot) == 1);

    scatter(Cdf_time(run_onset_idx), 2*run_onset_binary(run_onset_idx,c2plot),'*','MarkerEdgeColor', 'r');
    plot(Cdf_time, C_df(:,c2plot),'k');
    hold off
    
    %no run epochs
    subplot(2,1,2);
    hold on;
    
    stem(Cdf_time(onset_idx), 2*onset_binary(onset_idx,c2plot),'bo');
    
    norun_onset_idx =find(norun_onset_binary(:,c2plot) == 1);
    
    scatter(Cdf_time(norun_onset_idx), 2*norun_onset_binary(norun_onset_idx,c2plot),'*','MarkerEdgeColor', 'r');
    plot(Cdf_time, C_df(:,c2plot),'k');
    hold off
    
end

end



