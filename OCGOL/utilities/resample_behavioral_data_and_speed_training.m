function [Behavior] = resample_behavioral_data_and_speed_training(Behavior, Imaging, options)
% Resample behavioral data to match imaging acqusition frequency

%width of moving mean filter for Cossart (also filter for speed)
mov_wind = options.moving_window;

%% Select restricted vs. non-restricted data

%Imaging and behavior
if options.restrict==true
    
    %imaging traces
    %C_df = Imaging.trace_restricted;
    %imaging time (restricted)
    Cdf_time = Imaging.time_restricted;
    %behavior time (restricted)
    time = Behavior.restricted.time;
    %cumulative position (restricted)
    cum_position = Behavior.restricted.cumulativeposition;
    %normalized positions (0-1 each lap)
    norm_position = Behavior.restricted.normalizedposition;
    %relative position to start of lap (cm)
    position = Behavior.restricted.position;
    %index of lap at each behavior index
    lapNb = Behavior.restricted.lapNb;
    
elseif options.restrict==false
    %imaging traces
    %C_df = Imaging.trace;
    %imaging time
    Cdf_time = Imaging.time;
    %behavioral time
    time = Behavior.time;
    %cumulative position since start of recording
    cum_position = Behavior.cumulativeposition;
    %normalized position
    norm_position = Behavior.normalizedposition;
    %position (cm) relative to start of each laps
    position = Behavior.position;
    %get non-restricted version of this!!
    %lapNb = Behavior.restricted.lapNb;
end

%put behavioral time and position into 1 matrix together
time_position=[time position];

%% Resample behavioral data to match imaging frequency
%resample - change resolution while keeping the size the same

%take all the imaging times and see how many counts fall within each
%behavioral 
%bin will be the same size as Cdf_time - represents ths bin indices

%N is bin counts - how many events fall between each edges of time
%bin take each imaging time value and get corresponding bin that it
%corresponds to in the time (behavioral) vector
%histc return N that is the same size as time (edges)
%idea is that the imaging time will only fall within 1 bin of the much
%higher sampled signal (behavior)
%tells which in which behavioral bin (timepoint) each imaging timepoint
%belongs
%time (behavioral) - edges of bins
%use discretize since histc is no longer recommended
%[N,bin2]=histc(Cdf_time,time);

%this yields the same result as histc
%edges = time (input)
%bin - indices of bins in time (behavior - time range) that contain a
%given imaging time point in
[bin,edges] = discretize(Cdf_time,time);

%If missing behavior time - should never run, but valid check
%should have an assignment to a behavioral timepoint
%for all the imaging timepoints

%find bins that were not assigned (no corresponding behavioral recording to
%imaging recording)
notAssignedBins = find(bin == 0);

%give warning
if ~isempty(notAssignedBins)
    disp('Missing behavior recording data !!!')
end

%% Check time difference (s) between imaging times and resampled behav times
 
%imaging time - resampled beahvior time
time_diff_resampled = abs(Cdf_time - time(bin));
%shift behavior bin by 1 forward
time_diff_shift_bin_plus = abs(Cdf_time - time(bin+1));
%shift behavior bin by 1 backward
time_diff_shift_bin_minus = abs(Cdf_time - time(bin-1));

%mean diffs
mean_time_diffs(1) = mean(time_diff_resampled);
mean_time_diffs(2) = mean(time_diff_shift_bin_plus);
mean_time_diffs(3) = mean(time_diff_shift_bin_minus);
%min diffs
min_time_diffs(1) = min(time_diff_resampled);
min_time_diffs(2) = min(time_diff_shift_bin_plus);
min_time_diffs(3) = min(time_diff_shift_bin_minus);
%max diffs
max_time_diffs(1) = max(time_diff_resampled);
max_time_diffs(2) = max(time_diff_shift_bin_plus);
max_time_diffs(3) = max(time_diff_shift_bin_minus);

%print maximum difference during resampling
fprintf('The maximum time difference after resampling is %f seconds\n', max_time_diffs(1));


%% Resample the behavior data based on bin assignment

%combined time and position (cm since start)
res_time_position = time_position(bin,:);

%cumulative position since start of recording (cm)
res_cum_position = cum_position(bin);

%normalized position
res_norm_position = norm_position(bin);

%resample position from start of laps (cm)
res_position = position(bin);

%resampled lap idx - complete lap idx assigned to each timepoint
res_lapNb = lapNb(bin);

%resample time
res_time=time(bin);

%% Speed of animal and 

%insert dt value here (seconds)
Imaging.dt = 0.033;

% Mean measured framerate of imaging (Hz - frames/second)
avg_fr = 1/Imaging.dt;

%Speed
speed=[0;diff(res_cum_position)]*avg_fr;

%average speed using mean filter
speed = movmean(speed,[mov_wind]);



%% Save to Behavior struct

%resampled behavioral time to match imaging rate
Behavior.resampled.time = res_time;

%resampled behavioral lap index to match imaging rate
Behavior.resampled.lapNb = res_lapNb;

%cumulative position (absolute since beginning to recording session)
Behavior.resampled.cumulativeposition = res_cum_position;

%normalized position
Behavior.resampled.normalizedposition = res_norm_position;

%position (cm) from lap start
Behavior.resampled.position = res_position;

%spatial binning - import from spatial info script
%Behavior.spatialBins = 

%speed
Behavior.speed = speed;

end

