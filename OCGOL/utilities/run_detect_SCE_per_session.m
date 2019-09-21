function [SCE] = run_detect_SCE_per_session(imaging_traces, position_norm,run_epoch)
%Detect SCEs for each session with given type of trials


%% 3rd order Savitzky-Golay filter frame size 500ms (~15 frames)

%order of filter
sg_order = 3;
%window of filter
sg_window = 15;
%filter along first dimension (across rows)
filtered_traces = sgolayfilt(imaging_traces,sg_order,sg_window,[],1);

%plot side by side (sample ROI
ROI = 1;
figure;
hold on
plot(imaging_traces(:,ROI), 'k')
plot(filtered_traces(:,ROI)-1,'r');

%% Plot traces as line plot (vs imagesc)
figure;
subplot(1,2,1)
hold on;
%ylim([-0.5 1])
title('Non-filtered calcium traces')
stepSize = 2;
step = 0;
for ii=1:30%size(filtered_traces,2)
    plot(imaging_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end
subplot(1,2,2)
hold on;
%ylim([-0.5 1])
title('Savitzky-Golay filtered calcium traces')
stepSize = 2;
step = 0;
%first 30 ROIs
for ii=1:30size(filtered_traces,2)
    plot(filtered_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end

%ylabel('Normalized Position');
%plot(norm_position(st_idx:end_idx)-step,'r');

%% Detect events based on threshold
%for each cell
%sum of the median value with 3x interquartile range calculated within:
%sliding window -2/+2 s 
%2s = 60 frames
win_width = 60;

%moving median
med_traces = movmedian(filtered_traces,[win_width win_width],1);

%create blank iqr vector for 1 neuron (for all neurons - ROIx time)
iqr_range = zeros(size(filtered_traces,2),size(filtered_traces,1));
%interquartile range (start at first index will be win_width+1 etc
for ii=1:size(filtered_traces,1)-2*win_width %(set range here and add edges detection in the future
    iqr_range(:,win_width+ii) = iqr(filtered_traces(ii:((2*win_width)+ii),:),1);
end

%set threshold for detection later on
event_thres = 3*iqr_range + med_traces';
%set to time x ROI format
event_thres = event_thres';

%% Plot trace, median and irq range - works
figure
hold on
stepSize = 2;
step = 0;
for ii =1:30
    ROI=ii;
    %plot all
    %trace
    plot(filtered_traces(:,ROI) - step,'k');
    %sliding median
    plot(med_traces(:,ROI) - step, 'r');
    
    %sliding interquartile range
    %plot(iqr_range, 'g');
    plot((3*iqr_range(ROI,:) + med_traces(:,ROI)') - step, 'b');
    step = step - stepSize;
end

%position
plot(position_norm-step+2,'r');

%% Detect events above 3x IQR interval

%get logical of all points where event exceeds threshold
thres_traces = filtered_traces > event_thres;

figure;
imagesc(thres_traces')

%% Artifact correction - set 60 frames early and 60 frames late to 0 in thres_traces 
%edge artifact from using filtering and event detection
%use win_width size to clear artifacts from the edges

%save_original for comparison purposes
thres_traces_orig = thres_traces;

%set start frames of filtered events
thres_traces([1:win_width],:) = 0;
thres_traces([(size(thres_traces,1)-win_width+1):size(thres_traces,1)],:) = 0;


%plot to see that there are no obvious artifacts
figure
subplot(1,2,1)
imagesc(thres_traces_orig')
hold on
title('Original')
subplot(1,2,2)
imagesc(thres_traces')
hold on
title('Edge artifact cleared');

%% Select only events in noRun intervals

%select binary interval for select processing interval
run_binary_interval = run_epoch;

%remove run events
noRun_thres_traces = thres_traces & ~repmat(logical(run_binary_interval),1,size(thres_traces,2));

figure;
subplot(3,1,1)
imagesc(noRun_thres_traces')
subplot(3,1,2)
imagesc(~repmat(logical(run_binary_interval),1,size(thres_traces,2))')
subplot(3,1,3)
hold on
title('No run interval')
ylim([0 2]);
xlim([1 size(thres_traces,1)]);
plot(~run_binary_interval,'r')

%% Create an onset matrix (not filtered by delay separation)

diff_thres = diff(double(noRun_thres_traces),1,1);

%only get onsets
diff_thres = diff_thres == 1; 

%get index of each onset and check if event is within 
%for each ROI
for rr = 1:size(diff_thres,2)
    event_idx{rr} = find(diff_thres(:,rr) == 1);
end

%reconstruct matrix with onsets alone
%create blank 
noRun_event_matrix = zeros(size(diff_thres,1),size(diff_thres,2));
%reconstruct filtered events for each ROI
for rr=1:size(event_idx,2)
    noRun_event_matrix(event_idx{rr},rr) = 1;
end

%add zero to first frame to account for offset
noRun_event_matrix = [zeros(1,size(noRun_event_matrix,2));noRun_event_matrix];

%plot
figure;
imagesc(noRun_event_matrix')

%% Run SCE shuffle - use event matrix above as input
%make option to set # of shuffles
%also return mean, std, and 
[sce_threshold] = sce_detection_shuffle(noRun_event_matrix,run_binary_interval,thres_traces);
%previous thres = 8.33 sync events on frame

%% Select individual onsets separated at least 1s - WORKS!

%frames - 1s = 30 frames
min_dur = 30;

%get onsets and offsets
diff_thres = diff(double(noRun_thres_traces),1,1);

%only get onsets
diff_thres = diff_thres == 1; 

%get index of each onset and check if event is within 
%for each ROI
for rr = 1:size(diff_thres,2)
    event_idx{rr} = find(diff_thres(:,rr) == 1);
end

%run diff, if diff < min_dur, remove event
%do this iteratively for each ROI
%set iterative flag
check_events = 1;

for  rr = 1:size(diff_thres,2)
    diff_events{rr} = diff(event_idx{rr});
    
    while check_events == 1
        
        %find first diff less than frame space duration, remove recalulate diff
        temp_dur_flag = find(diff_events{rr} < min_dur,1);
        
        %if dur_flag not empty
        if ~isempty(temp_dur_flag)
            %remove that +1 event
            event_idx{rr}(temp_dur_flag+1) = [];
            %recalulate diff on updated event list
            diff_events{rr} = diff(event_idx{rr});
            %keep flag on
            check_events = 1;
        else
            %reset flag
            check_events = 0;
        end
    end
    %reset flag for next ROI
    check_events = 1;
end

%reconstruct events
%create blank 
dur_filtered_event = zeros(size(diff_thres,1),size(diff_thres,2));
%reconstruct filtered events for each ROI
for rr=1:size(event_idx,2)
    dur_filtered_event(event_idx{rr},rr) = 1;
end

%filter out events with in 1s of one another - lots of code when  more than
%1 events that are clustered events
%sum_dur_mat = movsum(diff_thres, [0 min_dur-1],1);

%% Plot final filter
figure;
subplot(2,1,1)
hold on
xlim([1 size(thres_traces,1)])
stepSize = 2;
step = 0;
for ii=1:30%size(filtered_traces,2)
    plot(noRun_thres_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end
%check onset and duration separation filter
subplot(2,1,2)
hold on
xlim([1 size(thres_traces,1)])
stepSize = 2;
step = 0;
for ii=1:30%size(filtered_traces,2)
    plot(diff_thres(:,ii)-step, 'k', 'LineWidth', 1.5)
    plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
    step = step - stepSize;
end

%% Separate plot to check sync events
figure
hold on
xlim([1 size(thres_traces,1)])
stepSize = 2;
step = 0;
for ii=1:100%size(filtered_traces,2) %first 100 ROIs to check
    plot(diff_thres(:,ii)-step, 'k', 'LineWidth', 1.5)
    plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
    step = step - stepSize;
end

%% Pad one index with 0 that's lost with derivative processing above

final_events = [zeros(1,size(dur_filtered_event,2));dur_filtered_event];

%% Detect SCE by running moving sum across the filtered traces and 

%minimum cells that must particupate in SCE
min_cell_nb = 5;

%width of sync events - 200ms = 6 frames
%5 will give 2 behind, center and 2 ahead
%6 will gives 3 beind, cetner
sync_window_width = 6;

%collapse all the no run SCE events 
summed_events_SCE = sum(final_events,2);

%count sync events within time window witdh
sce_event_count = movsum(summed_events_SCE,sync_window_width);

%make copy and save as original
sce_event_count_orig = sce_event_count;

%plot
figure;
hold on
ylabel('Synchronous event count')
ylim([0 25])
plot(sce_event_count)
%minimum cell involvement line
plot([1 size(thres_traces,1)],[min_cell_nb min_cell_nb],'k--')
%threshold determined by shuffle
plot([1 size(thres_traces,1)],[sce_threshold sce_threshold],'r--')

%% Find neurons involved in each SCE
%within 200 ms  = 6 frames
%500 ms = 15 frames

%get of SCE frame indices; use 5 as a cutoff for now - use shuffle-based
%threshold later
%threshold or min number - whichever is greater
if sce_threshold > 5
    sync_idx = find(sce_event_count >= sce_threshold);
    %number of calcium events associated with each SCE
    sync_event_count = sce_event_count(sync_idx);
else
    sync_idx = find(sce_event_count >= 5);
    sync_event_count = sce_event_count(sync_idx);
end

%combine sync_idx and sync_event_count into 1 matrix
sync_comb = [sync_idx sync_event_count];

%sync endpoints
sync_end = find(diff(sync_idx)>1);
%add 1 as starting index
sync_start = [1; sync_end+1];
%end end index to sync end
sync_end = [sync_end; size(sync_idx,1)];
%combined sync_start and sync_end
sync_range = [sync_start, sync_end];

%size of sync range equal to number of SCEs
nbSCE = size(sync_range,1);


%use all the indices for now, reconsider later (i.e. exclude some/all of
%the neighboring intervals (i.e. filter here in the future)
%proposal: for each segment of frame neighboring ROIs, select 1 index which
%has the most ROIs involved in SCE

%for each SCE, identify the neurons involved (based on window width of 6
%frames)
%defined above

%for each index, take frame range (3 idx before until 2 idx after)
%take frame slide of processed activity, sum and include ROI idx
for ss=1:size(sync_idx,1)
    SCE_ROIs{ss} = find(sum(final_events(sync_idx(ss)-3:sync_idx(ss)+2,:),1) > 0);
end

%numbers of ROIs participating in each SCE
size_dim2 = @(x) size(x,2);
SCE_ROI_nb = cell2mat(cellfun(size_dim2,SCE_ROIs,'UniformOutput',false));

%add to combined matrix
sync_comb = [sync_comb, SCE_ROI_nb'];

%% Output_struct
SCE.sync_comb = sync_comb;
SCE.nbSCE = nbSCE;
%SCE.final_events_shuffle = final_events_shuffle;
SCE.SCE_ROIs = SCE_ROIs;
SCE.sce_threshold = sce_threshold;
SCE.sync_range = sync_range;
SCE.sync_idx = sync_idx;
SCE.sce_event_count = sce_event_count;
SCE.sce_event_count_orig = sce_event_count_orig;
SCE.summed_events_SCE = summed_events_SCE;

end

