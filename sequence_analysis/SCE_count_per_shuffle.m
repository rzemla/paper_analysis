function [sce_event_count,final_events] = SCE_count_per_shuffle(noRun_thres_traces)
%function description
%input: event onset matrix permute during no run epochs
%output: return SCE event count


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
% figure;
% subplot(2,1,1)
% hold on
% xlim([1 size(thres_traces,1)])
% stepSize = 2;
% step = 0;
% for ii=1:30%size(filtered_traces,2)
%     plot(noRun_thres_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
%     step = step - stepSize;
% end
% %check onset and duration separation filter
% subplot(2,1,2)
% hold on
% xlim([1 size(thres_traces,1)])
% stepSize = 2;
% step = 0;
% for ii=1:30%size(filtered_traces,2)
%     plot(diff_thres(:,ii)-step, 'k', 'LineWidth', 1.5)
%     plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
%     step = step - stepSize;
% end

%% Separate plot to check sync events
% figure
% hold on
% xlim([1 size(thres_traces,1)])
% stepSize = 2;
% step = 0;
% for ii=1:100%size(filtered_traces,2) %first 100 ROIs to check
%     plot(diff_thres(:,ii)-step, 'k', 'LineWidth', 1.5)
%     plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
%     step = step - stepSize;
% end

%% Pad one index with 0 that's lost with derivative processing above

final_events = [zeros(1,size(dur_filtered_event,2));dur_filtered_event];

%% Detect SCE by running moving sum across the filtered traces and 

%minimum cells that must particupate in SCE
min_cell_nb = 5;

%width of sync events - 200ms = 6 frames
%5 will give 2 behind, center and 2 ahead
%6 will gives 3 beind, cetner - verified
sync_window_width = 6;

%collapse all the no run SCE events 
summed_events_SCE = sum(final_events,2);

%count sync events within time window witdh
sce_event_count = movsum(summed_events_SCE,sync_window_width);

%plot
if 0
figure;
hold on
ylim([0 40])
ylabel('Synchronous event count')
plot(sce_event_count(1:1000))
%5 cell threshold
plot([1 size(final_events,1)],[min_cell_nb min_cell_nb],'r--')
end

end

