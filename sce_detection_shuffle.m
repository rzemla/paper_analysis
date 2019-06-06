function [outputArg1,outputArg2] = sce_detection_shuffle(noRun_thres_traces,run_binary_interval,thres_traces)

%NOT WORKING RIGHT

%% Shuffle detected events within run intervals
%take detected events (binary in event space before delay filtering)
%randomly scramble across the entire no run interval 

%indices where the animal is not running
noRun_int_idx = find(run_binary_interval == 0);
%indices where it is running
run_int_idx = find(run_binary_interval == 1);

%events only in no run interval that were selected
noRun_thres_traces;

%expand noRun idxs for all ROIs
noRun_idx_allROI = repmat(noRun_int_idx,1,size(noRun_thres_traces,2))';

%preallocate cells for N shuffles
shuffle_nb = 100;
shuffled_cell = cell(1,shuffle_nb);

%preallocate each cell with 0 matrix to speed up
for cc=1:shuffle_nb
    shuffled_cell_idx{cc} = zeros(size(noRun_thres_traces,2),size(noRun_int_idx,1));
end



% tic; %70 seconds for 100 shuffles
% for cc=1:shuffle_nb
%     %shuffle events in no run inverval - each ROI
%     for rr=1:size(noRun_thres_traces,2)
%         shuffled_cell{cc}(rr,:) = randperm(size(noRun_int_idx,1));
%     end
%     disp(cc);
% end
% toc;

tic; %23 seconds (>3 fold speed-up compared to above code)
for cc=1:shuffle_nb
    [~, shuffled_cell_idx{cc}] = sort(rand(size(noRun_thres_traces,2),size(noRun_int_idx,1)),2);
    disp(cc);
end
toc;

%preallocate
for cc=1:shuffle_nb
    shuffled_cell_frame{cc} = noRun_idx_allROI;
end

%for each shuffle
for cc=1:shuffle_nb
    %for each ROI
    for rr=1:size(noRun_thres_traces,2)
        shuffled_cell_frame{cc}(rr,:) = noRun_idx_allROI(rr,(shuffled_cell_idx{cc}(rr,:)));
    end
    disp(cc);
end

%reconsitute a shuffled event matrix - use algorithm below that is faster
%preallocate using copy of original matrix
% shuffled_events = noRun_thres_traces';
% for rr=1:size(noRun_thres_traces,2)
%     shuffled_events(rr,noRun_int_idx) = noRun_thres_traces(shuffle_idx(rr,:),rr);
% end

%preallocate cells for shuffle events
shuffle_nb = 100;
shuffled_events = cell(1,shuffle_nb);

%assign the existing matrix to the matrix holders
for cc=1:shuffle_nb
    shuffled_events{cc} = noRun_thres_traces';
end

%permute the no run epochs based on the shuffle above
for cc=1:shuffle_nb
    shuffled_events{cc}(:,noRun_int_idx) = noRun_thres_traces(shuffled_cell_frame{cc});
    %shuffled_events{cc}(:,noRun_int_idx) = shuffled_events{cc}(shuffled_cell{cc});
end

%transpose the matrices back to the row,column arrange of inputs
for cc=1:shuffle_nb
    shuffled_events{cc} = shuffled_events{cc}';
end


%compare original vs shuffled events matrix in subplot
figure;
subplot(2,1,1)
imagesc(noRun_thres_traces')

subplot(2,1,2)
imagesc(shuffled_events{cc}')


%%%% DETECT SCEs BELOW AND GET THRESHOLD %%%%

%original matrix
original_input = noRun_thres_traces;


noRun_thres_traces = shuffled_events{10};

%noRun_thres_traces = original_input;

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

%plot
figure;
hold on
ylabel('Synchronous event count')
plot(sce_event_count)
plot([1 size(thres_traces,1)],[min_cell_nb min_cell_nb],'r--')

end

