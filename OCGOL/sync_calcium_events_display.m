%% Load in animal of interest

options.register = 0;

%lab workstation
%input directories to matching function
%path_dir = {'G:\OCGOL_training\I56_RLTS_041019\5A5B'};
path_dir = {'G:\OCGOL_training\I56_RLTS_041019\ABrand_no_punish_041619'};
%cross session directory
%crossdir = 'G:\OCGOL_training\I56_RLTS_041019\crossSession';

%load place cell variables for each session
%get mat directories in each output folder
for ii=1:size(path_dir,2)
    %get matfile names for each session
    matfiles{ii} = dir([path_dir{ii},'\output','\*.mat']);
end
%load in place cell variables (and others later)
for ii = 1:size(path_dir,2)
    %add event variables
    session_vars{ii} = load(fullfile(matfiles{ii}.folder,matfiles{ii}.name),'Place_cell', 'Behavior',...
        'Behavior_split_lap','Behavior_split','Events_split','Events_split_lap',...
        'Imaging_split', 'Imaging');
end

%% %Plot velcity, time, and lap indices

speed = session_vars{1}.Behavior.speed;
time = session_vars{1}.Imaging.time_restricted;
run_onsets = session_vars{1}.Events_split{3}.Run;
norun_onsets = session_vars{1}.Events_split{3}.NoRun;
norm_position = session_vars{1}.Behavior.resampled.normalizedposition;
run_epoch_binary = session_vars{1}.Behavior.run_ones;
traces = session_vars{1}.Imaging.trace_restricted;

%SI score
A_tuned = session_vars{1, 1}.Place_cell{1, 4}.Spatial_Info.significant_ROI;
B_tuned = session_vars{1, 1}.Place_cell{1, 5}.Spatial_Info.significant_ROI;
AorB_tuned = A_tuned | B_tuned;

%binary onset logicals
run_binary = run_onsets.run_onset_binary;
norun_binary = norun_onsets.norun_onset_binary;
%combine binary onsets:
combined_binary = run_binary + norun_binary;

%total (restricted) frames: 26023

%sum all the events from all cells
summed_events_all = sum(combined_binary,2);

%run moving sum window over 6 frames (~200 ms window) +1 (~230 ms = 7 frames total))
%5 frames = 170 ms
windowed_events = movsum(summed_events_all,[3, 3]);

%take neurons tuned to A
tunedA_binary = combined_binary(:,AorB_tuned);
%sum events across both run and no-run epochs
windowA_tuned =  movsum(sum(tunedA_binary,2),[3, 3]);

figure;
for rr = 20
    subplot(3,1,1)
    hold on;
    %speed
    plot(time,speed,'k');
    ylabel('Speed [cm/s]');
    xlabel('Time [s]');
    subplot(3,1,2)
    %position (norm)
    hold on
    ylabel('Normalized Position');
    plot(time,norm_position,'k');

    subplot(3,1,3)
    hold on
    ylabel('Summed event count');
    %sample ROI during run 
    plot(time,windowA_tuned,'k');
    
    %overlay no-run binary
    plot(time,10*(~run_epoch_binary),'r')
    
    %plot(time,run_binary(:,rr), 'g');
    %plot(time,norun_binary(:,rr), 'r');
    hold off

end
%% Plot summed events in separate figure

figure;
hold on;
ylabel('Summed event count');
%sample ROI during run
plot(time,windowA_tuned,'k');
%overlay no-run binary
plot(time,10*(~run_epoch_binary),'r')

%% Get indices of max sync activity

%find indices where summed events exceed 5 in no-run epochs
windowA_tuned_noRun = windowA_tuned;
%zero out run epochs
windowA_tuned_noRun(logical(run_epoch_binary)) = 0;

%more than threshold of summed events
thres_idx_restricted = find(windowA_tuned_noRun >= 8);
%time of restricted indices
thres_idx_times  = time(thres_idx_restricted);

%work on selecting only no run epochs during B runs


%time of entire stack
global_time = session_vars{1}.Imaging.time_restricted;

%find idx of putative SCAs in stack
% for ii=1:size(thres_idx_times,1)
% frames_idxs(ii) = find(thres_idx_times(ii))

[~,idx_global,~] = intersect(global_time,thres_idx_times,'stable');

%translated to abs frame count

%% Plot run event deleted version
figure;
hold on;
ylabel('Summed event count');
%sample ROI during run
plot(time,windowA_tuned_noRun,'k');
%overlay no-run binary
plot(time,10*(~run_epoch_binary),'r')
%% On frame domain
figure;
hold on;
ylabel('Summed event count');
%sample ROI during run
plot(windowA_tuned_noRun,'k');
%overlay no-run binary
plot(10*(~run_epoch_binary),'r')

%% Look at 3500, 5800, 6900, 33175 frame in global reference

%find neurons with onset aroudn 33175 - 33000 - 33200

if 0
    
    %all in no run period
    onsets_arnd_SCA = combined_binary(st_idx:end_idx,:) >= 1;
    
    run_binary_segment = run_epoch_binary(st_idx:end_idx);
    %take out no run segment
    onsets_arnd_SCA_noRun  = onsets_arnd_SCA(~run_binary_segment,:);
    %keep neurons with events
    onset_ROIs_log = sum(onsets_arnd_SCA_noRun)>1;
    
    onset_ROIs_idx = find(onset_ROIs_log ==1);
    
    [~, I] = max(onsets_arnd_SCA,[],1)
end

%% Sort neurons by max dF/F in given region

% onset_ROIs_log (override with all tuned ROI to both or either trial)
onset_ROIs_log = AorB_tuned;

%6238 - start point for preplay!
%11511 - nice replay of B!
%19175 - nice preplay in B!

%start and end points of sorted
st_evt_sort = 19175;
end_evt_sort = st_evt_sort+10;
plot_range = [1000, 2000];

%start idx (absolute)
st_idx =st_evt_sort-plot_range(1);
%end idx (absolute)
end_idx = st_evt_sort+plot_range(2);

%input_dFF_matrix = traces(2360+700:2900+700,onset_ROIs_log)';

%get absolute idx's of neurons
input_neuron_idxs = find(onset_ROIs_log ==1);

%1500-2000 +700
input_events_matrix = combined_binary(st_evt_sort:end_evt_sort,onset_ROIs_log)';

%sort by event onset
%for each ROI
for rr=1:size(input_events_matrix,1)
    if  ~isempty(find(input_events_matrix(rr,:) > 0,1))
        loc_ROI(rr) = find(input_events_matrix(rr,:) > 0,1)
    else
        loc_ROI(rr) = 0;
    end
end

%sort by index
[M_sort,I] =sort(loc_ROI,'ascend');


%remove neurons that do not have an event in region of interest
I(M_sort == 0) = [];

%neuron idx (global) sorted
input_neuron_idxs_sorted = input_neuron_idxs(I);

figure;
imagesc(input_events_matrix(I,:))
hold on

%% Plot speed, position and dF/F trace of neurons prior to SCE in no-run epoch

%start_idx

figure
subplot(4,1,1)
hold on;
xlim([1,end_idx-st_idx])
%speed
plot(speed(st_idx:end_idx),'k');
%plot(time(st_idx:end_idx),speed(st_idx:end_idx),'k');
ylabel('Speed [cm/s]');
xlabel('Time [s]');

subplot(4,1,2)
%position (norm)
hold on
xlim([1,end_idx-st_idx])
%xlim([st_idx,end_idx])
ylabel('Normalized Position');
plot(norm_position(st_idx:end_idx),'k');
%plot(time(st_idx:end_idx),norm_position(st_idx:end_idx),'k');
hold off

%dF/F
subplot(4,1,[3 4])
imagesc(traces(st_idx:end_idx,input_neuron_idxs_sorted)')
hold on
axis normal;
caxis([0 1]);
ylabel('Neuron #');
colormap(gca,'jet')
hold off

%% Plot traces as line plot (vs imagesc)
figure;
hold on;
stepSize = 2;
step = 0;
for ii=1:size(input_neuron_idxs_sorted,2)
    plot(traces(st_idx:end_idx,input_neuron_idxs_sorted(ii))-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end

ylabel('Normalized Position');
plot(norm_position(st_idx:end_idx)-step,'r');

%% Plot all ROIs

st_idx_all = 1;
end_idx_all = 36000;
figure
subplot(4,1,1)
hold on;
xlim([1,end_idx_all-st_idx_all])
%speed
plot(speed(st_idx_all:end_idx_all),'k');
%plot(time(st_idx:end_idx),speed(st_idx:end_idx),'k');
ylabel('Speed [cm/s]');
xlabel('Time [s]');

subplot(4,1,2)
%position (norm)
hold on
xlim([1,end_idx_all-st_idx_all])
%xlim([st_idx,end_idx])
ylabel('Normalized Position');
plot(norm_position(st_idx_all:end_idx_all),'k');
%plot(time(st_idx:end_idx),norm_position(st_idx:end_idx),'k');
hold off

%dF/F
subplot(4,1,3)
imagesc(traces(st_idx_all:end_idx_all,AorB_tuned)')
hold on
axis normal;
caxis([0 1]);
ylabel('Neuron #');
colormap(gca,'jet')
hold off

%events
subplot(4,1,4)
imagesc(combined_binary(st_idx_all:end_idx_all,AorB_tuned)')
hold on
axis normal;
caxis([0 1])
colormap(gca,'gray')
caxis(gca,[0 1])
ylabel('Neuron #');


%events
% subplot(4,1,4)
% imagesc(combined_binary(st_idx:end_idx,input_neuron_idxs_sorted)')
% hold on
% axis normal;
% caxis([0 1])
% colormap(gca,'gray')
% caxis(gca,[0 1])
% ylabel('Neuron #');


%active in no run
%31,38, 68

%sorted_ROI_list = session_vars{1, 1}.Place_cell{1, 3}.ROI_Spatial_tuning_sorted;

% 
% flash at 892/893; 7862, 
% 
% global
% flash at 895, 2145 (sustained), 32035, 33175
