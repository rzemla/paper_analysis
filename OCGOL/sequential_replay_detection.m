%% Sequential replay detection

%% Import/define variables

%calcium traces
traces =session_vars{1}.Imaging_split{4}.trace_restricted;

% %ROIs associated with SCEs
% SCE_ROIs 
% %indices of SCEs
% sync_idx

sce_nb=6;
%sync_idx(sce_nb) %--> nice replay
%SCE_ROIs{1}

%convert ROIs to logical
SCE_ROI_logi = false(1,size(final_events,2));
SCE_ROI_logi(SCE_ROIs{sce_nb}) = true; 


%% Extract calcium transient of each involved cell over 2s time windows
%2s window = 60 frames (30 before and 29 after (involved center point)
fr_range = 60;

%for each SCE
for ss=1:size(sync_idx,1)
    %use frames of sce as center for calcium traces
    SCE_traces{ss} = traces(sync_idx(ss)-(fr_range/2):sync_idx(ss)+((fr_range/2) - 1),SCE_ROIs{ss});
    %calculate reference as the median transient among cells involved
    ref_transient{ss} = median(SCE_traces{ss},2);
end


%% Plot involved traces as line plot 
figure;
%first 5 for show
for ss=1:5%size(sync_idx,1)
 subplot(2,1,1)
hold on
title('All transients involved in SCE')
xlim([1 fr_range])
stepSize = 2;
step = 0;   
    for ii=1:size(SCE_traces{ss},2)
        plot(SCE_traces{ss}(:,ii), 'LineWidth', 1.5)
        %plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
        %step = step - stepSize;
    end
    
    
    subplot(2,1,2)
    hold on
    title('Reference transient vs first transient')
    xlim([1 fr_range])
    %plot reference
    plot(ref_transient{ss},'k');
    %plot first transient in series
    plot(SCE_traces{ss}(:,1),'b')
    
    %pause;
    %clf
end

%% Calculate normalized covariance (correlation) for delay between -200ms and 200 ms) (cross-correlation)
%cross covariance lag in frames
fr_delay =6;

%time diff between frames (seconds)
dt = 0.0334;
time_lag = [0:dt:1.03];
time_lag = [-1*fliplr(time_lag),time_lag(2:end)];

%insert a 100x expanded timebase for max extrapolation from parbola fit
expanded_time =[0:dt/100:1.03];
expanded_time = [-1*fliplr(expanded_time),expanded_time(2:end)];


%for each SCE
for ss=1:size(sync_idx,1)
    %for each SCE transient signal
    for rr=1:size(SCE_traces{ss},2)
        %calculate normalized covariance for each trace against ref transient
        %set normalized biased or unbiased
        %xcorr
        %[r{ss}(rr,:) ,lags{ss}(rr,:)] = xcov(ref_transient{ss},SCE_traces{ss}(:,rr),fr_delay,'biased');
        %this function normalizes it as one would for a simple correlation
        %(may write custom in future to avoid using toolbox fxn)
        [r{ss}(rr,:) ,lags{ss}(rr,:),~] = crosscorr(ref_transient{ss},SCE_traces{ss}(:,rr),fr_delay);
    end
end

%get max value of covariance for each ROI in SCE and discard those with
%lower max than 0.6
%check to see if any %62 ROIs to excludce on first pass
for ss=1:size(sync_idx,1)
    %find ROIs with corr less than 0.6
    max_sync{ss} = find(max(r{ss},[],2) < 0.6)';
end
%number of neurons that were removed (considered noisy by this metric)
cell2mat(max_sync);

%manual way to get normalized covariance
%r{1, 6}(1,:)/std(ref_transient{ss})*std(SCE_traces{ss}(:,rr))

%for each SCE
for ss=1:size(sync_idx,1)
    %for each SCE transient signal
    for rr=1:size(SCE_traces{ss},2)
        %fit parabola to (2nd order polynomial to covariance within lag range)
        %insert x in second time domain
        p{ss}(rr,:) = polyfit(time_lag(25:37),r{ss}(rr,:),2);
    end
end

%generate parabola
%x1 = time_lag(25:37);

%extrapolate in expanded time domain
x1 = expanded_time;

for ss=1:size(sync_idx,1)
    %for each SCE transient signal
    for rr=1:size(SCE_traces{ss},2)
        y1{ss}(rr,:) = polyval(p{ss}(rr,:),x1);
    end
end

for ss=1:size(sync_idx,1)
    %for each SCE transient signal
    for rr=1:size(SCE_traces{ss},2)
        
        %find max time for each ROI in SCE
        [~,idx_max{ss}(rr)] = max(y1{ss}(rr,:));
        max_onset{ss}(rr) = x1(idx_max{ss}(rr));
    end
end

%plot
% figure
% hold on
% %plot original fit range
% plot(time_lag(25:37),r,'k')
% %time extrapolated fit
% plot(x1,y1,'r')
% %maximum point
% scatter(max_onset, y1(idx_max),'b');


%% Plot example with RUN sequence points and order by SCE onset
%start frames of id'd SCE
sync_idx;
%onsets of each SCE event
max_onset;
%ROIs identified in SCE
SCE_ROIs;
%neurons that are part of SCE and involved in RUN sequence
neurons_participating{6};
%number of the SCE
%nice preplay example: 169
%nice replay example: 262
sce_nb = 83;

%temporary generate times and speed for select input inverval
time_choice = session_vars{1}.Behavior_split{4}.resampled.time;
time = session_vars{1, 1}.Behavior.resampled.time;

[~,select_speed_idx,~] = intersect(time,time_choice,'stable');
speed = session_vars{1, 1}.Behavior.speed;
%overwrite
speed = speed(select_speed_idx);

%start and end points of sorted (10-15 frame range) - 500 ms = 15
st_evt_sort = sync_idx(sce_nb);
% st_evt_sort = 11511;
end_evt_sort = st_evt_sort+15;

plot_range = [1200, 1500];

%start idx (absolute)
st_idx =st_evt_sort-plot_range(1);
%end idx (absolute)
end_idx = st_evt_sort+plot_range(2);

%get overlap between run sequence neurons and SCE neurons
[run_sce_neurons, run_sce_idx,~] = intersect(SCE_ROIs{sce_nb},neurons_participating{sce_nb},'stable');

%only ROIs in SCE that are part of run sequence
run_SCE_ROIs = SCE_ROIs{sce_nb}(run_sce_idx);

%only SCE onsets of ROIS that are part of run sequence
run_SCE_onsets = max_onset{sce_nb}(run_sce_idx);

%sort only run sequence involved neurons
[~,I_sce_run] = sort(max_onset{sce_nb}(run_sce_idx),'ascend');

%sort SCE inputs by onset time
[~,I_sce] = sort(max_onset{sce_nb},'ascend');


figure;
imagesc(traces(st_idx:end_idx,SCE_ROIs{sce_nb}(I_sce))')
hold on;
title('All SCE ROIS traces sorted by onset time')

%
% Plot speed, position and dF/F trace of neurons prior to SCE in no-run epoch
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
plot(position_norm(st_idx:end_idx),'k');
%plot(time(st_idx:end_idx),norm_position(st_idx:end_idx),'k');
hold off

%dF/F
subplot(4,1,[3 4])
imagesc(traces(st_idx:end_idx,SCE_ROIs{sce_nb}(I_sce))')
hold on
axis normal;
caxis([0 1]);
ylabel('Neuron #');
colormap(gca,'jet')
hold off

% only with RUN sequence neurons
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
plot(position_norm(st_idx:end_idx),'k');
%plot(time(st_idx:end_idx),norm_position(st_idx:end_idx),'k');
hold off

%dF/F
subplot(4,1,[3 4])
imagesc(traces(st_idx:end_idx,run_SCE_ROIs(I_sce_run))')
hold on
title('RUN sequence neurons only');
axis normal;
caxis([0 1]);
ylabel('Neuron #');
colormap(gca,'jet')
hold off

% Correlate SCE activation onset time with median (try mean)

%median normalized position onset across laps (run epochs)
median_run_seq_onset
%indices of RUN sequence neurons
recurring_neuron_idx
%neuron idxs of those involved in SCE and RUN sequence
run_SCE_ROIs

%relative onsets of neurons in SCEs (that are also in RUN sequence)
run_SCE_onsets

[~,~,recur_idx_pos] = intersect(run_SCE_ROIs,recurring_neuron_idx,'stable');

%correlate SCE onsets with median position of firing
[rho,p] =  corr(run_SCE_onsets', median_run_seq_onset(recur_idx_pos)','Type','Spearman')
%if p less than 0.05 and positive --> forward replay; 
%if spearman correlation negative --> reverse replay;




