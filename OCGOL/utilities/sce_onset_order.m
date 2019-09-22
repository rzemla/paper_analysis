function [outputArg1,outputArg2] = sce_onset_order(session_vars,SCE)
%% Sequential replay detection

%% Import/define variables

ss=4;

%calcium traces
traces =session_vars{ss}.Imaging.trace_restricted;

% %convert ROIs to logical
% SCE_ROI_logi = false(1,size(final_events,2));
% SCE_ROI_logi(SCE_ROIs{sce_nb}) = true; 
%sce_nb=6;

%% Extract calcium transient of each involved cell over 2s time windows
%2s window = 60 frames (30 before and 29 after (involved center point)
fr_range = 60;

%for each SCE
for cc=1:SCE{ss}.nbSCE
    %use frames of sce as center for calcium traces
    SCE_traces{cc} = traces(SCE{ss}.sync_idx(SCE{ss}.sync_range(cc,1))-(fr_range/2):...
                                    SCE{ss}.sync_idx(SCE{ss}.sync_range(cc,1))+((fr_range/2) - 1),...
                                    SCE{ss}.SCE_unique_ROIs{cc});
    %calculate reference as the median transient among cells involved
    ref_transient{cc} = median(SCE_traces{cc},2);
end


%% Plot involved traces as line plot 
figure;
%first 5 for show
for cc=1:5%size(sync_idx,1)
 subplot(2,1,1)
hold on
title('All transients involved in SCE')
xlim([1 fr_range])
stepSize = 2;
step = 0;   
    for ii=1:size(SCE_traces{cc},2)
        plot(SCE_traces{cc}(:,ii), 'LineWidth', 1.5)
        %plot(dur_filtered_event(:,ii)-step, 'r', 'LineWidth', 1.5)
        %step = step - stepSize;
    end
    
    
    subplot(2,1,2)
    hold on
    title('Reference transient vs first transient')
    xlim([1 fr_range])
    %plot reference
    plot(ref_transient{cc},'k');
    %plot first transient in series
    plot(SCE_traces{cc}(:,1),'b')
    
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
for cc=1:SCE{ss}.nbSCE  
    %for each SCE transient signal
    for rr=1:size(SCE_traces{cc},2)
        %calculate normalized covariance for each trace against ref transient
        %set normalized biased or unbiased
        %xcorr
        %[r{ss}(rr,:) ,lags{ss}(rr,:)] = xcov(ref_transient{ss},SCE_traces{ss}(:,rr),fr_delay,'biased');
        %this function normalizes it as one would for a simple correlation
        %(may write custom in future to avoid using toolbox fxn)
        [r{cc}(rr,:) ,lags{cc}(rr,:),~] = crosscorr(ref_transient{cc},SCE_traces{cc}(:,rr),fr_delay);
    end
end

%get max value of covariance for each ROI in SCE and discard those with
%lower max than 0.6
%check to see if any %62 ROIs to excludce on first pass
for cc=1:SCE{ss}.nbSCE  
    %find ROIs with corr less than 0.6
    max_sync{cc} = find(max(r{cc},[],2) < 0.6)';
end
%number of neurons that were removed (considered noisy by this metric)
cell2mat(max_sync);

%manual way to get normalized covariance
%r{1, 6}(1,:)/std(ref_transient{ss})*std(SCE_traces{ss}(:,rr))

%for each SCE
for cc=1:SCE{ss}.nbSCE  
    %for each SCE transient signal
    for rr=1:size(SCE_traces{cc},2)
        %fit parabola to (2nd order polynomial to covariance within lag range)
        %insert x in second time domain
        p{cc}(rr,:) = polyfit(time_lag(25:37),r{cc}(rr,:),2);
    end
end

%generate parabola
%x1 = time_lag(25:37);

%extrapolate in expanded time domain
x1 = expanded_time;

for cc=1:SCE{ss}.nbSCE  
    %for each SCE transient signal
    for rr=1:size(SCE_traces{cc},2)
        y1{cc}(rr,:) = polyval(p{cc}(rr,:),x1);
    end
end

for cc=1:SCE{ss}.nbSCE  
    %for each SCE transient signal
    for rr=1:size(SCE_traces{cc},2)
        
        %find max time for each ROI in SCE
        [~,idx_max{cc}(rr)] = max(y1{cc}(rr,:)); 
        max_onset{cc}(rr) = x1(idx_max{cc}(rr));
    end
end


%% Plot example with RUN sequence points and order by SCE onset
%start frames of id'd SCE
%sync_idx;
%onsets of each SCE event
%max_onset;
%ROIs identified in SCE
%SCE_ROIs;
%neurons that are part of SCE and involved in RUN sequence
%neurons_participating{6};
%number of the SCE
%nice preplay example: 169
%nice replay example: 262

figure;
for cc=4:79
sce_nb = cc;

%temporary generate times and speed for select input inverval
time = session_vars{ss}.Imaging.time_restricted;
speed = session_vars{ss}.Behavior.speed;

%# of frames before and after onset of SCE
plot_range = [500, 1500];

%absolute time frame of SCE start (all trials) 
start_SCE_frame = SCE{ss}.sync_idx(SCE{ss}.sync_range(sce_nb,1));

%start and end points of sorted (10-15 frame range) - 500 ms = 15
st_evt_sort = start_SCE_frame;
% st_evt_sort = 11511;
end_evt_sort = st_evt_sort+15;


%start idx (absolute)
st_idx =st_evt_sort-plot_range(1);
%end idx (absolute)
end_idx = st_evt_sort+plot_range(2);

%sort all SCE involved neurons
[~,I_sce] = sort(max_onset{sce_nb},'ascend');


imagesc(traces(st_idx:end_idx,SCE{ss}.SCE_unique_ROIs{sce_nb}(I_sce))')
hold on;
title('All SCE ROIS traces sorted by onset time')
colormap('jet')
pause
clf;
end


%get overlap between run sequence neurons and SCE neurons
%[run_sce_neurons, run_sce_idx,~] = intersect(SCE_ROIs{sce_nb},neurons_participating{sce_nb},'stable');

%only ROIs in SCE that are part of run sequence
%run_SCE_ROIs = SCE_ROIs{sce_nb}(run_sce_idx);

%only SCE onsets of ROIS that are part of run sequence
%run_SCE_onsets = max_onset{sce_nb}(run_sce_idx);

%sort only run sequence involved neurons
%[~,I_sce_run] = sort(max_onset{sce_nb}(run_sce_idx),'ascend');

%sort SCE inputs by onset time
%[~,I_sce] = sort(max_onset{sce_nb},'ascend');




%time = session_vars{1, 1}.Behavior.resampled.time;

%[~,select_speed_idx,~] = intersect(time,time_choice,'stable');

%overwrite
%speed = speed(select_speed_idx);

%%
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

%% Correlate SCE activation onset time with median (try mean)

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

end

