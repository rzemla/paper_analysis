%% Cell activation onset for detected RUN sequences (Vilette Neuron 2015)

%indices of recurring neurons
recurring_neuron_idx

onset_input = session_vars{1, 1}.Imaging_split{1, 4}.trace_restricted;

%make raw copy
pca_input_raw = onset_input;

%s=2s filtering gaussician signal

%% Smooth calcium traces with gaussian with sigma = 2s

%5 sigma gaussian kernel
%5 seconds -  150 frames
%2 seconds -  60 frames
options.sigma_filter =  60;
gaussFilter = define_Gaussian_kernel(options);

for rr=1:size(onset_input,2)
    onset_input_sm(:,rr) =conv(onset_input(:,rr),gaussFilter, 'same');
end

%% Select run sequence neurons
%raw
traces = onset_input(:,recurring_neuron_idx);
%smoothed
run_seq_traces = onset_input_sm(:,recurring_neuron_idx);

%% Take first derivative of each trace and find max

deriv_traces = diff(run_seq_traces,1,1);
%zero pad 1 time point a start
deriv_traces = [zeros(1, size(deriv_traces,2));deriv_traces];

%% Discard onsets with maximal derivative less than 5% dF/F s-1 (s)
%may need to translate the derivative to the raw input to evalulate this
%condition


%% Get indices of lap starts based on position diff
diff_pos = diff(position_norm);
%zero pad
diff_pos = [0; diff_pos];
%find transition idx when position goes from 1 --> 0
%these are the lap starts
idx_laps = find(diff_pos <= -0.95);
%frames onset and offset for each lap
lap_on_off = idx_laps;
lap_on_off = [1;lap_on_off];
%define off laps indices
off_laps = idx_laps-1;
off_laps = [off_laps; size(position_norm,1)];
lap_on_off = [lap_on_off,off_laps];


%% Plot norm position and start and end points
figure
hold on
%plot position
plot(position_norm,'k')
%start
stem(lap_on_off(:,1), ones(1,size(lap_on_off,1)),'g')
%end
stem(lap_on_off(:,2), ones(1,size(lap_on_off,1)),'r')

%% Find max of first derivative on each lap and get median (save both)

%find max of derivative for each RUN neuron inside lap
%for each lap
for ll=1:size(lap_on_off,1)
    %return max index within that segment
    [~,max_lap(ll,:)] = max(deriv_traces(lap_on_off(ll,1):lap_on_off(ll,2),:),[],1);
    max_pos{ll} = position_norm(lap_on_off(ll,1):lap_on_off(ll,2));
    %extract the position of the max derivative in each lap
    max_pos_val(ll,:) = max_pos{ll}(max_lap(ll,:));
end
%get median normalized position
median_run_seq_onset = median(max_pos_val,1);



%% Plot original, smoothed, traces and first derivative below
figure
%first 20
for rr=1:20
    
    subplot(2,1,1)
    hold on
    title('Original calcium trace and Gaussian smoothed');
    plot(traces(:,rr),'k');
    plot(run_seq_traces(:,rr),'r');
    hold off
    subplot(2,1,2)
    hold on
    title('first derivative');
    plot(deriv_traces(:,rr),'b');
    
    pause(0.05)
    clf
end


%% Plot according to median activation sequence
[val_temp,idx_med_onset] = sort(median_run_seq_onset,'ascend');

figure;
imagesc(traces(:,idx_med_onset)');
hold on;
colormap('jet')
caxis([0 1.5])


