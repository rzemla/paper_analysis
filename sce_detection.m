
%input data (time x ROI)
input_data = traces;

input_data = traces(st_idx:end_idx,input_neuron_idxs_sorted);

%% 3rd order Savitzky-Golay filter frame size 500ms (~15 frames)

%order of filter
sg_order = 3;
%window of filter
sg_window = 15;
%filter along first dimension (across rows)
filtered_traces = sgolayfilt(input_data,sg_order,sg_window,[],1);

%plot side by side (sample ROI
ROI = 3;
figure;
hold on
plot(input_data(:,ROI), 'k')
plot(filtered_traces(:,ROI)-1,'r');


%% Plot traces as line plot (vs imagesc)
figure;
subplot(1,2,1)
hold on;
stepSize = 2;
step = 0;
for ii=1:size(filtered_traces,2)
    plot(filtered_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end

%ylabel('Normalized Position');
plot(norm_position(st_idx:end_idx)-step,'r');

subplot(1,2,2)
hold on;
stepSize = 2;
step = 0;
for ii=1:size(filtered_traces,2)
    plot(input_data(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end


