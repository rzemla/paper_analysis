
%% Input data
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
ROI = 1;
figure;
hold on
plot(input_data(:,ROI), 'k')
plot(filtered_traces(:,ROI)-1,'r');

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
%%
%try different ROIs:
ROI = 4;

%% Plot trace, median and irq range - works
figure
hold on
stepSize = 2;
step = 0;
for ii =1:16
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
plot(norm_position(st_idx:end_idx)-step+2,'r');

%% 

%% Plot traces as line plot (vs imagesc)
figure;
subplot(1,2,1)
hold on;
ylim([-0.5 1])
title('Savitzky-Golay filtered calcium traces')
stepSize = 2;
step = 0;
for ii=1:1%size(filtered_traces,2)
    plot(filtered_traces(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end

%ylabel('Normalized Position');
%plot(norm_position(st_idx:end_idx)-step,'r');

subplot(1,2,2)
hold on;
ylim([-0.5 1])
title('Non-filtered calcium traces')
stepSize = 2;
step = 0;
for ii=1:1%size(filtered_traces,2)
    plot(input_data(:,ii)-step, 'k', 'LineWidth', 1.5)
    step = step - stepSize;
end


