function [Events] = detect_events(options, Imaging,updated_dff)

%% Import parameters and C_df

%V2 - adds negative event detection
%V3 - integrates with OCGOL place cell analysis code

%whether to import complete lap restricted or all lap calcium signal and frame time  

% complete laps
if options.restrict==true && options.update_events == 0
    dff = Imaging.trace_restricted;
    dff_time = Imaging.time_restricted;
            
elseif options.restrict==true && options.update_events == 1
    %choose which dff to work with
    if strcmpi(options.dff_type, 'Jia/Danielson')
        dff =  updated_dff.F_df_exp; %Jia
    elseif strcmpi(options.dff_type, 'Rolling median')
        dff =  updated_dff.F_dff_exp; %Rolling median dFF
    end
    
    dff_time = Imaging.time_restricted;
    
elseif options.restrict==false && options.update_events == 1
    %choose which dff to work with
    if strcmpi(options.dff_type, 'Jia/Danielson')
        dff =  updated_dff.F_df_exp; %Jia
        disp('UP_J')
    elseif strcmpi(options.dff_type, 'Rolling median')
        dff =  updated_dff.F_dff_exp; %Rolling median dFF
        disp('UP_R')
    end
    
    dff_time = Imaging.time;
    
    % all laps
elseif options.restrict==false && options.update_events == 0
    dff = Imaging.trace;
    dff_time = Imaging.time;
    disp('START');
end

%minimum duration of the event to be considered significant (s)
mindurevent = options.mindurevent;

%imaging rate - used to calculate the event durations
imaging_rate = options.imaging_rate;

%which ROI to plot
ROI = options.ROI;

%for each dFF type here
%dff = dFF_types{1}; 
%transpose so that columns are ROIs
%dff = dff'; 
%second round
%dff = F_dff_exp_baseline_masked;


%if updated dff struct is not included in the call, make empty
if nargin < 3
    updated_dff = []; 
end

%standard deviation (default) - add variable selection in the future
std_on = 2;
std_off = 0.5;

%determine standard deviation and mean depending on whether it's an initialization
%call or update call

if options.update_events == 1
    if strcmpi(options.dff_type, 'Jia/Danielson')
        trace_std = updated_dff.std_df;
        trace_mean = updated_dff.mean_df;
        
    elseif strcmpi(options.dff_type, 'Rolling median')
        trace_std = updated_dff.std_dff;
        trace_mean = updated_dff.mean_dff;
    end
else
    trace_std = std(dff);
    trace_mean = mean(dff);
end

%expand the std onset threshold to matrix along time
stdThresOn= repmat((trace_std*std_on +trace_mean),size(dff,1),1);
stdThresOnNeg = repmat((trace_mean - trace_std*std_on),size(dff,1),1);

%find all points greater than the std_on
%[ROI,time_idx_lc] = find(dff >= stdThresOn);

%find time indices where the dF/F exceeds the threshold
%for each ROI
for rr=1:size(dff,2)
   %postive events
   onset{rr} = find(dff(:,rr) >= stdThresOn(:,rr));
   %negative events
   onset_neg{rr} =find(dff(:,rr) <= stdThresOnNeg(:,rr));
end

%plot all initially detected peaks
figure;
plot(dff(:,ROI))
hold on;
refline(0,stdThresOn(1,ROI))
stem(onset{ROI}, 2*ones(size(onset{ROI},1),1),'r');

%trace mean
hMean = refline(0,trace_mean(ROI));
hMean.Color = 'c';

%negative std deflection
hnegLine = refline(0,stdThresOnNeg(1,ROI));
hnegLine.Color = 'y';

%negative onsets
stem(onset_neg{ROI}, 2*ones(size(onset_neg{ROI},1),1),'b');

%expand the std offset threshold to matrix along time (postive offset)
stdThresOff = repmat((trace_std*std_off +trace_mean),size(dff,1),1);

% expand the std offset threshold to matrix along time (negative offset)
stdThresOffNeg = repmat((trace_mean - trace_std*std_off),size(dff,1),1);


%make empty storage cells
%find time indices where the dF/F exceeds the threshold
%for each ROI
offset = cell(1,size(dff,2));
offset_abs = cell(1,size(dff,2));
offset_abs_unique = offset_abs;

offset_neg = cell(1,size(dff,2));
offset_abs_neg = cell(1,size(dff,2));
offset_abs_unique_neg = offset_abs_neg;

%fill with empty vectors that correspond to length of each onset vector

for rr =1:size(dff,2)
    %positive events
    offset{rr} = zeros(1,size(onset{rr},1));
    %negative events
    offset_neg{rr} = zeros(1,size(onset_neg{rr},1));
end

%most time consuming - try to optimize - parfor seemed to speed it up 30+
%fold
tic;
%for each ROI
parfor rr=1:size(dff,2)
    disp(rr)
    %for each detected onset positive event
    for ee= 1:size(onset{rr},1)
        offset_temp = find(dff(onset{rr}(ee):end,rr) <= stdThresOff(1,rr),1);
        if isempty(offset_temp)
            offset{rr}(ee) = NaN;
        else
            offset{rr}(ee) = offset_temp;
        end
    end
    
    for ee= 1:size(onset_neg{rr},1)
        offset_temp = find(dff(onset_neg{rr}(ee):end,rr) >= stdThresOffNeg(1,rr),1);
        if isempty(offset_temp)
            offset_neg{rr}(ee) = NaN;
        else
            offset_neg{rr}(ee) = offset_temp;
        end
    end
    
    %add idx to account for find positive offset
    offset_abs{rr} = offset{rr}' + (onset{rr} - 1);
    
    %add idx to account for find positive offset
    offset_abs_neg{rr} = offset_neg{rr}' + (onset_neg{rr} - 1);
    
end

toc; 

%for each id'd peak above the threshold, look for the first value that
%falls below offset STD
% for ii=1:size(lc,2)
%    off_idx(ii) = find(dff(lc(ii):end) < (std_off*trace_std  + trace_mean),1);
%    
%    %correct to absolute frame index
%    %to each idx add the previous onset index
%    off_idx_abs(ii) = off_idx(ii) + lc(ii)-1;
%    
% end

%plot onset and offset
% figure;
% plot(dff(:,ROI))
% hold on;
% hOn = refline(0,stdThresOn(1,ROI));
% hOn.Color = 'r';
% 
% hOff = refline(0,stdThresOff(1,ROI));
% hOff.Color = 'k';
% 
% %onset
% stem(onset{ROI}, 2*ones(size(onset{ROI},1),1),'r');
% %offset
% stem(offset_abs{ROI}, 2*ones(size(offset_abs{ROI},1),1),'k');

%remove duplicate offsets (to get only one onset association with one
%offset

for rr=1:size(dff,2)
    %positive events
    offset_abs_unique{rr} = unique(offset_abs{rr});
    %negative events
    offset_abs_unique_neg{rr} = unique(offset_abs_neg{rr});
end

%remove onsets within actual events 
%take the first onset after the last offset
%every offset should have one corresponding onset for calcium event

%for all offsets
%for each ROI
for rr=1:size(dff,2)
    
    %if not events, set to empty
    if isempty(offset_abs_unique{rr})  || (size(offset_abs_unique{rr},1)-1) == 0
        onset_unique{rr} = NaN;
    else
    %for each unique positive offset
    for ee=1:(size(offset_abs_unique{rr},1)-1)
        %first onset idx after offset
        on_temp = find(onset{rr} > offset_abs_unique{rr}(ee),1);
        if isempty(on_temp)
            onset_unique{rr}(ee+1) = NaN;
        else
            onset_unique{rr}(ee+1) = onset{rr}(on_temp);
        end
    end
    end
    
    %if not events, set to empty
    if isempty(offset_abs_unique_neg{rr})  || (size(offset_abs_unique_neg{rr},1)-1) == 0
        %disp(rr)
        onset_unique_neg{rr} = NaN;
    else
        %for each unique positive offset
        for ee=1:(size(offset_abs_unique_neg{rr},1)-1)
            %first onset idx after offset
            on_temp = find(onset_neg{rr} > offset_abs_unique_neg{rr}(ee),1);
            %disp(rr);
            if isempty(on_temp)
                onset_unique_neg{rr}(ee+1) = NaN;
            else
                onset_unique_neg{rr}(ee+1) = onset_neg{rr}(on_temp);
            end
        end
    end
end

%add the first onset
for rr=1:size(dff,2)
    %positive events
    if ~(numel(onset_unique{rr}) == 1 && isnan(onset_unique{rr}))
        onset_unique{rr}(1)=onset{rr}(1);
    end
    
%     %negative events
    if ~(numel(onset_unique_neg{rr}) == 1 && isnan(onset_unique_neg{rr}))
        onset_unique_neg{rr}(1)=onset_neg{rr}(1);
    end
end

%add the first onset
% onset_idx(1)=lc(1);

%onset offset matrix
for rr=1:size(dff,2)
    %positive events
    onset_offset{rr} = [onset_unique{rr}',offset_abs_unique{rr}];
    %negative events
    onset_offset_neg{rr} =  [onset_unique_neg{rr}',offset_abs_unique_neg{rr}];
end

%remove NaN events from onset_offset_f
for rr=1:size(dff,2)
    %positive events
     nan_temp = logical(sum(isnan(onset_offset{rr}),2));
     onset_offset{rr}(nan_temp,:) = [];
     
     %negative events
      nan_temp_neg = logical(sum(isnan(onset_offset_neg{rr}),2));
     onset_offset_neg{rr}(nan_temp_neg,:) = [];
end

%onset offset matrix
%onset_offset_fr = [onset_idx',off_idx_unique'];

%plot onset and offset
figure;
plot(dff(:,ROI))
hold on;
hOn = refline(0,stdThresOn(1,ROI));
hOn.Color = 'r';

hOff = refline(0,stdThresOff(1,ROI));
hOff.Color = 'k';

if ~isempty(onset_offset{ROI})
    
    %positive onset
    stem(onset_offset{ROI}(:,1), 2*ones(size(onset_offset{ROI}(:,1),1),1),'r');
    %positive offset
    stem(onset_offset{ROI}(:,2), 2*ones(size(onset_offset{ROI}(:,2),1),1),'k');
end

if ~isempty(onset_offset_neg{ROI})
    %negative onset
    stem(onset_offset_neg{ROI}(:,1), 2*ones(size(onset_offset_neg{ROI}(:,1),1),1),'b');
    %negative offset
    stem(onset_offset_neg{ROI}(:,2), 2*ones(size(onset_offset_neg{ROI}(:,2),1),1),'k');
end

%% Calculate the amplitude (in sigma) and duration of events

%preallocate
duration_pos = cell(1,size(dff,2));
duration_neg = cell(1,size(dff,2));
amplitude_pos_dFF = cell(1,size(dff,2));
amplitude_pos_sigma = cell(1,size(dff,2));
amplitude_neg_dFF = cell(1,size(dff,2));
amplitude_neg_sigma = cell(1,size(dff,2));
 
%duration of each event (s)
for rr=1:size(dff,2)
    if ~isempty(onset_offset{rr})
        duration_pos{rr} = (onset_offset{rr}(:,2) - onset_offset{rr}(:,1))*(1/imaging_rate);
    end
    
    if ~isempty(onset_offset_neg{rr})
        duration_neg{rr} = (onset_offset_neg{rr}(:,2) - onset_offset_neg{rr}(:,1))*(1/imaging_rate);
    end
end

%amplitude of each event 
%take the max value of each event (from 2sigma to 0.5sigma) and scale to sigma
for rr=1:size(dff,2)
    if ~isempty(onset_offset{rr})
        %for each event
        for ee=1:size(onset_offset{rr})
            amplitude_pos_dFF{rr}(ee,1) = max(dff(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr));
            amplitude_pos_sigma{rr}(ee,1) = amplitude_pos_dFF{rr}(ee,1)./((trace_std(rr)));
        end
    end
    
    if ~isempty(onset_offset_neg{rr})
        %for each event
        for ee=1:size(onset_offset_neg{rr})
            amplitude_neg_dFF{rr}(ee,1) = min(dff(onset_offset_neg{rr}(ee,1):onset_offset_neg{rr}(ee,2),rr));
            amplitude_neg_sigma{rr}(ee,1) = amplitude_neg_dFF{rr}(ee,1)./((trace_std(rr)));
        end
    end
    
end

%% Bin positive and negative events into 250 ms bins

%reshape and convert to matrices ampltudes and durations for positive and
%negative events
duration_pos_reshaped = cell2mat(reshape(duration_pos,[],1));
duration_neg_reshaped = cell2mat(reshape(duration_neg,[],1));

amplitude_pos_reshaped = cell2mat(reshape(amplitude_pos_sigma,[],1));
amplitude_neg_reshaped = cell2mat(reshape(amplitude_neg_sigma,[],1));

%find max amplitude and duration
maxDuration = max([duration_pos_reshaped; duration_neg_reshaped]);
maxAmplitude = max([amplitude_pos_reshaped; amplitude_neg_reshaped]);

%0 to 5s every 250 ms
edgesTime = 0:0.25:ceil(maxDuration);

%bin amplitude of event into 0.5 sigma bins
edgesSigma = 0:0.5:ceil(maxAmplitude);

centers = 0.25:0.25:5;
%points for smooth exponential
centers_exp = 0.25:0.05:5;

%for fixed sd below
sigmaThres = [2 3 4];

%time edges for normal/fixed exp fitting
edges = 0:0.25:5;

%bin data in both duration and amplitude
%positive events
[N2d_pos,~,~,bin2dTime_pos,bin2dSigma_pos] = histcounts2(duration_pos_reshaped,amplitude_pos_reshaped,edgesTime,edgesSigma);

%negative events
[N2d_neg,~,~,bin2dTime_neg,bin2dSigma_neg] = histcounts2(duration_neg_reshaped,-amplitude_neg_reshaped,edgesTime,edgesSigma);

%negative to positive event ratio
event_ratio_2d = N2d_neg./N2d_pos;

%find bins with FPR < 0.05
event_ratio_2d_select = zeros(size(event_ratio_2d,1),size(event_ratio_2d,2));
event_ratio_2d_select(find(event_ratio_2d < 0.05)) = 1;

%score the events (binary mask to include vs exclude)
%event_ratio - rows - time bins; columns - ampltiude bins
%for check each event for significance
timeSigmaBin = [bin2dTime_pos, bin2dSigma_pos];

%positive events only - should match index of onset_offset events in cells
include_events = logical(zeros(size(bin2dTime_pos,1),1));

for ii = 1:size(bin2dTime_pos,1)
    if(event_ratio_2d_select(timeSigmaBin(ii,1),timeSigmaBin(ii,2)) == 1)
        include_events(ii) = 1;
    end
end


%% Filter positive events based on time

%if on last iteration of the algorithm, filter the events
if options.iterationNb == options.it
    
    disp('Excluding events based on duration threshold');    
    %find indices of events below the duration threshold
    eventsBelowTimeThresIdx = find(duration_pos_reshaped < mindurevent);
    
    %exclude them in the logical filter
    include_events(eventsBelowTimeThresIdx) = 0;
end

%% Display the fraction of events included
eventsFormat = 'Events included: %d / %d';
eventsOutput = sprintf(eventsFormat, length(find(include_events == 1)), size(bin2dTime_pos,1));
disp(eventsOutput);

%% Remove the discarded events from the original events

%set to nan the onset/offset frames in the original cell 
%preallocate
onset_offset_filtered = onset_offset; 
onset_offset_cleared = onset_offset;

amplitude_pos_dFF_filtered = amplitude_pos_dFF;
amplitude_pos_dFF_cleared = amplitude_pos_dFF;

duration_pos_filtered = duration_pos;
duration_pos_cleared = duration_pos;

amplitude_pos_sigma_filtered = amplitude_pos_sigma;
amplitude_pos_sigma_cleared = amplitude_pos_sigma;

%cell(size(onset_offset,1),size(onset_offset,2));
%get size of each cell in original onset/offset cell
for rr=1:size(onset_offset,2)
	event_size(rr) = size(onset_offset{rr},1);
end

%should equal to the include_events vector
isequal(sum(event_size), size(include_events,1));

%acumulation of event size for ROI index switch
cum_event_size_end = cumsum(event_size);
cum_event_size_start_inc = cum_event_size_end+1;
cum_event_size_start = [1, cum_event_size_start_inc(1:(end-1))];

%for each ROI
for rr=1:size(onset_offset,2)
       %check that the range is not zero before calculating the events
       if ~isempty(cum_event_size_start(rr):cum_event_size_end(rr))
    
       %preserve the order by inserting a nan in their place
       onset_offset_filtered{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = nan;
       %lose the indices by setting elements to empty
       onset_offset_cleared{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = [];
       
       %filter amplitudes and durations
       
       %amplitude in dFF
       amplitude_pos_dFF_filtered{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = nan;
       amplitude_pos_dFF_cleared{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = [];
       
       %duration in units sigma
       amplitude_pos_sigma_filtered{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = nan;
       amplitude_pos_sigma_cleared{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = [];
       
       %duration (s)
       duration_pos_filtered{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = nan;
       duration_pos_cleared{rr}(~include_events(cum_event_size_start(rr):cum_event_size_end(rr)),:) = [];
       
       end

end

%% Create a cell with corresponding imaging times

onset_offset_time = onset_offset_cleared;

%for each ROI
for rr=1:size(dff,2)
    %for each event
    for ee= 1:size(onset_offset_cleared{rr},1)
        onset_offset_time{rr}(ee,:) = dff_time(onset_offset_cleared{rr}(ee,:));
    end
end

%% for each sigma threshold, select events (old)
for ii=1:size(sigmaThres,2)
    
    %for positive events
    amp_pos_events_idx{ii} = find(amplitude_pos_reshaped >= sigmaThres(ii));
    amp_pos_events{ii} = amplitude_pos_reshaped(amp_pos_events_idx{ii});
    dur_pos_event{ii} = duration_pos_reshaped(amp_pos_events_idx{ii});
    
    %for negative events
    amp_neg_events_idx{ii} = find(amplitude_neg_reshaped <= -sigmaThres(ii));
    amp_neg_events{ii} = amplitude_neg_reshaped(amp_neg_events_idx{ii});
    dur_neg_event{ii} = duration_neg_reshaped(amp_neg_events_idx{ii});
    
    %counts in each bin for positive events
    [N_dur_pos{ii},edges,bin] = histcounts(duration_pos_reshaped(amp_pos_events_idx{ii}),edges);
    
    %counts in each bin for negative events
    [N_dur_neg{ii},edges,bin] = histcounts(duration_neg_reshaped(amp_neg_events_idx{ii}),edges);

end

%% Fit decaying exponentials to the FPR values

%generates errors

%for each sigma level
% for ii =1:3
%     %get FPR ratio at each time bin
%     fpr_ratio{ii} = N_dur_neg{ii}./N_dur_pos{ii};
%     %only include no nan values in fit
%     nonan_logical{ii} = ~isnan(fpr_ratio{ii});
%     %fit decaying exponential
%     fe= fit(centers(nonan_logical{ii})',fpr_ratio{ii}(nonan_logical{ii})','exp1');
%     %save the model parameters at each sigma
%     exp_model_param{ii} = [fe.a, fe.b];
%     %figure;
%     %plot the exponential fit at each sigma level
%     %plot(fe,centers(nonan_logical{ii})',fpr_ratio{ii}(nonan_logical{ii})');
% end

%figure;
%plot(centers_exp,exp_model_param{1}(1)*exp(exp_model_param{1}(2)*centers_exp))

%% Plot FPR and exponential fit at each sigma amplitude
%add histogram distributions as well
% 
% figure;
% colorAssn = {'r','g','b'};
% 
% hold on
% xlim([0 4]);
% ylim([0 0.25]);
% 
% for ii=1:size(sigmaThres,2)
%     p{ii} = scatter(centers,N_dur_neg{ii}./N_dur_pos{ii},colorAssn{ii});
%     e{ii} = plot(centers_exp,exp_model_param{ii}(1)*exp(exp_model_param{ii}(2)*centers_exp),colorAssn{ii});
% end

%plot reference line at p = 0.05
% refline(0,0.05)
% 
% hold off
% legend([p{1},p{2},p{3}],'2 \sigma','3 \sigma','4 \sigma')

%% Plot histograms - add this later


%% Create onset_binary and onset_ones for compatibility from filtered events (excluded)


onset_binary = zeros(size(dff,1),size(dff,2));
onset_ones = zeros(size(dff,1),size(dff,2));

%onset ones - 1's out all the frames during which the event is ongoing
%for each ROI
%onset binary - 1's onset at the onset of the event
for rr=1:size(dff,2)
    %for each event
    for ee= 1:size(onset_offset_cleared{rr},1)
        onset_ones(onset_offset_cleared{rr}(ee,1):onset_offset_cleared{rr}(ee,2),rr) = 1;
        onset_binary(onset_offset_cleared{rr}(ee,1),rr) = 1;
    end
end


%% Display only filtered events

%cleared events plotted and negative events

figure;
plot(dff(:,ROI))
hold on;
hOn = refline(0,stdThresOn(1,ROI));
hOn.Color = 'r';

hOff = refline(0,stdThresOff(1,ROI));
hOff.Color = 'k';

if ~isempty(onset_offset_cleared{ROI})
    
    %positive onset
    stem(onset_offset_cleared{ROI}(:,1), 2*ones(size(onset_offset_cleared{ROI}(:,1),1),1),'r');
    %positive offset
    stem(onset_offset_cleared{ROI}(:,2), 2*ones(size(onset_offset_cleared{ROI}(:,2),1),1),'k');
end

if ~isempty(onset_offset_neg{ROI})
    %negative onset
    stem(onset_offset_neg{ROI}(:,1), 2*ones(size(onset_offset_neg{ROI}(:,1),1),1),'b');
    %negative offset
    stem(onset_offset_neg{ROI}(:,2), 2*ones(size(onset_offset_neg{ROI}(:,2),1),1),'k');
end


%% Save ampltiudes, durations, and onset and offset frames into struct

%for filtered events
Events.onset_offset = onset_offset_cleared;
Events.onset_offset_time = onset_offset_time;
Events.duration = duration_pos_cleared;
Events.properties.peak = amplitude_pos_dFF_cleared;
Events.properties.peak_sigma = amplitude_pos_sigma_cleared;

%these are already calculated based on cleared onset_offset
Events.onset_binary = onset_binary;
Events.onset_ones = onset_ones;
Events.options = options;

%2 and 0.5 sigma threshold above which event is considered
Events.std_thresholdOn = stdThresOn(1,:);
Events.std_thresholdOff = stdThresOff(1,:);

%standard deviation of each trace
Events.STD_noise = trace_std; 

%calcium trace time
Events.dff_time = dff_time;

%onset_offset_time




end

