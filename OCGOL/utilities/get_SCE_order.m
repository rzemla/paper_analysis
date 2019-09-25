function [SCE] = get_SCE_order(traces,SCE,ss)


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

%remove noisy SCE events from  r, SCE_traces, SCE{ss}.SCE_unique_ROIs{cc}
%make copies of original values associated with each ROI
r_orig = r;
SCE_traces_orig = SCE_traces;
SCE_unique_ROIs_noise_rm  = SCE{ss}.SCE_unique_ROIs;

%remove selected ROIs from each variable 
for cc=1:SCE{ss}.nbSCE  
    r{cc}(max_sync{cc},:) = [];
    SCE_unique_ROIs_noise_rm{cc}(max_sync{cc}) = [];
    SCE_traces{cc}(:,max_sync{cc}) = [];
end

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

%% Sort all SCE involved neurons

%for each SCE in the session
for cc=1:SCE{ss}.nbSCE
    [~,sce_sort_order{cc}] = sort(max_onset{cc},'ascend');
    %re-sort the SCE-associated ROIs based on this order
     SCE_unique_ROIs_sorted{cc} = SCE_unique_ROIs_noise_rm{cc}(sce_sort_order{cc});
end

%% Add to SCE struct for particular session

SCE{ss}.SCE_unique_ROIs_sorted = SCE_unique_ROIs_sorted;
SCE{ss}.sce_sort_order = sce_sort_order;


end

