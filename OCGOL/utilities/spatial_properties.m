function [Place_cell]= spatial_properties(Behavior, Events, Imaging, options)

%define the spatial properties and inputs for subsequent spatial analysis


%% Import data

%imaging period
dt = Imaging.dt;

%uses detrended baseline with matlab msbackadj function
if Events.options.restrict==true
    C_df=Imaging.trace_restricted;
    
elseif Events.options.restrict==false
    C_df=Imaging.trace; 
end

%RUN ACTIVITY
%import downsampled binary run activity from run_epoch script
%1' indicate that the animal is running during that time - RZ
run_ones=Behavior.run_ones;

%downsampled time and position taken from run epochs only - RZ
runtime=Behavior.resampled.run_time;

%matches the number of indices of the taken from run epochs only - RZ
run_position=Behavior.resampled.run_position; 

%sampled lap index (lap index during run epochs only)
run_lapNb = Behavior.resampled.run_lapNb;

%EVENT ACTIVITY
%columns of onset and offset indices when animal was running and a significant event
%occured
run_onset_offset_ones=Events.Run.run_onset_ones;

%event onsets during run activity for each ROI
run_onset_binary=Events.Run.run_onset_binary;

%run_binary=Behavior.runbinary;

%% Convolution kernel for Gaussian smoothing

[gaussFilter] = define_Gaussian_kernel(options);

%% Mask C_df for significant events detected - set non significant frames to nan
%turn into logical for masking non run, not sig dF/F traces for Barthos PF
%detection code
mask_nonsig_dFF = ~logical(run_onset_offset_ones);

%take lap restricted C_df and set nonsig no run dF/F parts of trace to nan
%make copy
C_df_masked = C_df; 
%set nan (every timepoint that is not part of a significant run event
C_df_masked(mask_nonsig_dFF) = nan;

%% Bin data spatially - binned by run only

%Time and position for RUNNING epochs ONLY:
%Bin normalized running position
run_position_norm = Behavior.resampled.run_position_norm;

%Bin running position:
%for each number of bins, bin the normalized position during run epochs
for ii=1:length(options.Nbin)
    %histcounts(input vector, # of bins) 
    %count_bin - # of values in each bin 
    %edges - vector of the edges for each bin = #of bins + 1 
    %bin - which bin each normalized running position belongs to 
    [count_bin{ii},edges{ii},bin{ii}] = histcounts(run_position_norm, options.Nbin(ii));
end

%% Mask with only run epoch bin assignments

%export bin mask that corresponds to index of resticted frames
binFrameMask = zeros(1,size(C_df,1));
%all run indices 
runIdx = find(run_ones ==1);
%assign the run time points their respecitive spatial bins (100 bins = 8)
binFrameMask(runIdx) = bin{8};
%set all else to nan
binFrameMask(find(binFrameMask==0)) = nan;


% # of event onsets for each bin
% dF/F for each bins
%% 1)only run epochs

%run_onset_binary - event onsets during run activity for each ROI
%frames x ROIs
%for all ROIs, take frames that match when animal is running
%run_ones==1 --> turns double vector to logical and can use to select
%indices in matrices
%vector length should match match runtime and run_position
%not binned at this time, but used in binning analysis
run_onset_bin = run_onset_binary(run_ones==1,:);

%significant events (contiguous) during run epochs only
run_onset_offset_runEpochs = run_onset_offset_ones(run_ones==1,:);

%do the same for the C_df trace (including non-sig event dF/F)
run_C_df = C_df(run_ones==1,:);

%take out only run epochs from sig event run filtered dF/F traces
run_sig_C_df = C_df_masked(run_ones==1,:);


%plot and compare both C_df matrices
figure
subplot(1,2,1) 
imagesc(run_C_df');
title('Without nonsig events removed')
subplot(1,2,2) 
imagesc(run_sig_C_df');
title('With nonsig events removed')
hold off

%% 2)value for each bin onset
%for each ROI, find in which bin did each event onset occur when the animal
%was running
%bin - which bin each position/time bin the animal belongs to
%bin - is a cell of vectors which contain the # of the bin the frame
%belongs to

%for each binning size
%bb - each binning
%rr - each ROI
for bb=1:length(options.Nbin)
    %for each ROI
    for rr=1:size(run_onset_bin,2)
        %find the corresponding bin #s of each event that occurred the run
        %epochs
        %run_onset_bin(:,n)==1 - turns to logical run frames where an event
        %onset occurred
        %returns the bin values of where the running onset occured for each
        %ROI
        %RUNNING RELATED ONSET TRANSIENTS FOR EACH ROI BY BIN
        onset_bin{bb}{rr}=bin{bb}(run_onset_bin(:,rr)==1);
        %for contigous events (frames during which a sig event occurred)
        onset_bin_conti{bb}{rr}=bin{bb}(run_onset_offset_runEpochs(:,rr)==1);
    end
end

%return the corresponding lap of the bin on which the event occurred
for bb=1:length(options.Nbin)
    %for each ROI
    for rr=1:size(run_onset_bin,2)
        %find associated absolute lap index for each run event 
        onset_lap_bin{bb}{rr} = run_lapNb(run_onset_bin(:,rr)==1);
    end
end

%2)value for each bin dF/F

%for each bin range
for bb=1:length(options.Nbin)
    %for each bin in bin range
    for binN=1:options.Nbin(bb)
        %for each ROI
        for rr=1:size(run_onset_bin,2)
            %for calcium activity for each ROI during running epoch, find
            %indices of each bin and take dF/F activity in that bin
            %for each binning range, for each bin in that range, return
            %dF/F values in chronological order
            %run_C_df --> frames( when animal is running) x ROI # 
            %find(bin{i}==binN - find all frames associated with specific
            %bin --> move C_df values in those running bins for each ROI
            %into separate cells
            %takes dF/F values regardless of significance in this case
            %one mean dF/F per bin range, per bin, per ROI
            dF_bin{bb}{binN}(:,rr)=run_C_df(find(bin{bb}==binN),rr);
            %with non_sig dF/F removed
            dF_sig_bin{bb}{binN}(:,rr)=run_sig_C_df(find(bin{bb}==binN),rr);
        end
        %take mean of dF/F values for each bin range, in each bin for each ROI
        %takes mean over each column which corresponds to each ROI
        %over each interation of bin range and bin
        dF_map{bb}(binN,:)=nanmean(dF_bin{bb}{binN});
        %with non_sig dF/F removed
        dF_sig_map{bb}(binN,:)=nanmean(dF_sig_bin{bb}{binN});
    end
end

%smooth the significant running dF/F overneighboring bins
%smooth each map with a mean neighboring point moving window average
%size of window for input into PF finder
smoothWin_PF = 3;

dF_sig_map_100_smooth = [];
%for each ROI - use 100 bin map
for rr = 1:size(dF_sig_map{8},2)
    %smooth each ROI trace with moving average filter
    dF_sig_map_100_smooth(:,rr) = smooth(dF_sig_map{8}(:,rr),smoothWin_PF,'moving');
end


%% OCCUPANCY CALCULATION
%imaging period imported
avg_fr = dt;

%for each bin range
for bb=1:length(options.Nbin)
    %how many frames in each bin when animal is running for each bin range
    %multiple each frame count in each in by dt of each frame
    occupancy_time_s{bb}=count_bin{bb}*avg_fr; % time spend (in sec)
    %probability that the animal is the given bin during running
    %frames in each bin/total # of frames
    proba_bin{bb}=count_bin{bb}/sum(count_bin{bb}); % time spend / total time
end

%same regardless of how you bin the data - using first bin here
%should equal total # of run frames * avg_fr
%length(runtime)*avg_fr
%NOT USED IN CODE BELOW
total_occupancy=sum(count_bin{1})*avg_fr; %total time
    

% Position for each bin dF/F
%Find center position of each bin in cm (1 value per bin) based on running
%epochs
%for each bin range
for i=1:length(options.Nbin)
    %for each bin
    for binN=1:options.Nbin(i)
        %from position of animal during run epochs, take positions for
        %identified bin locations
        %bin - which bin each frame belongs to
        %find those that correspond to current bin (binN) and move to cell
        pos_bin{i}{binN}=run_position(find(bin{i}==binN),1);
        %take average of all the positions in each bin and assign to each
        %bin in vector location for bin range
        mean_pos_bin{i}(binN)=mean(pos_bin{i}{binN});
    end
end

%% Extract mean dF/F in each bin by lap for given ROI - copy of the above code 

%lap indices present in the analysis
lapIdx =unique(run_lapNb);

for i=1:length(options.Nbin)
    %for each bin
    for binN=1:options.Nbin(i)
        %for each lap present
        for ll=1:size(lapIdx,1)
            %for each ROI
            for n=1:size(run_onset_bin,2)

                %for calcium activity for each ROI during running epoch, find
                %indices of each bin and take dF/F activity in that bin
                %(frames)
                
                %for each binning range, for each bin in that range, return
                %dF/F values in chronological order
                
                %run_C_df --> frames( when animal is running) x ROI #
                
                %find(bin{i}==binN - find all frames associated with specific
                %bin --> move C_df values in those running bins for each ROI
                %into separate cells
                
                %takes dF/F values regardless of significance in this case
                %one mean dF/F per bin range, per bin, per ROI
                %added by lap - SLOW - try preallocating array
                dF_bin_lap{i}{ll}{binN}(:,n)=run_C_df(find(bin{i}==binN & run_lapNb==lapIdx(ll)),n);
                
            end
            
        %take mean of dF/F values for each bin range, in each bin for each ROI
        %takes mean over each column which corresponds to each ROI
        %over each interation of bin range and bin
        dF_lap_map{i}{ll}(binN,:)= mean(dF_bin_lap{i}{ll}{binN},1); 
        
        %dF_map{i}{ll}(binN,:)= mean(dF_bin{i}{ll}{binN},1);  
        end

    end
end

%% Calculate the occupancy each lap for each bins (100 bins)

%for each lap
for ll=1:size(dF_bin_lap{8},2)
    %for each bin
    for bb=1:size(dF_bin_lap{8}{ll},2)
        %in # of frames
        lap_occupancy{ll}(bb,1) = size(dF_bin_lap{8}{ll}{bb},1);
        %in seconds
        lap_occupancy_s{ll}(bb,1) = lap_occupancy{ll}(bb,1).*avg_fr;
    end
    
    lap_occupancy_s_mat(ll,:) = lap_occupancy_s{ll}';
end

%% Construct raster of the mean dF/F across 100 bins for each lap for each ROI
%and store in a cell

%for each ROI
for rr=1:size(C_df,2)

    %for each lap
    for ll=1:size(dF_lap_map{8},2)
        
        %only mean dF/F
        dF_lap_map_ROI{rr}(ll,:) = dF_lap_map{8}{ll}(:,rr);

    end
    
    %oocupancy normalized mean dF/F
    dF_lap_map_ROI_norm{rr} = dF_lap_map_ROI{rr}./lap_occupancy_s_mat;
end

%smooth each map with a mean neighboring point moving window average
%size of window
smoothWin = options.smooth_span;

%for each ROI
for rr=1:size(C_df,2)
    %for each lap
    for ll=1:size(dF_lap_map{8},2)
        dF_lap_map_ROI_smooth{rr}(ll,:) = smooth(dF_lap_map_ROI{rr}(ll,:)',smoothWin,'moving');
    end
end


figure
imagesc(dF_lap_map_ROI_smooth{options.c2plot})
hold on;
caxis([-0.2 2.5])
colormap('jet')

%Save to struct

%non-occupancy normalized (100 bins)
Place_cell.dF_lap_map_ROI = dF_lap_map_ROI;

%non-occupancy normalized (100 bins) smoothed with average moving filter
%(Sheffield 2015)
Place_cell.dF_lap_map_ROI_smooth =  dF_lap_map_ROI_smooth;

%occupancy normalized (mean dF/F/s)
Place_cell.dF_lap_map_ROI_norm = dF_lap_map_ROI_norm;

%% Rate maps: 
%total number of onsets that occurred in a location across all laps divided 
%by the amount of time the animal spent there during run epochs - RZ

%bin divided by the time the animal spent in that bin (during running
%epochs)

% onset map / occupancy
%Smooth rate map:
%Onset counts (onset map) and occupancy times in each bin were are independently smoothed 
%by convolving with a Gaussian smoothing kernel
% kernel = 3 bins Danielson et al.
% try other filter...

%Smooth dF/F (Dombeck et al. 2010 used moving average span=3bins)
%onset bin - for each bin range, for each ROIs, bins in which it fired

%for each bin range
for i=1:length(options.Nbin)
    %for each ROI
    for n=1:size(run_onset_bin,2)
        %for each bin position
        for binN=1:options.Nbin(i)
            %for each ROI in each bin (and bin range), find # of onsets in
            %each bin - RZ
            %find() won't return NaNs because value always being comapred
            %to integer bin #; numel will not generate NaN counts
            %fill in each element of 2D matrix as bin # x ROI
            onset_map{i}(binN,n)=numel(find(onset_bin{i}{n}==binN));
            
            %number of frames during which significant events were occured
            %in each bin
            onset_map_conti{i}(binN,n)=numel(find(onset_bin_conti{i}{n}==binN));
        end
        %filter
        %using 2D Gaussian filter function - see if result the same if
        %convolving using a 1D Gaussian kernel - RZ
        %aligns quite well with event count using imgaussfilt - RZ
        %smooth each ROI over all bins with given bin range - Gaussian
        %filter
        %onset_map_sm{i}(:,n)=imgaussfilt(onset_map{i}(:,n),options.sigma_filter);
        %replace with convoltuon filter here
        onset_map_sm{i}(:,n) = conv(onset_map{i}(:,n),gaussFilter, 'same');
        
        %moving average filter - take symmetrical neighboring values (bins) around
        %values of interest and divide by span (moving window)
        dF_map_sm{i}(:,n)=smooth((dF_map{i}(:,n)),options.smooth_span);
    end
    %occupancy time is also smoothed with a Gaussian filter - WHY?
    %see if you get the same result with using conv against with Gaussian
    %occupancy_time_s_sm{i} = imgaussfilt(occupancy_time_s{i},options.sigma_filter);
    %replace with convolution filter here
    occupancy_time_s_sm{i} = conv(occupancy_time_s{i},gaussFilter, 'same');
    
    %gaussian unsmoothed rate map divide by unsmoothed occupancy - RZ
    %occupancy diminishes the rate the location of the animals b/c it spends
    %most of the time running above the speed threshold there
    rate_map{i}=onset_map{i}./occupancy_time_s{i}';
    %both smoothed with imgaussfilt fxn - RZ
    %rate_map_sm{i}=onset_map_sm{i}./occupancy_time_s_sm{i}';
    rate_map_sm{i}=onset_map_sm{i}./occupancy_time_s{i}';
    %rate map for sorting and display purposes
    rate_map_sm_plot{i} = onset_map_sm{i}./occupancy_time_s_sm{i}';
    
end

%Gaussian smoothing of rate map (post loop) - RZ (not sure why MD smoothed transient
%count and occupancy time separately
for i=1:length(options.Nbin)
    %for each ROI
    for n=1:size(run_onset_bin,2)
        rate_map_sm_whole{i}(:,n) = conv(rate_map{i}(:,n),gaussFilter, 'same');
    end
end

%% Overall firing rate - overall rate and mean rate
%Mean firing rate - this is based on the smoothed rate (not actual)! - RZ
%rate map not normalized yet

%for each bin range
for i=1:length(options.Nbin)
    %sum smoothed rate across all bin and divide by total time of run
    %smoothed occupancy varies per bin range; smoothed rate map spills to
    %neighboring bins
    overall_rate_sm{i}=sum(onset_map_sm{i})/sum(occupancy_time_s_sm{i});
    %RZ non smoothed rate
    overall_rate{i} = sum(onset_map{i})/sum(occupancy_time_s{i});
    
    mean_rate_sm{i} = mean(rate_map_sm{i});
    mean_rate{i} = mean(rate_map{i});
    
    %Skaggs definition - smaller but same scale as mean_rate
    mean_rate_Sk{i} = sum(rate_map{i}.*proba_bin{i}');
end

%% Normalize rate maps (smoothed)
%Normalize rate map and dF map - RZ - min/max range
%for each bin range
for i=1:length(options.Nbin)
    %for each ROI
    for n=1:size(rate_map_sm{i},2)
        %normalize each ROI against min and max value in each bin range
        %using previous smoothed rate map
        
        norm_rate_map_sm{i}(:,n)=( rate_map_sm{i}(:,n) - min(rate_map_sm{i}(:,n)) )/( max(rate_map_sm{i}(:,n)) - min(rate_map_sm{i}(:,n)) );
        norm_dF_map_sm{i}(:,n)=( dF_map_sm{i}(:,n) - min(dF_map_sm{i}(:,n)) )/( max(dF_map_sm{i}(:,n))-min(dF_map_sm{i}(:,n)) );
    end
end

%% Spatial tuning curve (=normalized rate map) = STC 
%store # of spatial bins used to construct tuning curve
bin_STC=find(options.bin_spatial_tuning==options.Nbin);

%check that the value is correct
if isempty(bin_STC)
    disp('change value of "options.bin_spatial_tuning" to a value present in "options.Nbin"')
    return
end

%select the curve from the normalized and smoothed rate map with
%appropriate bin range
%based on rate map
ST_curve = norm_rate_map_sm{bin_STC};

%non-normalized rate map (Danielson makes no mention of ROI-based rate map
%normalization
ST_curve_nonNorm = rate_map_sm{bin_STC};

%based on dF map (based on Dombeck)
ST_dF_curve=norm_dF_map_sm{bin_STC};

%TEMPORARY - CHANGE IN THE FUTURE - RZ

%% Split the onset map along the center bin for place field analysis code (place_field_finder)
%shifts the onset onset map 50 bins down 
%bin 1 is now at position 50
rate_map_center_shift = circshift(rate_map{8},50);

%add the center bin shifted map to the edges of the normal onset map
rate_map_extended = [rate_map_center_shift(1:50,:); rate_map{8}; rate_map_center_shift(51:100,:)];

%plot shifted onset map
figure;
imagesc(rate_map_extended');
hold on
title('Extended onset map');
%lap start
plot([51 51], [0 size(rate_map_extended,2)],'r');
%lap end
plot([151 151], [0 size(rate_map_extended,2)],'r');
hold off

%% Sort spatial tuning curve by max rate
%STC = spatial tuning curve

%based on norm smoothed curve - should preserve the maxima - RZ
%M - returns max (1) or NaNs
%I - bin where max was found
[M,I]=max(ST_curve); %find max in each column and save index

%make matrix with each column being the numeric index of ROI, bin where the normalized rate is max, and
%spatial tuning curve for each ROI (row) and each bin (column)
ord_STC=[(1:size(ST_curve,2))' I' ST_curve'];

%sort each row (ROI) according to the bin where it has max firing rate (column 2)
%1st column - sorted indices
%2nd column - bins in sorted order
%3-102 column - ROI x spatial bin - sorted tuning curve
ROI_ord_STC=sortrows(ord_STC,2);

%STC taken out from sort matrix (without ROI indices and bins)
STC_sorted=(ROI_ord_STC(:,3:end));

%Sort mean dF/F by max rate (dF_map is not a rate map)

%make matrix with each column being the numeric index of ROI, bin where the normalized rate is max, and
%spatial tuning curve for each ROI (row) and each bin (column)
ord_STC_dF=[(1:size(ST_dF_curve,2))' I' ST_dF_curve'];

%sort ROIs according to bin where max dF rate occurred - RZ
ROI_ord_STC_dF = sortrows(ord_STC_dF,2);

%take out the dF_map (sorted according to the rate map) - RZ
STC_dF_sorted=(ROI_ord_STC_dF(:,3:end));

%Remove NaN - find is any value is NaN 
%isnan - returns 1 where the index is NaN value
%any - find which values are non-zero along the chosen dimension
    %2- along each row (for each ROI across all bins)
%~take the indices which are non-NaN by taking inverse
%use as indices to STC_sorted and take all spatial bins
STC_sorted_nonan=STC_sorted(~any(isnan(STC_sorted),2),:);

%non-NaN logical map (1's is keep values)
ROI_sorted_nonNan = ~any(isnan(STC_sorted),2);

%how many NaN indices were found - RZ
nb_ROI_Nans = size(STC_sorted,1) - size(STC_sorted_nonan,1);

%take out column with ROI order and place in separate variable (contains
%NaNs
ROInb_STC_sorted = ROI_ord_STC(:,1);

%ROI ordered without ROIs that are NaNs
ROInb_STC_sorted_noNan = ROI_ord_STC(ROI_sorted_nonNan,1);

%% Spatial information - RZ - review both formulas and confirm

%for each bin range
for i=1:length(options.Nbin)
    %spatial information according to (Danielson et al. 2016)
    %smoothed rate map  - check if smoothing is appropriate for this
    %divided by the smoothed mean rate with NaN rates
    %mean rate for each ROI across all bins
    %what fraction of the rate occurs in each bin given total rate across
    %laps
    %this division should repmat mean_rate automatically - verify - RZ
    rate_ratio{i} = (rate_map_sm{i}./mean_rate_sm{i});
    
    %rate_ratio_noSm{i} = (rate_map{i}./mean_rate{i});
    rate_ratio_noSm{i} = (rate_map{i}./mean_rate{i});
    
    %spatial information metric:
    %log - returns natural log of each element in the array
    %take whole rate map matrix - element-wise multiply by the rate of relative rate of activity for that ROI in each bin
    %times how likely is it that the animal is in that bin
    %try constructing spatial info maps in the future - RZ
    spatial_info_bin{i}=rate_map_sm{i}.*reallog(rate_ratio{i}).*proba_bin{i}';
    
        
    %spatial_info_bin_noSm{i} = rate_map{i}.*reallog(rate_ratio_noSm{i}).*proba_bin{i}';
    %sum the spatial info scores across all the bins ,excluding NaN values
    %SI - spatial information score for each bin
    %vector for SI scores for each ROI for each bin range 
    %each column - spatial information score for each ROI
    %each row spatial information score for each bin range -RZ
    %spatial_information score increases as # of bins increases
    SI(i,:)=nansum(spatial_info_bin{i});
    
    %use non-smoothed values -RZ
    %SI(i,:)=nansum(spatial_info_bin_noSm{i});
    
    %(Skaggs) - all as function of relative activity in each bin and
    %occupancy of that bin
    spatial_info_bin_S{i}=rate_ratio{i}.*log2(rate_ratio{i}).*proba_bin{i}';
    
    %SIS - spatial information Skaggs
    SIS(i,:)=nansum(spatial_info_bin_S{i});
end

%% Save structure

%TODO - add spatial info calculated with nonsmoothed parameters

Place_cell.options=options;

%how likely that the animal is in that bin
Place_cell.proba_per_bin=proba_bin;
%spatial tuning curve - normalized and smoothed rate map
Place_cell.Spatial_tuning_curve=ST_curve;
%spatial tuning curve - non-normalized; only smoothed (smoothed event
%map)/unsmoothed occupancy time
Place_cell.Spatial_tuning_curve_no_norm=ST_curve_nonNorm;

%normalized and smoothed average dF values per bin
Place_cell.Spatial_tuning_dF=ST_dF_curve;
%sorted spatial tuning curve according to bin where max normalized STC is
Place_cell.Spatial_tuning_curve_sorted=STC_sorted;
%sorted dF tuning curve according to bin where max normalized STC is
Place_cell.Spatial_tuning_dF_sorted=STC_dF_sorted;
%ROI indices according to spatial tuning maximums
Place_cell.ROI_Spatial_tuning_sorted=ROInb_STC_sorted;

%how many Ca  events per bin for each bin range
Place_cell.Spatial_Info.event_map=onset_map;
%mean dF/F value in each bin (unsorted), non-normalized, non-smoothed
Place_cell.Spatial_Info.mean_dF_map=dF_map;
%mean dF/F value in each bin smoothed with moving average filter
Place_cell.Spatial_Info.mean_dF_map_smooth=dF_map_sm;

%non-smoothed occupancy
Place_cell.Spatial_Info.occupancy_map=occupancy_time_s;
%non-smoothed rate map
Place_cell.Spatial_Info.rate_map=rate_map;
%extended onset map for place field analysis
Place_cell.Spatial_Info.extended_rate_map = rate_map_extended;
%smoothed rate map (spatial info)
Place_cell.Spatial_Info.rate_map_smooth=rate_map_sm;
%smoothed rate map (plotting/sorting)
Place_cell.Spatial_Info.rate_map_sm_plot = rate_map_sm_plot;

%average rate for each ROI
Place_cell.Spatial_Info.overall_rate=mean_rate;
%probability that animal is in a given bin for each bin range 
Place_cell.Spatial_Info.proba_bin=proba_bin;
%spatial info score for each bin range (Danielson)
Place_cell.Spatial_Info.Spatial_Info=SI;
%spatial info score for each bin range (Skaggs)
Place_cell.Spatial_Info.Spatial_Info_Skaggs=SIS;
%binary event onsets for each running frame for each ROI 
Place_cell.Spatial_Info.Run_onset_bin=run_onset_bin;

% which bin each running frame belongs to for each bin range
Place_cell.Bin=bin;
% add the bin frames mask that say what bin each frames is in during run
% epochs
Place_cell.binFrameMask = binFrameMask;

% what position (in cm) each bin belongs to
%bin # and corresponding location (cm) for each bin range
Place_cell.Position_Bin=mean_pos_bin;

%Barthos place field finder input dF/F traces (100 bins with mean dF/F from only sig run events) 
Place_cell.dF_sig_map_100_smooth = dF_sig_map_100_smooth;

%number of frames in which sig calcium transients are occuring across 100
%spatial bins
%in frames
Place_cell.onset_map_conti_fr = onset_map_conti{8};
%in secodns
Place_cell.onset_map_conti_s = onset_map_conti{8}.*avg_fr;

%Spatialinfo.timebin=time_bin_sec;
%Spatialinfo.positionbin=position_bin;
%Spatialinfo.onset_bin=onset_bin;

end