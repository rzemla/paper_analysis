function [Place_cell] = place_field_finder_gaussian(Place_cell,options)

%find place fields based on Zaremba et al. 2017
%gaussian kernel with sigma = 3 spatial bins

%use rate map - number of event onsets/ occupancy across all laps
gSigma = options.gSigma;
%rename option variable to use for input for define_Gaussian_kernel fxn
options.sigma_filter = options.gSigma;

%if to exclude based on center of the located peak
centerExclude = options.centerExclude;

%peak distance separation (in bins)
%not used if centerExclude is set to 0
peakDistance = options.peakDistance;

%fraction of each peak to define as left and right edge of each Gaussian
%fit
peakFraction = options.peakFraction;

%how large the area has to be of other fields in comparison to the the
%place field with the largest area;
gaussAreaThreshold = options.gaussAreaThreshold;

%whether to plot the individual place fields
plotFields = options.plotFields;

%offset necessary from using extended rate map
offset =50;

%% Gaussian kernel smoothing filter
%using convolution with custom window
%size of window
gaussFilter = define_Gaussian_kernel(options);

%% Smooth extended_rate_map

ex_rate_map = Place_cell{1}.Spatial_Info.extended_rate_map;

%smooth extended rate map for each ROI
for rr = 1:size(ex_rate_map,2)
    ex_rate_map_sm(:,rr) = conv(ex_rate_map(:,rr),gaussFilter, 'same');
end

%% Find local maxima - (maxima must be between bin 51 and 150)

%[pks,lc] = findpeaks(ex_rate_map_sm(:,1),'MinPeakDistance', peakDistance);
%without min separation
%[pks,lc] = findpeaks(ex_rate_map_sm(:,1));

%find local maxima for each ROI
for rr=1:size(ex_rate_map,2)
    [pks{rr},lc{rr},w{rr}] = findpeaks(ex_rate_map_sm(:,rr),'WidthReference','halfprom',...
        'MinPeakHeight',0.1,'Threshold',0);
end


%% Discard peaks that are within (<=) 15 bins of edge of extended rate map
for rr=1:size(ex_rate_map,2)
    discard_pks = lc{rr} <= 15;
    lc{rr}(discard_pks) =[];
    pks{rr}(discard_pks) =[];
    w{rr}(discard_pks) =[];
end

%% Remove pks with width < 1.5 (~3 bins)
for rr=1:size(ex_rate_map,2)
    discard_pks = w{rr} <= 1.5;
    lc{rr}(discard_pks) =[];
    pks{rr}(discard_pks) =[];
    w{rr}(discard_pks) =[];
end

%% Exclude peaks that are beyond the boundaries of the complete laps (not in extention zone)
for rr=1:size(ex_rate_map,2)
    discard_pks = lc{rr} < 51 | lc{rr} > 150;
    lc{rr}(discard_pks) =[];
    pks{rr}(discard_pks) =[];
    w{rr}(discard_pks) =[];
end


%% Plot to visualize maxima - skip; move advanced plotting beloq
if 0
    figure;
    for rr=36%1:size(ex_rate_map,2)
        hold on;
        %plot smoothed rate map
        plot(ex_rate_map_sm(:,rr),'k');
        %plot start and end point of bin
        stem([51,150],[1 1],'b');
        %plot peaks
        stem(lc{rr},pks{rr},'r')
        %plot edges based on width return
        %start
        stem(lc{rr}-(w{rr}./2),ones(size(lc{rr},1),1),'g');
        %end
        stem(lc{rr}+(w{rr}./2),ones(size(lc{rr},1),1),'g');
        pause;
        clf;
    end
end

%% Select ROIS signficant tuned by SI score

SI_tuned_ROIs = find(Place_cell{1, 1}.Spatial_Info.significant_ROI == 1);
%do for all and then select tuned ones at the end

%% Fit single term gaussian to id'd local maxima

%number of bins to extend plotting of fit gauss beyond the intial width
%from findpeaks
gauss_extend = 10;
tic;
%for rr = SI_tuned_ROIs
%for all ROIs
for rr=1:size(pks,2)

    if ~isempty(pks{rr})
        
        %for each id'd peak
        for peak_nb =1:size(pks{rr})
            %x, y input
            %for each local maximum
            loc_range{rr}{peak_nb} = [round(lc{rr}(peak_nb)-(w{rr}(peak_nb)./2)):round(lc{rr}(peak_nb)+(w{rr}(peak_nb)./2))];
            curve_range{rr}{peak_nb} = ex_rate_map_sm(loc_range{rr}{peak_nb},rr);
            
            %fit gaussian to peak
            [f{rr}{peak_nb},gof{rr}{peak_nb},output{rr}{peak_nb}] = fit(loc_range{rr}{peak_nb}', ex_rate_map_sm(loc_range{rr}{peak_nb},rr),'gauss1');
            gauss_fit{rr}{peak_nb} = f{rr}{peak_nb}.a1*exp(-(([loc_range{rr}{peak_nb}]-f{rr}{peak_nb}.b1)./f{rr}{peak_nb}.c1).^2);
            gauss_fit{rr}{peak_nb} =f{rr}{peak_nb}.a1*exp(-(([loc_range{rr}{peak_nb}(1)-gauss_extend:loc_range{rr}{peak_nb}(end)+gauss_extend]-f{rr}{peak_nb}.b1)./f{rr}{peak_nb}.c1).^2);
            gauss_fit_bin_range{rr}{peak_nb} = [loc_range{rr}{peak_nb}(1)-gauss_extend:loc_range{rr}{peak_nb}(end)+gauss_extend];
        end

    end
end
toc;

%% Plot - split into 2 subplots with smoothed rate and non-smoothed rate
%skip for now - add as option later for display
if 0
    figure;
    for rr =542%SI_tuned_ROIs%1:100
        subplot(2,1,1)
        if ~isempty(pks{rr})
            %for each id'd peak
            for peak_nb =1:size(pks{rr})
                %plot range
                hold on
                title(num2str(rr));
                ylim([0 0.8])
                %plot extended rate map (smoothed)
                plot(ex_rate_map_sm(:,rr),'k')
                
                %plot extended rate map (non-smoothed)
                %plot(ex_rate_map(:,rr),'k-');
                %plot peak and width ends
                %plot peak center
                stem(lc{rr}(peak_nb),pks{rr}(peak_nb),'r')
                %plot width around peak
                %start
                stem(lc{rr}(peak_nb)-(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
                %end
                stem(lc{rr}(peak_nb)+(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
                %plot gaussian fit to ROI peaks
                plot([loc_range{rr}{peak_nb}(1)-gauss_extend:loc_range{rr}{peak_nb}(end)+gauss_extend],gauss_fit{rr}{peak_nb},'m')
            end
        end
        
        %plot start and end point of bin
        stem([51,150],[1 1],'b');
        %cutoff transient rate refline
        cutoff_line = refline(0,0.1);
        cutoff_line.Color = [0.5 0.5 0.5];
        cutoff_line.LineStyle = '--';
        %second cut off line
        cutoff_line2 = refline(0,0.05);
        cutoff_line2.Color = [0.5 0.5 0.5];
        cutoff_line2.LineStyle = '--';
        
        subplot(2,1,2)
        
        hold on
        if ~isempty(pks{rr})
            %for each id'd peak
            for peak_nb =1:size(pks{rr})
                %plot range
                hold on
                title(num2str(rr));
                ylim([0 0.8])
                %plot extended rate map (smoothed)
                %plot(ex_rate_map_sm(:,rr),'k')
                
                %plot extended rate map (non-smoothed)
                plot(ex_rate_map(:,rr),'k-');
                %plot peak and width ends
                %plot peak center
                stem(lc{rr}(peak_nb),pks{rr}(peak_nb),'r')
                %plot width around peak
                %start
                stem(lc{rr}(peak_nb)-(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
                %end
                stem(lc{rr}(peak_nb)+(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
                %plot gaussian fit to ROI peaks
                plot([loc_range{rr}{peak_nb}(1)-gauss_extend:loc_range{rr}{peak_nb}(end)+gauss_extend],gauss_fit{rr}{peak_nb},'m')
            end
        end
        
        %plot start and end point of bin
        stem([51,150],[1 1],'b');
        %cutoff transient rate refline
        cutoff_line = refline(0,0.1);
        cutoff_line.Color = [0.5 0.5 0.5];
        cutoff_line.LineStyle = '--';
        %pause
        %clf
    end
end


%% Check intersecting gaussian fit curves and merge

%add minimum smoothed transient rate for cross (set 0 for now adjust later)
%may change this on percentage of highest/lowest curve or decide this value
%empirically
minCurveCross = 0.05;

%check sequential peaks

%if there is a peak (gauss fit that spills beyond end lap, shift by bin lap length
%(100 bins) - do this separately

%threshold for edge of gaussian curve as to where end of field is
%fixed values of percentage of the height of the id'd peak - decide

%%% use ROI 17 as a starting point to test out merging
%check all neighboring peaks first
%then look for anyone with curve extending beyong start or end, if so shift
%and compare to peak next to it
rr=17; %single peak merger test
rr=36; %edge merger test (forward shift peak)
%rr = ;%(rearward shift peak)

%for all spatially tuned ROIs
%for rr =SI_tuned_ROIs
%for all ROIs
%preallocate merge_end to 0
merge_end = zeros(1,size(gauss_fit,2));

for rr=1:size(gauss_fit,2)
    %if more than 1 peak
    if size(gauss_fit{rr},2) > 1
        %for each peak comparison
        for pp=1:size(gauss_fit{rr},2)-1
            int_pt{rr}{pp} = InterX([gauss_fit_bin_range{rr}{pp};gauss_fit{rr}{pp}],[gauss_fit_bin_range{rr}{pp+1};gauss_fit{rr}{pp+1}]);
        end
        
        %create logical with intersection
        %for each comparison, check if there is an intersection and is above
        %minimum value
        for cc=1:size(int_pt{rr},2)
            if isempty(int_pt{rr}{cc})
                %no crossing of curves
                merge_middle{rr}(cc) = 0;
            else
                %if above minimuj threshold for crossing
                if int_pt{rr}{cc}(2) > minCurveCross
                    merge_middle{rr}(cc) = 1;
                else %if does not exceed threshold - no merge flag
                    merge_middle{rr}(cc) = 0;
                end
            end
            
        end
        
        %check if first or last peak crosses edge of lap bin (<51 or >150) -
        %only 1 test should be true when tested
        %first peak
        if sum(gauss_fit_bin_range{rr}{1} < 51) | sum(gauss_fit_bin_range{rr}{1} > 150)
            %shift forward by 100 bins and check intersection with last peak
            int_edge{rr} = InterX([gauss_fit_bin_range{rr}{1}+100;gauss_fit{rr}{1}],[gauss_fit_bin_range{rr}{end};gauss_fit{rr}{end}]);
            %assign 1 flag
            crossed_edge(rr) = 1;
            %if intersects,merge (update logical with intersection
            %last peak
        elseif (gauss_fit_bin_range{rr}{end} < 51) | sum(gauss_fit_bin_range{rr}{end} > 150)
            %shift backward by 100 bins and check intersection with first peak
            int_edge{rr} = InterX([gauss_fit_bin_range{rr}{end}-100;gauss_fit{rr}{end}],[gauss_fit_bin_range{rr}{1};gauss_fit{rr}{1}]);
            %assign -1 flag
            crossed_edge(rr) = -1;
        else
            %set to empty for conditional check below
            int_edge{rr} = [];
        end
        
        %check if there is intersection at edges and if yes, then set merge
        %flag
        if ~isempty(int_edge{rr})
            %check flag for which end curve is shifted
            if crossed_edge(rr) == 1
                %check if minimum for curve crossing is met
                if int_edge{rr}(2) > minCurveCross
                    %set merge to +1
                    merge_end(rr) = 1;
                end
            elseif crossed_edge(rr) == -1
                if int_edge{rr}(2) > minCurveCross
                    %set merge to -1
                    merge_end(rr) = -1;
                end
            end
        end
    end
end

%% Get endpoints of each gaussian based on either fraction of peak or set threshold 
%use fraction of peak for now; implement fixed threshold later

%fraction of peak to use as endpoint
fracPeak = 0.2;

for rr=1:size(gauss_fit,2)
    %for each peak if not empty
    if ~isempty(gauss_fit{rr})
        for pp=1:size(gauss_fit{rr},2)
        %gives the bin edge and value at thresgold
        gauss_fit_edges{rr}{pp} = InterX([gauss_fit_bin_range{rr}{pp}; ones(1,size(gauss_fit_bin_range{rr}{pp},2))*fracPeak*max(gauss_fit{rr}{pp})],...
            [gauss_fit_bin_range{rr}{pp};gauss_fit{rr}{pp}]);
        end
    end
end
%% Set merge input such into groups of which #'ed peaks to merge
%alternative is to re-fit a single term gaussian into the smoothed space
%now defined by the new endpoints - likely to be broader and overestimate
%field

%for middle peaks
%flag for whether previous peak was merged
prevMerge = 0;

for rr=1:size(merge_middle,2)
    
    %if more than 1 peak
    if size(gauss_fit{rr},2) > 1
        %for each merge comparison
        for mm=1:size(merge_middle{rr},2)
            if merge_middle{rr}(mm) == 0 && prevMerge == 0
                merge_peak_nb{rr}{mm} = mm;
                prevMerge = 0;
            elseif merge_middle{rr}(mm) == 1 && prevMerge == 0
                merge_peak_nb{rr}{mm} = [mm,mm+1];
                prevMerge = 1;
                %remove previous single peak if current merge
                if mm >1
                     if size(merge_peak_nb{rr}{mm-1},1) ==1 & merge_peak_nb{rr}{mm-1} == mm
                         %delete
                         merge_peak_nb{rr}{mm-1} = [];
                     end
                end
            elseif merge_middle{rr}(mm) == 1 && prevMerge == 1
                merge_peak_nb{rr}{end} = [merge_peak_nb{rr}{end},mm+1];
                prevMerge = 1;
           elseif merge_middle{rr}(mm) == 0 && prevMerge == 1 
                merge_peak_nb{rr}{mm} = mm+1;
                prevMerge = 0;
            end
        end
        %outside check - if last was no merge, add single peak to end
        if merge_middle{rr}(mm) == 0 && prevMerge == 0
            %if not already added
            if merge_peak_nb{rr}{end} ~= size(merge_middle{rr},2)+1
            merge_peak_nb{rr}{end+1} = size(merge_middle{rr},2)+1;
            end
        end
        %add last peak if standalone
    %if 1 peak
    elseif size(gauss_fit{rr},2) == 1
        merge_peak_nb{rr} = 1;
        %if no peaks
    elseif isempty(gauss_fit{rr})
        merge_peak_nb{rr} = [];
    end
    %reset prevMerge flag for next ROI
    prevMerge = 0;
end

%for edges
for rr=1:size(merge_end,2)
    if merge_end(rr) == 1 || merge_end(rr) == -1 %forward shift and integrate with end
        %designate mere peaks
        merge_peak_nb_edge{rr} = [size(gauss_fit{rr},2),1];
    else
        merge_peak_nb_edge{rr} = [];
    end
end

%make combined merging of middle and edge points into 1 cell
%take first and end point (1) (end)
%see if there are single or combined vectors
%if both single, delete, both, add shared
%if either 1 has more than 1 element; add to bigger 1, delete smaller
%if both more than one elements, merge into 1

for rr =1:size(merge_peak_nb_edge,2)
    %if there are peaks to merge
    if ~isempty(merge_peak_nb_edge{rr})
        if size(merge_peak_nb{rr}{end},2)==1 & size(merge_peak_nb{rr}{1},2)==1
            %set first index to merge
            merge_peak_nb{rr}{1} = merge_peak_nb_edge{rr};
            %delete last
            merge_peak_nb{rr}{end} = [];
        
        elseif size(merge_peak_nb{rr}{end},2)>1 & size(merge_peak_nb{rr}{1},2)>1
            %merge into first position
            merge_peak_nb{rr}{1} = [merge_peak_nb{rr}{end}, merge_peak_nb{rr}{1}];
            %delete last
            merge_peak_nb{rr}{end} = [];
        else %if either end contains standalone peak
            %first has multiple, take last and end to first, delete last
            if size(merge_peak_nb{rr}{1},2)>1
               merge_peak_nb{rr}{1} = [merge_peak_nb{rr}{end}, merge_peak_nb{rr}{1}];
               %delete last
               merge_peak_nb{rr}{end} = [];
            else size(merge_peak_nb{rr}{end},2)>1
               merge_peak_nb{rr}{end} = [merge_peak_nb{rr}{end},merge_peak_nb{rr}{1}];
               %delete first
               merge_peak_nb{rr}{1} = [];
            end
        end
    end
end

%% Set place field centers, edges

%check for peaks
for rr =1:size(merge_peak_nb_edge,2)
    if ~isempty(merge_peak_nb{4})
        %check if cell (more than 1 intial peaks)
        if iscell(merge_peak_nb)
        %single peak
        else
            
        end
    else
        placeField.edge{rr} = [];
        placeField.width{rr} = [];
        placeField.centers{rr} = [];
    end
    
end

iscell(merge_peak_nb{6})

%% Merge (take endpoints of fitting Gaussians) and set place field endpoints


%% Save to output struct



%both work well - use InterX for now
%one function
inter_point = InterX([gauss_fit_bin_range{17}{2};gauss_fit{17}{2}],[gauss_fit_bin_range{17}{3};gauss_fit{17}{3}])
%another intersection function
%[inter_point_2(1),inter_point_2(2)] = intersections(gauss_fit_bin_range{17}{2},gauss_fit{17}{2},gauss_fit_bin_range{17}{3},gauss_fit{17}{3})

rr =17;
figure;
if 1
%plot intersection point
    if ~isempty(pks{rr})
        %for each id'd peak
        for peak_nb =1:size(pks{rr})
            %plot range
            hold on
            title(num2str(rr));
            ylim([0 0.8])
            %plot extended rate map (smoothed)
            plot(ex_rate_map_sm(:,rr),'k')
            
            %plot extended rate map (non-smoothed)
            %plot(ex_rate_map(:,rr),'k-');
            %plot peak and width ends
            %plot peak center
            stem(lc{rr}(peak_nb),pks{rr}(peak_nb),'r')
            %plot width around peak
            %start
            stem(lc{rr}(peak_nb)-(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
            %end
            stem(lc{rr}(peak_nb)+(w{rr}(peak_nb)./2),pks{rr}(peak_nb)*ones(size(lc{rr}(peak_nb),1),1),'g');
            %plot gaussian fit to ROI peaks
            plot([loc_range{rr}{peak_nb}(1)-gauss_extend:loc_range{rr}{peak_nb}(end)+gauss_extend],gauss_fit{rr}{peak_nb},'m')
        end
    end
    
    %plot start and end point of bin
    stem([51,150],[1 1],'b');
    %cutoff transient rate refline
    cutoff_line = refline(0,0.1);
    cutoff_line.Color = [0.5 0.5 0.5];
    cutoff_line.LineStyle = '--';
    %plot test intersection point
    stem(inter_point(1),inter_point(2),'r--','LineWidth', 2)
end

%want to know at what value the curves intersect - try function
%deal with merging of curves that 'go around' track
%


%%% OLD CODE BELOW HERE %%% 
%%%%% DELETE ONCE DONE WITH NEW CODE ABOVE %%%%%

%%
%{
for ii=1:size(lc,1)
    %make sure that the peak is not nan
    if ~isnan(pks(ii))
        %fit Gaussian into each peak
        [f{ii},gof{ii},output{ii} ] = fit([edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)]', convG(edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)),'gauss1');
        %gaussian curve from fit (need x inputs)
        gauss_fit{ii} = f{ii}.a1*exp(-(([edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)]-f{ii}.b1)./f{ii}.c1).^2);
        %get area (bins)
        gauss_area(ii) = trapz(gauss_fit{ii});
    else
        f{ii} = nan;
        gof{ii} = nan;
        output{ii} = nan;
        gauss_fit{ii} = nan;
        gauss_area(ii) = nan;
        
    end
    
end


%% Debugging info
%work on 101
%where merging is performed [36 55 94 101 144, 161, 217 233]

%problems - 233, 217, 208, 199, 182, 174, 161, 144, 137, 101, 94, 89, 55,
% 69, 99 -merge example
% 78 - mult peak example
% 69 - merge and multi
%  64, 59, 299, 298, 294, 285, 281, 269, 215 - multi

%TODO - DONE
%make center of Gaussian fit, the center of the PF
%not necessary as the peak points towards the higher inital value after
%smoothing
%%
%vector for holding number of PFs for given ROI
placeFieldNb = zeros(size(Place_cell{1}.Tuned_ROI_mask,2),size(Place_cell,2));

%for each trial type
for tt=1:size(Place_cell,2)

    %find all spatially tuned ROI (by spatial info criterion)
    sigROI = find(Place_cell{tt}.Spatial_Info.significant_ROI==1);
    
    %number of ROIS
    ROInb = size(Place_cell{tt}.Spatial_Info.significant_ROI,2);
    %make blank cell container for fieldWidths and centerLoc (prevent skipped cells)
    fieldWidth{tt} = cell(1,ROInb);
    centerLoc{tt} = cell(1,ROInb);
    
    %for each neuron selected according to spatial info criterion
    for rr=1:size(sigROI,2)
        
        ROI = sigROI(rr);
        
        disp('ROI: ')
        disp(ROI)
        %normalized map
        %rate_map_norm_smooth = Place_cell{1}.Spatial_tuning_curve(:,ROI);
        
        %non-normalized map (from 100 bins)
        %rate_map = Place_cell{tt}.Spatial_Info.rate_map{8}(:,ROI);
        
        %extended map
        rate_map = Place_cell{tt}.Spatial_Info.extended_rate_map(:,ROI);
        
        %convolve the rate map with a Gaussian filter
        convG = conv(rate_map,gaussFilter, 'same');
        
        %both filters match each other around window size 27 since imgaussfilt
        %using window size on each side defined by the equation above
        
        clear pks lc
        %find local maxima
        %restrict the search to between the start and end bin
        
        [pks,lc] = findpeaks(convG(51:150),'MinPeakDistance', peakDistance);
        
        %conditional for when the peak lies at the very edge of lap
        %add 10 to the edge
        if isempty(pks)
            [pks,lc] = findpeaks(convG(41:160),'MinPeakDistance', peakDistance);
            pks = pks(1);
            lc = lc(1) - 10;
        end
        
        %add offset to the peak location value given the search
        %restriction
        lc = lc + 50;
        
        %make a copy of the original peaks and locations
        pks_original = pks;
        lc_original  = lc;
        
        %% Exclude ROIs that are location within start or end 10 bins - NOT used with extended map
        %bypass exclusion and exclude ROIs with edges that are not defined
        %make sure that centerExclude is set to 0
%         
%         if centerExclude ==1
%             nearEdgeIdx = find(lc < 10 | lc > 90);
%             %set them to nan
%             if ~isempty(nearEdgeIdx)
%                 pks(1:size(pks,1)) = nan;
%                 lc(1:size(pks,1)) = nan;
%             end
%         end
        
        %% find edge threshold of each peak
        %calculate the threshold for the onset/offset of each field on each side
        peakThres = pks*peakFraction;
        
        clear edgeIdx
        %edgeIdx = [];
        %find left and right edge for each peak - works
        for ii=1:size(lc,1)
            %each row is id'd peak
            %check that peak is not at the edge of the track
            if ~isnan(pks(ii))
                %left edge
                if ~isempty(find(convG(fliplr(1:lc(ii))) < peakThres(ii),1))
                    edgeIdx(ii,1) = lc(ii) - find(convG(fliplr(1:lc(ii))) < peakThres(ii),1) +1 ;
                else
                    edgeIdx(ii,1) = nan;
                    pks(ii) = nan;
                    %lc(ii) = nan;
                end
                %right edge
                if ~isempty(find(convG(lc(ii):end) < peakThres(ii),1))
                    edgeIdx(ii,2) = lc(ii) + find(convG(lc(ii):end) < peakThres(ii),1) -1 ;
                else
                    edgeIdx(ii,2) = nan;
                    pks(ii) = nan;
                    %lc(ii) = nan;
                end
            else
                edgeIdx(ii,1) = nan;
                edgeIdx(ii,2) = nan;
            end

        end
        
        %% Merge peaks that overlap
        %look at the edges and check if neighboring edges run past each other
        
        %turn this into iterative process - while
        %flag for merging 
        merging = 1;
        mergeEdge = 0;
        mergePks = [];
        
        %while merging ==1 
        %if more than one peak (PF)
        if size(edgeIdx,1) > 1
            
            %for all peaks
            for ii =1:size(edgeIdx,1)
                %check that edges are not nan
                if ~isnan(edgeIdx(ii,1))
                    
                    %if first peak - check if right edge exceeds center of next peak
                    if ii == 1
                        if (edgeIdx(ii,2) > lc(ii+1))
                            %marker that right edge of first peak exceeds
                            %center of next peak
                            mergePks(ii,:) = [0 1];
                            merging = 1;
                            
                        else
                            mergePks(ii,:) = [0 0];
                            merging = 0;
                        end
                        
                        %if last peak - check if left edge is below center of peak before
                    elseif ii == size(edgeIdx,1)
                        if edgeIdx(ii,1) < lc(ii-1)
                            mergePks(ii,:) = [1 0];
                            merging = 1;
                        else
                            mergePks(ii,:) = [0 0];
                            merging = 0;
                        end
                        
                    else %- check if left edge falls below previous peak or right edge
                        %falls above next  peak
                        %left edge below previous peak
                        if edgeIdx(ii,1) < lc(ii-1)
                            mergePks(ii,:) = [1 0];
                            merging = 1;
                            %right edge above next peak
                        elseif  (edgeIdx(ii,2) > lc(ii+1))
                            mergePks(ii,:) = [0 1];
                            merging = 1;
                        else
                            mergePks(ii,:) = [0 0];
                            merging = 0;
                        end
                    end
                end
                
            end
            
        end
        %end of merge process
        
        %end of while loop
        %end
        
        %save merge matrix
        mergeRR{rr} = mergePks;
        
        %set merge flag to 0
        merged = 0;
        
        %find if any row (corresponding to peak) has a value of 1
        for ii=1:size(mergePks,1)
            [~,cTemp] = find(mergePks(ii,:) == 1);
            if ~isempty(cTemp)
                %if right edge flagged
                if cTemp == 2
                    %check if left edge of next peak is flagged
                    if mergePks(ii+1,1) ~= 1
                        %remove the peak if neighbor not flagged
                        pks(ii) = nan;
                        lc(ii) = nan;
                        mergeEdge = 1;
                    elseif mergePks(ii+1,1) == 1
                        %merge with previous one depending on which one is
                        %larger
                        %set merge flag
                        merged = 1;
                        if pks(ii+1) > pks(ii)
                            %remove (set to nan) location and associated peak
                            lc(ii) = nan;
                            pks(ii) = nan;
                            %change the right edge of the next peak
                            edgeIdx(ii+1,1) = edgeIdx(ii,1);
                        else
                            %remove (set to nan) location and associated peak
                            lc(ii+1) = nan;
                            pks(ii+1) = nan;
                            %change the left edge of the next peak
                            edgeIdx(ii,2) = edgeIdx(ii+1,2);
                        end
                    end
                    
                elseif cTemp == 1 && merged ==0
                    %check if left edge of next peak is flagged
                    if mergePks(ii-1,2) ~= 1
                        %remove the peak if neighbor not flagged
                        pks(ii) = nan;
                        lc(ii) = nan;
                        %set merge edge flag
                        
                    elseif mergePks(ii-1,2) == 1
                        %merge with previous one depending on which one is
                        %larger
                        %set merge flag
                        merged = 1;
                        if pks(ii-1) > pks(ii)
                            %remove (set to nan) location and associated peak
                            lc(ii) = nan;
                            pks(ii) = nan;
                            %change the right edge of the next peak
                            edgeIdx(ii-1,2) = edgeIdx(ii,2);
                        else
                            %remove (set to nan) location and associated peak
                            lc(ii-1) = nan;
                            pks(ii-1) = nan;
                            %change the left edge of the previous peak
                            edgeIdx(ii,2) = edgeIdx(ii-1,1);
                        end
                    end
                end
            end
        end
        
        %reshape(mergeRR{1, 55},1,[])
        clear mergePks
        
        %check if either end lap peaks have been merged together
        
        %make copy of edges for curveFit
        edgeIdx_curvefit = edgeIdx;
        
        %subtract offset to the edges
        edgeIdx_offset = edgeIdx - offset;
        
        %check if both end edge is less than 1 or greater than 100
        [r_edge,~] = find(edgeIdx_offset(1,1) < 1 && edgeIdx_offset(end,2) > 100);
        
        %if overlapping, use the widest peak and remove the less wide
        edgeIdx_diff = abs(edgeIdx(:,2)- edgeIdx(:,1));
        
        %perform merging if both sides overlap
        if ~isempty(r_edge)
            if edgeIdx_diff(1)> edgeIdx_diff(end)
                lc(end) = nan;
                pks(end) = nan;
                edgeIdx(1,1) = edgeIdx(end,1);
                %edges for Gaussian fitting
                edgeIdx_curvefit(1,1) = edgeIdx(end,1);
                edgeIdx_curvefit(1,2) = edgeIdx(1,2) + 100;
                
            elseif edgeIdx_diff(1)< edgeIdx_diff(end)
                lc(1) = nan;
                pks(1) = nan;
                edgeIdx(end,2) = edgeIdx(1,2);
                %edges for Gaussian fitting
                edgeIdx_curvefit(end,2) = edgeIdx(1,2) + 100;
            end
            %check to see if any side if hanging over and addjust edge to
            %minimum
        elseif (edgeIdx_offset(1,1) < 1) || (edgeIdx_offset(end,2) > 100)
            %check which edge peak is out of range
            if (edgeIdx_offset(1,1) < 1)
                %find the lowest value from max of previous peak until max
                %of current peak
                %[~, minIdx ] = min(convG(lc(end):(lc(1)+100)));
                [~, minIdx ] = find(fliplr(convG(lc_original(end):150)') < peakThres(1),1);
                %if can't find the min below threshold, take the lowest
                %value on that interval
                if isempty(minIdx)
                    [~, minIdx ] = min(fliplr(convG(lc_original(end):150)'));
                end          
                
                edgeIdx(1,1) = 150 - minIdx;
                edgeIdx_curvefit(1,1) = edgeIdx(1,1);
                edgeIdx_curvefit(1,2) = edgeIdx(1,2)+100;
                
            elseif (edgeIdx_offset(end,2) > 100) 
                [~, minIdx ] = find(convG(51:lc_original(1)) < peakThres(end),1);
                if isempty(minIdx)
                    [~, minIdx ] = min(convG(51:lc_original(1)));
                end
                
                edgeIdx(end,2) = minIdx + 50;
                edgeIdx_curvefit(end,2) = edgeIdx(end,2)+100;
            end
        end

        %% fit local maximum in rate map with a Gaussian and get area of Gaussian
        %iterate this and fit over every peak of neuron
        
        clear f gof output gauss_fit gauss_area
        
        for ii=1:size(lc,1)
            %make sure that the peak is not nan
            if ~isnan(pks(ii))
                %fit Gaussian into each peak
                [f{ii},gof{ii},output{ii} ] = fit([edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)]', convG(edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)),'gauss1');
                %gaussian curve from fit (need x inputs)
                gauss_fit{ii} = f{ii}.a1*exp(-(([edgeIdx_curvefit(ii,1):edgeIdx_curvefit(ii,2)]-f{ii}.b1)./f{ii}.c1).^2);
                %get area (bins)
                gauss_area(ii) = trapz(gauss_fit{ii});
            else
                f{ii} = nan;
                gof{ii} = nan;
                output{ii} = nan;
                gauss_fit{ii} = nan;
                gauss_area(ii) = nan;
                
            end
            
        end
        
        %% find the peak of each id'd gaussian and make that the place field center
        clear lc_max lc_corr
        %make copy for output
        lc_corr = lc_original;
                
        for gg=1:size(gauss_fit,2)
            if ~isnan(gauss_fit{gg})
                [~,lc_max(gg)] = max(gauss_fit{gg});
                %add offset
                lc_max(gg) = lc_max(gg)+edgeIdx_curvefit(gg,1) -1;
%                 disp('max')
%                 disp(lc_max);
%                 disp('orig')
%                 disp(lc_original);
                %add offset to original detected peak
                %add additional offset that is introduced for Gaussian
                %fitting
                if lc_max(gg) < 151
                    lc_corr(gg) = lc_max(gg);
                else
                    lc_corr(gg) = lc_max(gg) - 100;
                end
            else
                lc_max(gg) = nan;
            end
        end

        %% discard any peaks with area less than 50% of largest Gaussian/peak by area
        
        %find max area peak for given ROI
        %do not search if all gaussian areas are nan
        if ~(sum(isnan(gauss_area)) == size(gauss_area,2))
            [maxArea,idxMaxArea] = max(gauss_area);
        end
        
        %divide all areas by max area
        gaussAreaNorm = gauss_area/maxArea;
        
        %find indices of fields which are less than the area of the max - threshold
        %defined in gaussAreaThreshold
        belowCriteriaPeaks = find(gaussAreaNorm < gaussAreaThreshold);
        
        lc_filtered = [];
        %if peaks found
        if ~isempty(belowCriteriaPeaks)
            %set them to NaN
            lc_filtered = lc;
            lc_filtered(belowCriteriaPeaks) = nan;
            %shifted centers
            lc_max(belowCriteriaPeaks) = nan;
        else
            lc_filtered = lc;
        end
        
        %exclude peaks that have undefined edges
        noEdgePeaks = logical(sum(isnan(edgeIdx),2));
        lc_filtered(noEdgePeaks) = nan;
        
        %exclude peaks if one of the edge peaks with undefined edge is greater
        %(avoid bias by exclusion due to lack of edge detection
        
        %for each accepted peak, remove if another exists that is larger
        maxPeak = max(pks_original);
        
        %TEMPORARY SKIP
%         for ii=1:size(lc_filtered,1)
%             if ~isnan(lc_filtered(ii))
%                 %make the peak is at least 2x the size either of the excluded
%                 if  (sum(noEdgePeaks) ~=0) %( pks_original(ii) < maxPeak)  ((pks_original(ii) < maxPeak) &&
%                     
%                     if (pks_original(ii) < 2*max(pks_original(noEdgePeaks)))
%                         lc_filtered(ii) = nan;
%                     end
%                 end
%             end
%         end
        

        %% Plot edges of id'd place fields - check for code optimization
        
        if plotFields == 1
            
        figure;
        hold on;
        
        %display ROI #:
        title(['ROI #: ', num2str(ROI)]);
        
        %rate map
        rm = plot(rate_map(51:150),'b');
        
        %edges of each peak
        for ee=1:size(lc_max,2)
            if ~isnan(lc_max(ee))
                plot([edgeIdx(ee,1)-offset, edgeIdx(ee,1)-offset], [0 1], 'r');
                plot([edgeIdx(ee,2)-offset, edgeIdx(ee,2)-offset], [0 1], 'r');
            end
        end
        
        %rate map convolved with gaussian smoothing kernel
        conG = plot(convG(51:150), 'k');
        
        %plot the identified local maxima
        stem(lc_max-offset, ones(size(lc,1)),'k', 'LineWidth', 1.5)
        
        %add stars over those that are significant (meet the area criteria)
        
        for ii=1:size(lc_max,2)
            %if filter peaks locations are not empty
            if ~isempty(lc_max)
                if ~isnan(lc_max(ii))
                    scatter(lc_max(ii)-offset, 1, 'r*');
                end
            end
        end
        
        %fitted Gaussian curve
        %fit curve - plot each fitted gaussian
        for ii=1:size(edgeIdx,1)
            %
            plot(edgeIdx_curvefit(ii,1)-offset:edgeIdx_curvefit(ii,2)-offset,gauss_fit{ii}, 'LineWidth',2)
        end
        
        hold off
        legend([rm conG], 'Rate map', 'Convolved gaussian smoothing filter');

        end
        
        
        %% save 
        %number of place fields per neuron
        placeFieldNb(ROI,tt) = size(lc_max,2) - sum(isnan(lc_max));
        
        %field width 
        fieldWidth{tt}{ROI} = edgeIdx_diff(~isnan(lc_max))*2;
        
        %center location of each PF
        centerLoc{tt}{ROI} = lc_corr(~isnan(lc_max))-offset;
        
        %% clear variables
        
        clear edgeIdx edgeIdx_curvefit edgeIdx_diff edgeIdx_offset
        %% Plot results
%         
%         if plotFields == 1
%             %plot different filters with raw data
%             figure;
%             hold on;
%             
%             %display ROI #:
%             title(['ROI #: ', num2str(ROI)]);
%             
%             %rate map
%             rm = plot(rate_map(51:150),'b');
%             
%             %rate map convolved with gaussian smoothing kernel
%             conG = plot(convG(51:150), 'k');
%             
%             %plot the identified local maxima
%             stem(lc-50, ones(size(lc,1)),'k', 'LineWidth', 1.5)
%             
%             %add stars over those that are significant (meet the area criteria)
%             
%             for ii=1:size(lc,1)
%                 %if filter peaks locations are not empty
%                 if ~isempty(lc_filtered)
%                     if ~isnan(lc_filtered(ii))
%                         scatter(lc_filtered(ii)-50, 1, 'r*');
%                     end
%                 end
%             end
%             
%             %fitted Gaussian curve
%             %fit curve - plot each fitted gaussian
%             for ii=1:size(edgeIdx,1)
%                 %
%                 plot(edgeIdx(ii,1)-50:edgeIdx(ii,2)-50,gauss_fit{ii}, 'LineWidth',2)
%             end
%             
%             hold off
%             legend([rm conG], 'Rate map', 'Convolved gaussian smoothing filter');
%             
%         end
        
    end
    
    %% save to Place_cell structure
    Place_cell{tt}.Field.placeFieldNb = placeFieldNb(:,tt);
    Place_cell{tt}.Field.centerLoc = centerLoc{tt};
    Place_cell{tt}.Field.width = fieldWidth{tt};
end

%%old code
%define the Gaussian filter to convolve with
%using 2D gauss imaging filter
%default window size = 2*ceil(2*sigma)+1 --> 13 for sigma = 3;
%windowSize = 2*ceil(2*gSigma)+1;
%yImGf = imgaussfilt(rate_map,gSigma);d

%}

end
