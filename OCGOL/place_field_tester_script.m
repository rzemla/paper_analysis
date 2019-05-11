%place field tester script

%% 1) Use this approach b/c based on rate map and tuning criteria in use are based on rate map
%Zaremba 2017
% For all cells, rate maps were formed by dividing the number of transients
% initiated in each spatial bin by the occupancy of that bin.
% We calculated rate maps with 100 position bins and smoothed 
% with a Gaussian kernel (? = 3 bins). To define place fields for 
% cells that were identified as containing significant spatial information,
% we fit each local maximum in the rate map with a Gaussian, merged 
% overlapping putative fields and then discarded any with an area less than 
% 50% of the largest.

%take extended rate map and smooth with gaussian kernel (sigma =3)
%define the start and end point in the lap bins (51 and 150) (1-50) -
%repetition before; (151-200) - repetition after
%find maxima 
%only run on spatially tuned ROIS by TS for jow

%use rate map - number of event onsets/ occupancy across all laps
options.gSigma = 3;

%if to exclude based on center of the located peak
options.centerExclude = 0;

%peak distance separation (in bins)
%not used if centerExclude is set to 0
options.peakDistance = 9;

%fraction of each peak to define as left and right edge of each Gaussian
%fit
options.peakFraction = 0.2;

%how large the area has to be of other fields in comparison to the the
%place field with the largest area;
options.gaussAreaThreshold = 0.3;

%whether to plot the individual place fields
options.plotFields = 0;

%Problem ROIs:
%B: 26, 42, 57, 69, 71, 101, 159,161,167
%A: 42, 132, 160
%dbstop if error
[Place_cell] = place_field_finder_gaussian(Place_cell,options);


%% 2 Dombeck/Barthos - don't use unless working in dF/F space
%Place field finder - Dombeck/Harvey 2010
tic;
%imaging period
options.dt = Imaging.dt;
%number of shuffles
options.nbShuffle = 1000;

%detect place fields and return - real
for ii=1:3
    [Place_cell{ii}] =  placeFieldDetection_barthosBeta(Place_cell{ii});
    
    
    if 0
        %shuffled the dF/F traces
        %1 min for 1000 shuffles of indices
        trace_shuffled_mat = shuffle_dFF(Imaging_split{ii},options);
        
        %recalculate the inputs parameters for each shuffle run
        %preallocate shuffle matrix with # of place fields per ROI
        placeField_nb_shuffle = zeros(options.nbShuffle,171);
        PF_inputs = cell(1,options.nbShuffle);
        
        tic;
        parfor n=1:options.nbShuffle
            %generate inputs based on shuffle
            PF_inputs{n} = activity_map_generator_PF(Behavior_split{ii}, Events_split{ii}, trace_shuffled_mat(:,:,n),options);
            %detect place fields and return - shuffled
            [placeField_nb_shuffle(n,:)] =  placeFieldDetection_barthosBeta_shuffle(PF_inputs{n});
        end
        toc;
        
        %get p_value for each ROI
        pVal = sum(logical(placeField_nb_shuffle))/options.nbShuffle;
        
        %assign which ROIs have sig fields
        Place_cell{ii}.placeField.sig_ROI = pVal < 0.05;
    end
    
end
toc;


%% Plot place fields for ROI on top of mean dF/F and transient rate

SI_sig_all = find(Place_cell{1, 3}.Spatial_Info.significant_ROI  ==1);

figure
for ii=1:100
    ROI = randi(418,1,1);
    ROI = SI_sig_all(ROI);
    hold on
    ylim([-0.2 2])
    %mean dF/F
    plot(Place_cell{1, 3}.Spatial_Info.mean_dF_map{1, 8}(:,ROI),'k');
    %rate map
    plot(Place_cell{1, 3}.Spatial_Info.rate_map{1, 8}(:,ROI),'r');
    
    
    %mark start and end of each place field for that neuron
    for pp=1:size(Place_cell{3}.placeField.range{ROI},1)
        stem([Place_cell{3}.placeField.range{ROI}(pp,1),Place_cell{3}.placeField.range{ROI}(pp,2)], [1,1]);
    end
    pause;
    clf;
end
