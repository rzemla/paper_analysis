
function [Place_cell]=tuning_specificity_RZ_V2(Place_cell,Behavior,Events,options)


%% Extract the data
%position of the animal in frames when animal was running
run_position=Behavior.resampled.run_position;

%bin position data
%N bin is the bin assignment to each frame
if options.binPosition == 1
    [Ncount, edges, Nbin] =histcounts(run_position,200);
    run_position = Nbin;
else
    %position of the animal in frames when animal was running
    run_position=Behavior.resampled.run_position;
end

%time of during the frames when the animal was running
run_time=Behavior.resampled.run_time;

%not used in the script - RZ - useful later when breaking down by lap
if Events.options.restrict==1
    lap_start_stop= Behavior.restricted.lap;
elseif Events.options.restrict==0
    lap_start_stop= Behavior.lap;
end

%events (onset) that occured when the animal was running for each ROI across all
%frames
%only event onsets that occurred when the animal was running, 
%check by seeing if the events fall into running epochs - RZ
run_onset_binary=Events.Run.run_onset_binary;

%vector of 1's when animal met running criteria
run_ones=Behavior.run_ones;

%ANDing EACH ROI AGAINST EPOCH RUN VECTOR SHOULD YIELD THE SAME VECTOR
%AS THE EVENT ONSET OR SUM OF EXTRACTED VECTOR SHOULD BE THE SAME
%try for one cell - WORKS RZ
% x = sum(run_onset_binary(logical(run_ones),100))
% y = sum(run_onset_binary(:,100));
% isequal(x,y);

%% Assemble DATA

%Take time and position for running epochs ONLY and normalize and against min and max position (0-1):
run_position_norm=(((run_position - min(run_position)) / (( max(run_position) - min(run_position)))));

%Take out the events that occured during the running epochs and project onto
%frames that occured during the running 
run_onset_bin=run_onset_binary(run_ones==1,:);

%same as the following expression
%run_onset_bin_2=run_onset_binary(logical(run_ones),:);
%isequal(run_onset_bin,run_onset_bin_2);

%run_onset_bin=run_onset_binary;

%position of the mouse at onset time:
%add onset lap and trial type to this - RZ
%for each ROI
for u=1:size(run_onset_bin,2)
    %find all the positions on the track when the event onset occurred
    %normalized run position == run_onset_bin in length
    onset_position{u}=run_position_norm(run_onset_bin(:,u)==1);
end

%TRY COMPUTING THIS ON DIFFERENT POSITION BINS (approximation to 5cm, 2cm, 1cm,
%and 0.5 cm rather than exact position - see how the tuning vectors compare


%Fraction of frames acquired at position of onset
%for each ROI
for u=1:size(onset_position,2)
    %for each position where an onset occurred
    
    %RZ bug fix - if not position because of no running onsets
    if size(onset_position{u},1) == 0
        frame_onset{u} = [];
        nbframe_onset{u} = [];
        fraction_frames{u} = [];
        
    else
        %for each position of onset
    for uu=1:size(onset_position{u},1)
        %make sure that the cell fired at at least 1 position
        %find how many frames were captured at that location
        if isempty(onset_position{u})==0
            %for every position of onset, give frame indices where this
            %EXACT POSITION OCCURRED AGAIN - RZ - I think this is too exact
            %- should perform this on binned locations (200 binned - every
            %1cm/ 0.005 of normalized position or 400 binned every 0.5 cm- check Danielson 
            %example: for every onset event at that location, list all frame indices that were captured at that EXACT location 
            frame_onset{u}{uu}=find(run_position_norm==onset_position{u}(uu));
            %how many frames occurred across all laps at the location of each transient
            %onset
            nbframe_onset{u}(:,uu)= size(frame_onset{u}{uu},1);
            %fraction of frames for each event onset divided by total run frames
            fraction_frames{u}(:,uu)= nbframe_onset{u}(:,uu)/ (size(run_onset_bin,1));
        end
    end
    end
end

%% Compute spatial tuning vectors

%working with normalized position, so plugging it into e multiplied by 2pi
%will give the relative position along the unit circle

%for each ROI
for u=1:size(onset_position,2)
    
    %if no onsets
    if size(onset_position{u},1) == 0
        %disp(size(onset_position{u},1))
        spatial_tuning_vector{u}=NaN;
        spatial_tuning_vector_angle{u}=NaN;
        spatial_tuning_vector_magnitude{u}=NaN;
    else
        
        %for each running related transient of that ROI
        for uu=1:size(onset_position{u},1)
            %make sure that at least one frame was found for that onset
            if isempty(frame_onset{u})==0
                %get the vector associated with position and divide by the
                %occupancy of that position
                %this is correct way to divide complex by real #
                %get a occupancy normalized tuning vector for each
                %running-related transient
                spatial_tuning_vector{u}(uu)=(exp(i*onset_position{u}(uu)*2*pi))./(fraction_frames{u}(:,uu));
                %returns the angle of each complex vector associated with each
                %running-related transient between +pi and -pi 
                spatial_tuning_vector_angle{u}(uu)=angle(spatial_tuning_vector{u}(uu));
                
                %magnitude of each running-related vector
                spatial_tuning_vector_magnitude{u}(uu)=abs(spatial_tuning_vector{u}(uu));
                
                               
                %fill all values with NaNs if no frame onset for that ROI
            elseif isempty(frame_onset{u})
                
                spatial_tuning_vector{u}(uu)=NaN;
                spatial_tuning_vector_angle{u}(uu)=NaN;
                spatial_tuning_vector_magnitude{u}(uu)=NaN;
            end
            
        end
        
    end
    
end


%Normalize magnitude of tuning vector
%for each ROI
for u=1:size(spatial_tuning_vector,2)
    %normalzied spatial tuning vector
    %the max of each of the normalized tuning vectors (large values because
    %of very low occupancy values)
    %normalize each magnitude by the maximum magnitude
    %vector of magnitudes divided by the max magnitude among them
    spatial_tuning_vector_magnitude_norm{u}=spatial_tuning_vector_magnitude{u}./max(abs(spatial_tuning_vector_magnitude{u}(:)));
    
    %take normalized magnitude for each running transient and multiply by
    %associated angle (which is already from -pi to pi (covers entire
    %distance)
    spatial_tuning_vector_norm{u} = spatial_tuning_vector_magnitude_norm{u}.*exp(spatial_tuning_vector_angle{u}*sqrt(-1));
    
    if isempty(spatial_tuning_vector_norm{u});
        spatial_tuning_vector_norm{u}=NaN;
    end
end

%% Compute Tuning specificity - one vector (sum of other vectors associated with events for given ROI)

%for each ROI
for u=1:size(spatial_tuning_vector_norm,2)
    
    % compute weighted sum of cos and sin of angles
    sum_vector(u) = sum(spatial_tuning_vector_magnitude_norm{u}.*exp(1i*spatial_tuning_vector_angle{u}));
    
    % compute magnitude = tuning specificity
    tuning_specificity(u) = abs(sum_vector(u))./sum(spatial_tuning_vector_magnitude_norm{u});
    
    %compute angle
    sum_angle(u)= angle(sum_vector(u));
    
    % get vector from magnitude and angle
    vector_tuning_specificity(u)=tuning_specificity(u).*exp(sum_angle(u)*sqrt(-1));
end

%RZ - should be same to 1-circular variance of occupancy normalized
%transients - verified - gives exactly the same result
%need circular statistics toolbox for this function
for rr=1:size(spatial_tuning_vector_angle,2)
    tuning_specificity_circ_var(rr) = 1-circ_var(spatial_tuning_vector_angle{rr},1./fraction_frames{rr},[],2);
end

%% write results in structure

Place_cell.Tuning_Specificity.tuning_vector=spatial_tuning_vector_norm;
Place_cell.Tuning_Specificity.tuning_specificity=tuning_specificity;
Place_cell.Tuning_Specificity.tuning_vector_specificity=vector_tuning_specificity;
Place_cell.Tuning_Specificity.tuning_specificity_circ_var = tuning_specificity_circ_var;

end