function [outputArg1,outputArg2] = transient_rate_in_field(session_vars)
%input - cell of structs with animal data


%% Define input variables 

%try for 1 session first

%for each session
%for ss=1:size(session_vars,2)

%Place field edge data
placeField_edges = session_vars{1}.Place_cell{4}.placeField.edge;

%Event rate for 100 bins
event_rate = session_vars{1}.Place_cell{4}.Spatial_Info.rate_map{8};

%event map (not rate)% 100 bins
event_map = session_vars{1}.Place_cell{4}.Spatial_Info.event_map{8}; 

%Occupancy for 100 bins
occupancy = session_vars{1}.Place_cell{4}.Spatial_Info.occupancy_map{8};

%should equal event rate (it does)
rate_map = event_map./occupancy';
%isequal(rate_map,event_rate)

%% Find and store transient rate in each place field

%for each ROI
for rr=1:size(event_rate,2)
    %for each place field
    if ~isempty(placeField_edges{rr})
        for pp=1:size(placeField_edges{rr},1)
            %if place bins split at lap edge
            if placeField_edges{rr}(pp,1) < placeField_edges{rr}(pp,2)
                place_bins{rr}{pp} = placeField_edges{rr}(pp,1):placeField_edges{rr}(pp,2);
            else
                place_bins{rr}{pp} = [placeField_edges{rr}(pp,1):100,1:placeField_edges{rr}(pp,2)];
            end            
            
            field_event_rates{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr))./sum(occupancy(place_bins{rr}{pp}));
        end
    else
        
    end
end

%histogram of rates in id'd fields
figure;
hold on;
title('In-field transient rates for all neurons');
histogram(cell2mat(field_event_rates))

%% Recalculate centroid based on peak with highest transient rate
rr=3;
%tuning vectors for each ROI
tun_vectors{rr} = session_vars{1}.Place_cell{4}.Tuning_Specificity.tuning_vector{rr};
%original sum vector
tun_vector{rr} = session_vars{1}.Place_cell{4}.Tuning_Specificity.tuning_vector_specificity(rr);

%% Turn into function
%input edges and tuning vector
%output - adjusted tuning vector

adjust_tuning_vector()

%convert to tuning vector to bin value for each roi
rad_angles = angle(tun_vectors{rr});
%convert to degrees
deg_angles = rad2deg(rad_angles);
%convert to degree range from 0-360
neg_angles_logical = deg_angles < 0;
%add 360 to negative angles
deg_angles(neg_angles_logical) = 360+deg_angles(neg_angles_logical);
%convert angles to bin position
bins = discretize(deg_angles,1:360/100:360);

%check which vectors fall into the place field bin range
pf_event_keep =  bins >= placeField_edges{rr}(1,1) &   bins <=placeField_edges{rr}(1,2);

%updated vector with adjusted magnitude
pf_unit_vector_kp = sum(tun_vectors{rr}(pf_event_keep))/abs(sum(tun_vectors{rr}(pf_event_keep)));
pf_mag_factor_kp = abs(sum(tun_vectors{rr}(pf_event_keep)))/sum(abs(tun_vectors{rr}(pf_event_keep)));
%update vector
pf_vector = pf_unit_vector_kp*pf_mag_factor_kp;

%convert angle of vectors to bin position
figure;
compass(tun_vectors{rr})
hold on
%tuning vectors with multiple fields
compass(tun_vector{rr},'k');
%tuning vector based on place field with highest transient rate
compass(pf_vector,'r');


%% Plot rate map, place field edges

figure;
for rr=21:22
    hold on;
    %plot unsmoothed/non-normalized event rate
    plot(event_rate(:,rr), 'k');
    %plot identified place fields
    for pp=1:size(placeField_edges{rr},1)
        stem(placeField_edges{rr}(pp,:), [1,1], 'r')
    end
    
    pause
    clf
end



%%

end

