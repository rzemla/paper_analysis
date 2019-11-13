function [field_event_rates, field_total_events] = field_rate(event_map,occupancy,placeField_edges)

%QC checked

%% Get the bins associated with each place field for each ROI and divide by occupancy(s) in those bins

%normalized occupancy
norm_occupancy = occupancy./sum(occupancy);

%for each ROI
for rr=1:size(event_map,2)
    %if there is a place field
    if ~isempty(placeField_edges{rr})
        %for each place field
        for pp=1:size(placeField_edges{rr},1)
            %get the bin range for the rate calculation
            %if no split
            if placeField_edges{rr}(pp,1) < placeField_edges{rr}(pp,2)
                place_bins{rr}{pp} = placeField_edges{rr}(pp,1):placeField_edges{rr}(pp,2);
            else %if place bins split at lap edge
                place_bins{rr}{pp} = [placeField_edges{rr}(pp,1):100,1:placeField_edges{rr}(pp,2)];
            end            
            %event rate (number of events across bin range / time spent
            %across bin (s)
            field_event_rates{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr))./sum(occupancy(place_bins{rr}{pp}));
            
            %event rate relative to normalized occupancy (should yield same
            %outcome in terms relative order and max)
            field_event_rates_norm{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr))./sum(norm_occupancy(place_bins{rr}{pp}));
            
            %total events (not occupancy normalzied)
            field_total_events{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr));
        end
    else
        %set to empty to prevent clipping of vector
        field_event_rates{rr} = [];
        field_event_rates_norm{rr} = [];
        field_total_events{rr} = [];
    end
end


end

