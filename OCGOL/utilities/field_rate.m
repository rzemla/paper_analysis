function [field_event_rates, field_total_events] = field_rate(event_map,occupancy,placeField_edges)


%%
%for each ROI
for rr=1:size(event_map,2)
    %for each place field
    if ~isempty(placeField_edges{rr})
        for pp=1:size(placeField_edges{rr},1)
            %if place bins split at lap edge
            if placeField_edges{rr}(pp,1) < placeField_edges{rr}(pp,2)
                place_bins{rr}{pp} = placeField_edges{rr}(pp,1):placeField_edges{rr}(pp,2);
            else
                place_bins{rr}{pp} = [placeField_edges{rr}(pp,1):100,1:placeField_edges{rr}(pp,2)];
            end            
            %event rate (number of events across bin range / time spent
            %across bin (s)
            field_event_rates{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr))./sum(occupancy(place_bins{rr}{pp}));
            %total events (not occupancy normalzied)
            field_total_events{rr}(pp) = sum(event_map(place_bins{rr}{pp},rr));
        end
    else
        %set to empty to prevent clipping of vector
        field_event_rates{rr} = [];
        field_total_events{rr} = [];
    end
end

end

