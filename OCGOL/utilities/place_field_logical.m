function [pf_count_filtered_log] = place_field_logical(select_fields,options)


%% Count number of detected PFs based on filter for 5 lap-distinct events in each field (for each ROI) - global (no selection)

%for each trial type
for tt =1:size(options.selectTrial,2)
    for rr =1:size(select_fields{1}{tt},2)
        %count the place fields for each roi
        %row is trial type
        %column is each ROI
        pf_count_filtered(tt,rr)= sum(select_fields{1}{tt}{rr});
    end
end

%create logical where there is at least 1 field with at least 5
%lap-distinct events
pf_count_filtered_log = ~(pf_count_filtered == 0);


end

