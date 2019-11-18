function [pf_data] = extract_pf_width_count(idxs, trial_type, session_vars, select_fields, pf_count_filtered,binToCm)


%count the field bin width of task-selective neurons (includes both sig and
%in-sig field - must be filtered below)
bin_width_not_filt = session_vars{1}.Place_cell{trial_type}.placeField.width(idxs);

%task selective neurons (only fields that are significant) - filtered
field_count = pf_count_filtered(trial_type,idxs);

%select fields that match the place field crtieria (logicals corresponding
%to sign field 
field_select = select_fields{1}{trial_type}(idxs);

%count number of fields
for ii=1:3
    if ii < 3
        field_count_total(ii) = size(find(field_count == ii),2);  
        
    elseif ii ==3 %3 or more fields
        field_count_total(ii) = size(find(field_count >= ii),2);

    end
end

%% Get filtered bin width for each ROI
%for each ROI in class
for rr=1:size(bin_width_not_filt,2)
    %extract the relevant fields
    field_bin_width_filt{rr} = bin_width_not_filt{rr}(field_select{rr});
end

%% Convert bin width to cm from selective ROIs
%combine widths into single vector and multiply by conversion factor
width_cm = cell2mat(field_bin_width_filt)*binToCm;

%% Export as struct
pf_data.bin_width_not_filt = bin_width_not_filt;
pf_data.field_count = field_count;
pf_data.field_select = field_select;
pf_data.field_count_total = field_count_total;
pf_data.field_bin_width_filt = field_bin_width_filt;
pf_data.width_cm = width_cm;


end

