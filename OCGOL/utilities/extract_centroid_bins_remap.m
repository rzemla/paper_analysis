function [bin_center] = extract_centroid_bins_remap(cent_diff,task_remapping_ROIs, select_fields,partial_field_idx)

%% Assign parameters
sig_center_bins = cent_diff.all_sig_bin;


%% For each class of neurons extract bin centers from each trial type

%for each type
for tt=1:2

    %for each ROI for common ROIs
    for rr=1:size(task_remapping_ROIs.common,2)
        bin_center.common(tt,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.common(rr)};
    end

    %for each ROI for globar near ROIs
    for rr=1:size(task_remapping_ROIs.global_near,2)
        bin_center.global_near(tt,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.global_near(rr)};
    end

    %for each ROI for globar far ROIs
    for rr=1:size(task_remapping_ROIs.global_far,2)
        bin_center.global_far(tt,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.global_far(rr)};
    end

    %for each ROI for rate remapping ROIs
    for rr=1:size(task_remapping_ROIs.rate,2)
        bin_center.rate(tt,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.rate(rr)};
    end


    %for each ROI for rate remapping ROIs
    for rr=1:size(task_remapping_ROIs.partial,2)
        if tt==1
            bin_center.partial(1:2,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.partial(rr)};
        elseif tt==2
            bin_center.partial(3:4,rr) = sig_center_bins{1}{tt}{task_remapping_ROIs.partial(rr)};
        end
        
    end
end

%% Remove duplicate finds from partial remappers and create common and partial locations 

bin_center.partial(2,(bin_center.partial(1,:)-bin_center.partial(2,:)) == 0) = nan;
bin_center.partial(4,(bin_center.partial(3,:)-bin_center.partial(4,:)) == 0) = nan;

%isolate the common and partial fields centers from the remove duplicate remappers in matrix
%above and using far field idx from remapper seek script
%for each type
for tt=1:2
    
    for rr=1:size(task_remapping_ROIs.partial,2)
        if tt==1
            if sum(isnan(bin_center.partial(1:2,rr))) ==1
                bin_center.partial_com(tt,rr) = bin_center.partial(1,rr);
                %set partial to nan
                bin_center.partial_far(tt,rr) = nan;
                
            else
                if partial_field_idx(rr) ==1
                    bin_center.partial_com(tt,rr) = bin_center.partial(2,rr);
                elseif partial_field_idx(rr) ==2
                    bin_center.partial_com(tt,rr) = bin_center.partial(1,rr);
                end
                %assign partial field
                bin_center.partial_far(tt,rr) = bin_center.partial(partial_field_idx(rr),rr);
            end
            
        elseif tt==2
            if sum(isnan(bin_center.partial(3:4,rr))) ==1
                bin_center.partial_com(tt,rr) = bin_center.partial(3,rr);
                %set partial to nan
                bin_center.partial_far(tt,rr) = nan;
            else
                if partial_field_idx(rr) ==1
                    bin_center.partial_com(tt,rr) = bin_center.partial(4,rr);
                elseif partial_field_idx(rr) ==2
                    bin_center.partial_com(tt,rr) = bin_center.partial(3,rr);
                end
                %assign partial field
                bin_center.partial_far(tt,rr) = bin_center.partial(2+partial_field_idx(rr),rr);
            end
            
        end
        
    end
end


%%
end

