function [pf_vector] = adjust_tuning_vector(tun_vectors,tun_vector,placeField_edges,options)
%define the tuning vectors for each place field and outputs tuning vector
%for each field in the pf_vector cell

%QC checked

%% Set variables



%% Define which bin each tuning vector (event) lies in
%for each ROI
for rr=1:size(placeField_edges,2)
    %convert to tuning vector to bin value for each roi
    rad_angles{rr} = angle(tun_vectors{rr});
    %convert to degrees
    deg_angles{rr} = rad2deg(rad_angles{rr});
    %convert to degree range from 0-360
    neg_angles_logical{rr} = deg_angles{rr} < 0;
    %add 360 to negative angles
    deg_angles{rr}(neg_angles_logical{rr}) = 360 + deg_angles{rr}(neg_angles_logical{rr});
    %convert angles to bin position
    bins{rr} = discretize(deg_angles{rr},0:360/100:360);
end

%% Define updated tuning vector for each place field
%bin of each tuning vector

%for each ROI
for rr=1:size(placeField_edges,2)
    %for each place field
    %if at least 1 field detected
    if size(placeField_edges{rr},1) ~= 0
        
        for pp=1:size(placeField_edges{rr},1)
            if ~isempty(placeField_edges{rr}) %if not empty
                %check which vectors fall into the place field bin range
                %if place bins split at lap edge
                if placeField_edges{rr}(pp,1) < placeField_edges{rr}(pp,2)
                    pf_event_keep{rr}{pp} =  bins{rr} >= placeField_edges{rr}(pp,1) & bins{rr} <=placeField_edges{rr}(pp,2);
                else
                    %disp(rr)
                    pf_event_keep{rr}{pp} =  (bins{rr} >= placeField_edges{rr}(pp,1) & bins{rr} <= 100) | ...
                        (bins{rr} >= 1 & bins{rr} <= placeField_edges{rr}(pp,2));
                    
                end
                
                %updated vector with adjusted magnitude
                pf_unit_vector_kp{rr}(pp) = sum(tun_vectors{rr}(pf_event_keep{rr}{pp}))/abs(sum(tun_vectors{rr}(pf_event_keep{rr}{pp})));
                pf_mag_factor_kp{rr}(pp) = abs(sum(tun_vectors{rr}(pf_event_keep{rr}{pp})))/sum(abs(tun_vectors{rr}(pf_event_keep{rr}{pp})));
                %update vector
                pf_vector{rr}(pp) = pf_unit_vector_kp{rr}(pp)*pf_mag_factor_kp{rr}(pp);
                
            end
        end
    else
        %set empty to prevent cutting vector short on output
        pf_vector{rr} = [];
    end
    
end

%% Plot example - checked and works as expected

if options.pf.skipDisplay == 1
    figure;
    for rr=1:size(pf_vector,2)
        %tuning vectors with multiple fields
        compass(tun_vectors{rr})
        hold on
        title(num2str(rr));
        %original tuning vector
        compass(tun_vector(rr),'k');
        %tuning vector based on place field with highest transient rate
        compass(pf_vector{rr},'r');
        pause()
        clf
    end
end

end

