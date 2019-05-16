function [outputArg1,outputArg2] = adjust_tuning_vector(inputArg1,inputArg2)



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

end

