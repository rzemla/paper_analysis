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
event_map = session_vars{1}.Place_cell{1, 4}.Spatial_Info.event_map{8}; 

%Occupancy for 100 bins
occupancy = session_vars{1}.Place_cell{4}.Spatial_Info.occupancy_map{8};

%all A and all B (regardless if correct)
%session_vars{1}.Place_cell.
%end

%should equal event rate (it does_
rate_map = event_map./occupancy';
%isequal(rate_map,event_rate)

%% Assign transient rate for each place field


%% Plot rate map, place field edges

figure;
for rr=1:20
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

