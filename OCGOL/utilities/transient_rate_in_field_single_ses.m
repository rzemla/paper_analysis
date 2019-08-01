function [field_event_rates,pf_vector] = transient_rate_in_field_single_ses(session_vars,options)
%input - cell of structs with animal data

%% Define input variables 

%try for 1 session first

%for each session
for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        %Place field edge data
        placeField_edges{ss}{tt} = session_vars{ss}.Place_cell{tt}.placeField.edge;
        
        %Event rate for 100 bins
        event_rate{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.rate_map{8};
        
        %event map (not rate)% 100 bins
        event_map{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.event_map{8};
        
        %Occupancy for 100 bins (seconds, not normalized (0-1))
        occupancy{ss}{tt} = session_vars{ss}.Place_cell{tt}.Spatial_Info.occupancy_map{8};
    end
end

%should equal event rate (it does)
%rate_map = event_map./occupancy';
%isequal(rate_map,event_rate)

%% Find and store transient rate in each place field

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct (4,5) or only correct
    %(1,2)
    for tt=options.selectSes
        %find field event rate for each session and trial types
        field_event_rates{ss}{tt} = field_rate(event_map{ss}{tt},occupancy{ss}{tt},placeField_edges{ss}{tt});
    end
end

%which session
session_nb = 1;

%histogram of rates in id'd fields
figure;
subplot(2,1,1)
hold on;
xlim([0 1])
title('In-field transient rates for all neurons - A trials');
histogram(cell2mat(field_event_rates{session_nb}{options.selectSes(1)}))
subplot(2,1,2)
hold on
xlim([0 1])
title('In-field transient rates for all neurons - B trials');
histogram(cell2mat(field_event_rates{session_nb}{options.selectSes(2)}))

%% Recalculate centroid based on peak with highest transient rate

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        %tuning vectors for each ROI
        tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
    end
end

%% Turn into function
%input edges and tuning vector
%output - adjusted tuning vector

options.pf.skipDisplay = 0;

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=options.selectSes
        [pf_vector{ss}{tt}] = adjust_tuning_vector(tun_vectors{ss}{tt},tun_vector{ss}{tt},placeField_edges{ss}{tt},options);
    end
end

%% Plot rate map, place field edges - use for debug, final check

% figure;
% for rr=21:22
%     hold on;
%     %plot unsmoothed/non-normalized event rate
%     plot(event_rate(:,rr), 'k');
%     %plot identified place fields
%     for pp=1:size(placeField_edges{rr},1)
%         stem(placeField_edges{rr}(pp,:), [1,1], 'r')
%     end
%     
%     pause
%     clf
% end



%%

end

