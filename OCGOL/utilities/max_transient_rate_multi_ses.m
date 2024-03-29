function [bins,max_transient_peak] = max_transient_rate_multi_ses(session_vars,field_event_rates,pf_vector,options)
%QC checked


%% Variable descriptions
%field_event_rates = in field rates (events/s) for each id'd place field
%pf_vector = tuning vector associated with each id'd place field


%% Create a sort lists based on transient rate
%for each each peak find the idx of the one with highest in-field transient rate
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %for each ROI
        for rr =1:size(field_event_rates{ss}{tt},2)
            %if more than 1 field
            if size(field_event_rates{ss}{tt}{rr},2) > 1
                %get index of field with max rate
                [~,max_transient_peak{ss}{tt}(rr)] = max(field_event_rates{ss}{tt}{rr});
            %if single field, assign place field index of 1
            elseif size(field_event_rates{ss}{tt}{rr},2) == 1 %if 1
                max_transient_peak{ss}{tt}(rr) = 1;
            else
                max_transient_peak{ss}{tt}(rr) = nan;
            end
        end
    end
end

%% Get centroid for each neuron (bin of centroid)

%two approaches (use 1 for now) - can set up as conditional option
%1) use global vector for single field and place field vector
%for neurons with more than 1 vector
%2) use max field vectors regardless if 1 or multiple fields

%get tuning vectors from tuning specificity calculation
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %tuning vectors for each ROI
        tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
    end
end

%convert centroid to bin location
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %for each ROI
        for rr=1:size(pf_vector{ss}{tt},2)
            %if not empty
            if ~isempty(pf_vector{ss}{tt}{rr})
                %if single field, use original vector
                if size(pf_vector{ss}{tt}{rr},2) == 1
                    %select whether to use TS vector of adjusted vector for
                    %cells with single fields
                    if options.select_adj_vec == 1
                        sort_vector{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                    else %select TS calculated vector
                        sort_vector{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                    end
                    
                else %if more than 1
                    sort_vector{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                end
            else
                sort_vector{ss}{tt}(rr) = nan;
            end
        end
    end
end

%% Get the angle of each vector

%get degree angle of each and classify in 1-100 bin
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        sort_vec_angle{ss}{tt} = angle(sort_vector{ss}{tt});
    end
end


%% Convert vector angle to track bins (1-100)

for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %convert rad 2 deg and assign spatial bins
        deg_angles{ss}{tt} = rad2deg(sort_vec_angle{ss}{tt});
        %convert to degree range from 0-360
        neg_angles_logical{ss}{tt} = deg_angles{ss}{tt} < 0;
        %add 360 to negative angles
        deg_angles{ss}{tt}(neg_angles_logical{ss}{tt}) = 360 + deg_angles{ss}{tt}(neg_angles_logical{ss}{tt});
        %convert angles to bin position
        bins{ss}{tt} = discretize(deg_angles{ss}{tt},0:360/100:360);
    end
end

%% QC plot
if 0
    figure
    tt=1
    for rr=1:298
        subplot(1,2,1)
        hold on
        ylabel('Event count')
        plot(1:100,session_vars{1}.Place_cell{tt}.Spatial_Info.event_map{8}(:,rr))
        %plot max field center
        plot([bins{1}{tt}(rr) bins{1}{tt}(rr)], [0 2],'r')
        
        subplot(1,2,2,polaraxes)
        hold on
        polarplot([0 angle(sort_vector{1}{tt}(rr))],[0 abs(sort_vector{1}{tt}(rr))])
        pause
        clf
    end
end



end

