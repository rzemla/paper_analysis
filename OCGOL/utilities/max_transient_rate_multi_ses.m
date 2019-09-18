function [bins,max_transient_peak] = max_transient_rate_multi_ses(session_vars,field_event_rates,pf_vector,options)



%% Create a sort lists based on transient rate
%for each each peak find the idx of the one with highest in-field transient rate
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %for each ROI
        for rr =1:size(field_event_rates{ss}{tt},2)
            %if more than 1 field
            if size(field_event_rates{ss}{tt}{rr},2) > 1
                [~,max_transient_peak{ss}{tt}(rr)] = max(field_event_rates{ss}{tt}{rr});

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
                    sort_vector{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                else %if more than 1
                    sort_vector{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                end
            else
                sort_vector{ss}{tt}(rr) = nan;
            end
        end
    end
end

%get the bin location of sort tuning vector
%get degree angle of each and classify in 1-100 bin
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        sort_vec_angle{ss}{tt} = angle(sort_vector{ss}{tt});
    end
end

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



end

