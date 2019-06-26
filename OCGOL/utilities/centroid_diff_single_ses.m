function [cent_diff_AandB, pf_vector_max] = centroid_diff_single_ses(session_vars,tunedLogical, pf_vector,field_event_rates,options)

%input


%output:
%cent_diff_A
% pf_vector_max = tuning vector associated with max in-field transient rate


%% Define tuned combinations

switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            
            all_neurons{ss} = true(size(Atuned{ss}));
        end
   case 'ts' %spatial information 
        for ss =1:size(session_vars,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            
            all_neurons{ss} = true(size(Atuned{ss}));
        end
end


%% Extract tuning spec structs and assign reward location, place field bin centers
%tuning specificities struct for each session
for ss=1:size(session_vars,2)
    for tt=1:2
        %original TS
        ts{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity;
        %tuning vectors for each ROI
        %tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
        
        place_field_centers{ss}{tt} = session_vars{1}.Place_cell{tt}.placeField.center;
        %cm
        %rewardLoc(ii) = session_vars{ii}.Behavior.performance.reward_loc;
    end
end
    
%% Get index of max intensity for place field vector

for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=1:2
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

%% Maximum transient intensity tuning vectors (select) and bin centers (100 bin space)

%select the place field vectors based on the max intensity analysis
for ss=1:size(session_vars,2)
    %for each trial (A or B) regardless if correct
    for tt=1:2
        %for each ROI
        for rr=1:size(pf_vector{ss}{tt},2)
            %if not empty
            if ~isempty(pf_vector{ss}{tt}{rr})
                %if single field, use original vector
                if size(pf_vector{ss}{tt}{rr},2) == 1
                    pf_vector_max{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                    %place field bin center
                    place_field_centers_max{ss}{tt}(rr) = place_field_centers{ss}{tt}{rr};
                else %if more than 1
                    pf_vector_max{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                    place_field_centers_max{ss}{tt}(rr) = place_field_centers{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                end
            else
                pf_vector_max{ss}{tt}(rr) = nan;
                place_field_centers_max{ss}{tt}(rr) = nan;
            end
        end
    end
end

%% 

%pf_vector_max{1, 1}{1, 4}

%% Get each tuning vector (complex and cartesian) for each max place field vector

%for each session
for ss=1:size(session_vars,2)
    %for each trial
    for tt=1:2
        %for each ROI in that session
        for ROI = 1:size(pf_vector_max{ss}{tt},2)
            %get angle for that ROI from each session relative to 0
            %angle of each vector
            %ang{ss}{tt} = angle(pf_vector_max{ss}{tt});
            ang{ss}{tt}(ROI) = angle(pf_vector_max{ss}{tt}(ROI));
            %ang{ss}(ROI) = angle(ts{ss}.tuning_vector_specificity(ROI));
            
            %complex vector
            %u{ss}{tt} = pf_vector_max{ss}{tt};
            u{ss}{tt}(ROI) = pf_vector_max{ss}{tt}(ROI);
            %u{ss}(ROI) = ts{ss}.tuning_vector_specificity(ROI);
            %cartesian converted
            %uCar{ss}{tt} = [real(u{ss}{tt})', imag(u{ss}{tt})'];
            uCar{ss}{tt}{ROI} = [real(u{ss}{tt}(ROI)), imag(u{ss}{tt}(ROI))];
            %uCar{ss}{ROI} = [real(u{ss}(ROI)), imag(u{ss}(ROI))];
        end
    end
end

%% Calculate difference from each centroid to 

%for each session
%largest angle difference should be pi
% for ss=1:size(session_vars,2)
%     [theta{ss}] = centroid_angle_diff_btw_ROIs(uCar{ss},rewardLoc(ss));
% end
    %between A and B in first  session 
[theta{1}] = centroid_angle_diff_btw_ROIs(uCar{1}{1},uCar{1}{2},pf_vector_max);
 
    %between A and B in last session 
%[theta{2}] = centroid_angle_diff_btw_ROIs(uCar{2}{4},uCar{2}{5},pf_vector_max);

%% Plot the differences as histograms
figure('Position', [2900 500 1300 400]);
subplot(1,2,1)
hold on;
ylim([0 0.3])
ylabel('Normalized density');
xlabel('Angular difference [rad]');
xticks([pi/4 pi/2 3*pi/4,pi])
xticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
title('Centroid difference')
histogram(theta{1}(AandB_tuned{1}),15,'Normalization','probability','BinLimits',[0,pi]);
hold off

%plot cdf curves for each
subplot(1,2,2)
hold on

e1 = cdfplot(theta{1}(AandB_tuned{1}));
grid off

gca
xlabel('Angular difference [rad]');
xticks([pi/4 pi/2 3*pi/4,pi])
xticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
ylabel('Cumulative fraction');
hold off

%% Export centroid diffs for mutually tuned neurons

%get only mutually tuned ROIs (A), convert to matrix and combined into 1
place_field_centers_max_AandB{1}{1} = place_field_centers_max{1}{1}(AandB_tuned{1}); 
place_field_centers_max_AandB{1}{2} = place_field_centers_max{1}{2}(AandB_tuned{1});
%combined into single matrix
comb_place_field_max_AB = [place_field_centers_max_AandB{1}{1};place_field_centers_max_AandB{1}{2}];

%export as struct
cent_diff_AandB.angle_diff = theta{1}(AandB_tuned{1});
cent_diff_AandB.max_bin = comb_place_field_max_AB;

end

