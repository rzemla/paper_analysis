function [cent_diff,pf_vector_max] = centroid_diff_multi_ses(session_vars,tunedLogical, pf_vector,field_event_rates,select_fields,registered,options)

%% Import variables

matching_list = registered.multi.assigned_all;


%% Define tuned combinations

switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =options.sessionSelect
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
        for ss =options.sessionSelect
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

%% Extract tuning spec structs and assign reward location
%tuning specificities struct for each session
for ss=options.sessionSelect
    for tt=options.selectTrial
        %original TS
        ts{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity;
        %tuning vectors for each ROI
        %tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;

        place_field_centers{ss}{tt} = session_vars{ss}.Place_cell{tt}.placeField.center;
        %cm
        %rewardLoc(ii) = session_vars{ii}.Behavior.performance.reward_loc;
    end
end
    
%% Get index of max intensity for place field vector

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

%%  Maximum transient intensity tuning vectors (select) and bin centers (100 bin space)

%select the place field vectors based on the max intensity analysis
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %for each ROI
        for rr=1:size(pf_vector{ss}{tt},2)
            %if not empty
            if ~isempty(pf_vector{ss}{tt}{rr})
                %if single field, use original vector
                if size(pf_vector{ss}{tt}{rr},2) == 1
                    pf_vector_max{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                    place_field_centers_max{ss}{tt}(rr) = place_field_centers{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));

                    %pf_vector_max{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                    %place field bin center
                    %place_field_centers_max{ss}{tt}(rr) = place_field_centers{ss}{tt}{rr};
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

%% Bin location of corresponding centroid (signifcant PFs using select_fields as input)
%select the place field vectors based on the max intensity analysis
for ss=options.sessionSelect
    %for each trial (A or B) regardless if correct
    for tt=options.selectTrial
        %for each ROI
        for rr=1:size(pf_vector{ss}{tt},2)
            %if not empty
            if ~isempty(pf_vector{ss}{tt}{rr})
                %if single field, use original vector
                if size(pf_vector{ss}{tt}{rr},2) == 1
                    %pf_vector_max{ss}{tt}(rr) = tun_vector{ss}{tt}(rr);
                    %place field bin center
                    place_field_centers_select_all{ss}{tt}{rr} = place_field_centers{ss}{tt}{rr}(select_fields{ss}{tt}{rr}');
                else %if more than 1
                    %pf_vector_max{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                    place_field_centers_select_all{ss}{tt}{rr} = place_field_centers{ss}{tt}{rr}(select_fields{ss}{tt}{rr}');
                end
            else
                %pf_vector_max{ss}{tt}(rr) = nan;
                place_field_centers_select_all{ss}{tt}{rr} = nan;
            end
        end
    end
end


%% Get each tuning vector (complex and cartesian) for each max place field vector

%for each session
for ss=options.sessionSelect
    %for each trial
    for tt=options.selectTrial
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
%for each session
for ss=options.sessionSelect
    %between A and B in first  session
    [theta{ss}] = centroid_angle_diff_btw_ROIs(uCar{ss}{options.selectTrial(1)},uCar{ss}{options.selectTrial(2)},pf_vector_max{ss});
end

%% Plot the differences as histograms
figure('Position', [2900 500 1300 400]);
subplot(1,3,1)
hold on;
ylim([0 0.3])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('5A5B')
histogram(theta{1}(AandB_tuned{1}),15,'Normalization','probability','BinLimits',[0,pi]);
hold off
subplot(1,3,2)
hold on
ylim([0 0.3])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('Random AB')
histogram(theta{6}(AandB_tuned{6}),15,'Normalization','probability','BinLimits',[0,pi]);
hold off
%plot cdf curves for each
subplot(1,3,3)
hold on
e1 = cdfplot(theta{1}(AandB_tuned{1}));
e2 = cdfplot(theta{6}(AandB_tuned{6}));
gca
xlabel('Angular difference [rad]');
ylabel('Cumulative fraction');
legend([e1,e2], '5A5B','Random AB', 'Location','southeast');
hold off


%% Extract matching neurons tuning between any 2 session and fitlered by event criteria

ses_comp = [1,3];

select_match = sum(~isnan(registered.multi.matching_list_filtered.si_AB_filt_event_filt(:,ses_comp)),2)==2;

ses_match_theta{ses_comp(1)} = theta{ses_comp(1)}(registered.multi.matching_list_filtered.si_AB_filt_event_filt(select_match,ses_comp(1)));

ses_match_theta{ses_comp(2)} =theta{ses_comp(2)}(registered.multi.matching_list_filtered.si_AB_filt_event_filt(select_match,ses_comp(2)));


%% Plot the differences as histograms - only matching ROIs
figure('Position', [2900 500 1300 400]);
subplot(1,3,1)
hold on;
ylim([0 0.3])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('5A5B')
histogram(ses_match_theta{ses_comp(1)},15,'Normalization','probability','BinLimits',[0,pi]);
hold off
subplot(1,3,2)
hold on
ylim([0 0.3])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('Random AB')
histogram(ses_match_theta{ses_comp(2)},15,'Normalization','probability','BinLimits',[0,pi]);
hold off
%plot cdf curves for each
subplot(1,3,3)
hold on
e1 = cdfplot(ses_match_theta{ses_comp(1)});
e2 = cdfplot(ses_match_theta{ses_comp(2)});
gca
xlabel('Angular difference [rad]');
ylabel('Cumulative fraction');
legend([e1,e2], '5A5B','Random AB', 'Location','southeast');
hold off

[h,p] = ttest2(ses_match_theta{ses_comp(1)},ses_match_theta{ses_comp(2)})

%% Export centroid diffs for all neurons

%combined into single matrix
for ss=options.sessionSelect
    comb_place_field_max{ss} = [place_field_centers_max{ss}{options.selectTrial(1)};place_field_centers_max{ss}{options.selectTrial(1)}];
end

%ALL NEURONS
cent_diff.angle_diff = theta;
cent_diff.max_bin = comb_place_field_max;
%place field bin centers for all fields with significant fields (min 5
%event on distinct laps and PFs id'd using PF_finder
cent_diff.all_sig_bin = place_field_centers_select_all;

end

