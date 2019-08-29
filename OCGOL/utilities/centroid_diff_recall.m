function [outputArg1,outputArg2] = centroid_diff_recall(session_vars,tunedLogical, pf_vector,field_event_rates,registered,options)

%% Import variables

%use filtered matches
%matching_list = registered.multi.assigned_filtered;
matching_list =registered.multi.matching_list_filtered;

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

%% Extract tuning spec structs and assign reward location
%tuning specificities struct for each session
for ss=1:7
    for tt=1:2
        %original TS
        ts{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity;
        %tuning vectors for each ROI
        %tun_vectors{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector;
        %original sum vector
        tun_vector{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_vector_specificity;
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

%% Maximum transient intensity tuning vectors (select)

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
                else %if more than 1
                    pf_vector_max{ss}{tt}(rr) = pf_vector{ss}{tt}{rr}(max_transient_peak{ss}{tt}(rr));
                end
            else
                pf_vector_max{ss}{tt}(rr) = nan;
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

%for each session
for ss=1:size(session_vars,2)
    %between A and B in first  session
    [theta{ss}] = centroid_angle_diff_btw_ROIs(uCar{ss}{1},uCar{ss}{2},pf_vector_max);
    
    %between A and B in last session
    %[theta{ss}] = centroid_angle_diff_btw_ROIs(uCar{2}{1},uCar{2}{2},pf_vector_max);
    
end
%% Plot the differences as histograms

figure('Position', [2900 500 1300 400]);
subplot(1,3,1)
hold on;
ylim([0 0.6])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('5A5B')
histogram(theta{1}(AandB_tuned{1}),15,'Normalization','probability','BinLimits',[0,pi]);
hold off

subplot(1,3,2)
hold on
ylim([0 0.6])
ylabel('Normalized prob');
xlabel('Angular difference [rad]');
title('Random AB')
histogram(theta{4}(AandB_tuned{4}),15,'Normalization','probability','BinLimits',[0,pi]);
hold off

%plot cdf curves for each
subplot(1,3,3)
hold on
e1 = cdfplot(theta{1}(AandB_tuned{1}));
e2 = cdfplot(theta{2}(AandB_tuned{2}));
gca
xlabel('Angular difference [rad]');
ylabel('Cumulative fraction');
legend([e1,e2], '5A5B','Random AB', 'Location','southeast');
hold off

%% Select centroid diff based on min event/PF filtered match list

ses_comp = [1 6];

%AB SI criterion, at least 1 PF, and min 5 event filtered matching events
%(post-manual filter as well)
filtered_matching_list = matching_list.si_AB_filt_event_filt;

%generate matrix of matching componenets between any 2 days
list_ROI_matches_idx = find(sum(~isnan(filtered_matching_list(:,ses_comp)),2) == 2);

%day 2 day match list
d2d_match_list = filtered_matching_list(list_ROI_matches_idx,ses_comp);

%centroid diff between any 2 days
nanmean(theta{ses_comp(1)}(d2d_match_list(:,1)))
nanmean(theta{ses_comp(2)}(d2d_match_list(:,2)))

figure;
hold on
axis square
xlim([0 3.5])
ylim([0 3.5])
scatter(theta{ses_comp(1)}(d2d_match_list(:,1)), theta{ses_comp(2)}(d2d_match_list(:,2)))

%% Generate match list

% tuning_selection = AandB_tuned;
% %tuned to A or B on either sessions
% select_match_idx{1} = find(tuning_selection{1} ==1);
% select_match_idx{2} = find(tuning_selection{2} ==1);
% 
% %intersect with
% [tuned_match_idx{1},match_idx{1},~] = intersect(matching_list(:,1),select_match_idx{1},'stable');
% [tuned_match_idx{2},match_idx{2},~] = intersect(matching_list(:,2),select_match_idx{2},'stable');
% 
% %create not logical for nan exclusion from copied matrix assignement below
% include_log{1} = false(1,size(matching_list,1));
% include_log{1}(match_idx{1}) = 1; 
% %session 2 
% include_log{2} = false(1,size(matching_list,1));
% include_log{2}(match_idx{2}) = 1;
% 
% %make copy
% tuned_matching_ROI_list = matching_list;
% %nan first session that are not tuned and last session that are not
% %tuned
% tuned_matching_ROI_list(~include_log{1},1) = nan;
% tuned_matching_ROI_list(~include_log{2},2) = nan;
% 
% %which neurons to remove based on tuning criterion
% keep_ROI = sum(isnan(tuned_matching_ROI_list),2) == 0;
% 
% %retain only tuned and matched ROIs
% tuned_matching_ROI_list(~keep_ROI,:) = [];
% 
% %convert to respective logicals
% %ses 1
% select_ses_logi{1} = false(1,size(AandB_tuned{1},2));
% select_ses_logi{1}(tuned_matching_ROI_list(:,1)) = true;
% %ses 2
% select_ses_logi{2} = false(1,size(AandB_tuned{2},2));
% select_ses_logi{2}(tuned_matching_ROI_list(:,2)) = true;

%% Plot the differences as histograms - only matching ROIs
% figure('Position', [2900 500 1300 400]);
% subplot(1,3,1)
% hold on;
% ylim([0 0.3])
% ylabel('Normalized prob');
% xlabel('Angular difference [rad]');
% title('5A5B')
% histogram(theta{1}(select_ses_logi{1}),15,'Normalization','probability','BinLimits',[0,pi]);
% hold off
% subplot(1,3,2)
% hold on
% ylim([0 0.3])
% ylabel('Normalized prob');
% xlabel('Angular difference [rad]');
% title('Random AB')
% histogram(theta{2}(select_ses_logi{2}),15,'Normalization','probability','BinLimits',[0,pi]);
% hold off
% %plot cdf curves for each
% subplot(1,3,3)
% hold on
% e1 = cdfplot(theta{1}(select_ses_logi{1}));
% e2 = cdfplot(theta{2}(select_ses_logi{2}));
% gca
% xlabel('Angular difference [rad]');
% ylabel('Cumulative fraction');
% legend([e1,e2], '5A5B','Random AB', 'Location','southeast');
% hold off

%% ROIs that significantly TS tuned in each session

% for ss=1:size(session_vars,2)
%     theta_ts_sig{ss} = theta{ss}(ts{ss}.significant_ROI);
% end


%% OLD
%% Make matching matrix


%find session 1 matches
% [~,idx_match{1}, ~] = intersect(find(AandB_tuned{1}==1),matching_list(:,1));
% %find session 2 matches
% [~,idx_match{2}, ~] = intersect(find(AandB_tuned{2}==1),matching_list(:,2));
% 
% 
% %convert matching list to logical
% %session 1
% matching_log{1} = false(1,size(AandB_tuned{1},2));
% matching_log{1}(idx_match{1}) = true;
% %session 2
% matching_log{2} = false(1,size(AandB_tuned{2},2));
% matching_log{2}(idx_match{2}) = true;
% 
% find(matching_log{1} ==1)
% find(matching_log{2} ==1)


end

