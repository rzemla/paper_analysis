function [outputArg1,outputArg2] = TS_score_diff(session_vars,tunedLogical, registered)

%% Import variables
%list of matching ROIs between sessions
matching_list = registered.multi.assigned_all;

%use ts tuning for this analysis
options.tuning_criterion = 'ts';

%tuning specificty scores in A and B trials from each session
for ss=1:2
    for tt=4:5
        ts_scores{ss}{tt} = session_vars{ss}.Place_cell{tt}.Tuning_Specificity.tuning_specificity;
        
    end
end


%% Get tuning logicals
%make conditional here for si or ts tuned neurons
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


%% Taking matching ROI scores for all ROIs (regardless of stat sig. tuning)

%for A
ts_scores_match = [ts_scores{1}{4}(matching_list(:,1))', ts_scores{2}{4}(matching_list(:,2))'];

%tuned on day 2 to A
Atuned_all_idx{2} = find(Atuned{1, 2} ==1);

%find correspoinding indices in matching ROI list and get scores
[intersects,match_idx_Ad2,~] = intersect(matching_list(:,2),Atuned_all_idx{2},'stable');
%select corresponding TS scores
ts_scores_match_tuned_A_ses2 = [ts_scores{1}{4}(matching_list(match_idxAd2,1))', ts_scores{2}{4}(matching_list(match_idxAd2,2))'];



%plot scatter
figure;
hold on
histogram(ts_scores_match_tuned_A_ses2(:,1),'DisplayStyle','stairs','Normalization','probability')
hold on
histogram(ts_scores_match_tuned_A_ses2(:,2),'DisplayStyle','stairs','Normalization','probability')

figure; 
hold on
ecdf(ts_scores_match_tuned_A_ses2(:,1))
ecdf(ts_scores_match_tuned_A_ses2(:,2))
gca
xlabel('Tuning specificity');
ylabel('Cumulative fraction');

end

