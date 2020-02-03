function [mean_PV,sem_PV] = filter_convert_day_return_mean_sem_PV(exp_struct,excl_day_combined_day_nan,exp_type,day_range)

%% Extract PV correlation matrices

%define number of animals
nb_animal = sum(~cellfun(@isempty,excl_day_combined_day_nan(:,exp_type)));

%relative to day 1
%for each animal
for aa=1:nb_animal
    %for each day
    for dd=day_range
        %look for day and get session index that corresponds to that day
        ses2day_idx = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == dd);
        
        %if no assignmenet for that day
        if ~isempty(ses2day_idx)
            %animal x day PV matrix cell
            %for A trials
            day_PV_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A{1,ses2day_idx};
            %for B trials
            day_PV_mat.exp.B{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B{1,ses2day_idx};
        else
            day_PV_mat.exp.A{aa,dd} = [];
            day_PV_mat.exp.B{aa,dd} = [];
        end
        
    end
    
end

%% Make substitutions for Day 4 and Day 5 for short term recall animals
%only do for short term recall animal, i.e. exp_type #2
if exp_type == 2
    for aa=1:nb_animal
        
        %for day 4 substitution = Day 2 vs. Day 6
        %for day 5 substitution = Day 2 vs. Day 7
        for dd=4:5
            %if making day 4 substitution, get day 2 and day 6 ses
            if dd==4
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 2);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 6);
            elseif dd==5
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 2);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 7);
            end
            %if no assignmenet for that day (make sure both days exist in this case)
            if size(ses2day_idx,2) == 2
                %animal x day PV matrix cell
                %for A trials
                day_PV_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A{ses2day_idx(1),ses2day_idx(2)};
                %for B trials
                day_PV_mat.exp.B{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B{ses2day_idx(1),ses2day_idx(2)};
            else
                day_PV_mat.exp.A{aa,dd} = [];
                day_PV_mat.exp.B{aa,dd} = [];
                
            end
        end
    end
    
end


%% Get diagnonal and mean of PV for each day

%get diagnonal for A/B
day_PV_diag.exp.A = cellfun(@diag,day_PV_mat.exp.A,'UniformOutput',false);
day_PV_diag.exp.B = cellfun(@diag,day_PV_mat.exp.B,'UniformOutput',false);

%get mean of diagnonal for A
day_PV_diag_mean.exp.A = cell2mat(cellfun(@nanmean,day_PV_diag.exp.A,'UniformOutput',false));
day_PV_diag_mean.exp.B = cell2mat(cellfun(@nanmean,day_PV_diag.exp.B,'UniformOutput',false));

%get mean across animals of the mean of the diagnonals
day_PV_diag_mean_mean_exp.A = nanmean(day_PV_diag_mean.exp.A,1);
day_PV_diag_mean_mean_exp.B = nanmean(day_PV_diag_mean.exp.B,1);

%get sem across animals using correct count of animals per day
day_PV_diag_sem_mean_exp.A = nanstd(day_PV_diag_mean.exp.A,0,1)./sqrt(sum(~isnan(day_PV_diag_mean.exp.A),1));
day_PV_diag_sem_mean_exp.B = nanstd(day_PV_diag_mean.exp.B,0,1)./sqrt(sum(~isnan(day_PV_diag_mean.exp.B),1));

%% For export
mean_PV.A = day_PV_diag_mean_mean_exp.A;
mean_PV.B = day_PV_diag_mean_mean_exp.B;

sem_PV.A = day_PV_diag_sem_mean_exp.A;
sem_PV.B = day_PV_diag_sem_mean_exp.B;


end

