function [mean_PV,sem_PV,raw] = filter_convert_day_return_mean_sem_PV(exp_struct,excl_day_combined_day_nan,exp_type,day_range)

%% Extract PV correlation matrices relative to D1

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

%% Make substitutions for Day 4 and Day 5 for short term recall animals as a substitute for time distance relative to D1
%only do for short term recall animal, i.e. exp_type #2
if exp_type == 2
    for aa=1:nb_animal
        
        %for day 4 substitution = Day 2 vs. Day 6
        %for day 5 substitution = Day 2 vs. Day 7
        for dd=4:5
            %if making day 4 substitution, get day 2 and day 6 ses
            if dd==4
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 3);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 6);
            elseif dd==5
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 3);
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

%% Extract all day 2 day cross correlations for short-term recall animals for Figure 4F plots
%run only for short-term recall since assuming behavioral steady-state
%there

%if short term recall experiment
if exp_type == 2
    %extract PV correlation matrices for each animal that are day matched
    for aa=1:nb_animal
        %for each day
        for dd1 = day_range
            %for each day
            for dd2 = day_range
                %look for day and get session index that corresponds to that day
                ses2day_idx_d1 = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == dd1);
                ses2day_idx_d2 = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == dd2);
                
                %if both day indices are not empty
                if ~isempty(ses2day_idx_d1) && ~isempty(ses2day_idx_d2)
                    %animal x day PV matrix cell
                    %for A trials
                    crossday_PV_mat.exp(aa).A{dd1,dd2} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A{ses2day_idx_d1,ses2day_idx_d2};
                    %for B trials
                    crossday_PV_mat.exp(aa).B{dd1,dd2} = exp_struct.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B{ses2day_idx_d1,ses2day_idx_d2};
                else
                    crossday_PV_mat.exp(aa).A{dd1,dd2} = [];
                    crossday_PV_mat.exp(aa).B{dd1,dd2} = [];
                end
                
            end
        end
    end
    
    %set empty cells in cross-day for each animal to nan
    %return cell indices where the cells are empty
    %for each animal
    for aa=1:nb_animal
        
        %for A trials
        empty_log = cellfun(@isempty,crossday_PV_mat.exp(aa).A);
        %set the empty cells to nan
        crossday_PV_mat.exp(aa).A(empty_log) = {nan};
        
        %for B trials
        empty_log = cellfun(@isempty,crossday_PV_mat.exp(aa).B);
        %set the empty cells to nan
        crossday_PV_mat.exp(aa).B(empty_log) = {nan};
    end
    
    
    %extract each diagnonal and upward for each day-by-day match
    %for each animal
    for aa=1:nb_animal
        %for each off diagonal
        %for A
        crossday_PV_main_diag(aa).A = cellfun(@(x) diag(x,0),crossday_PV_mat.exp(aa).A,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_PV_main_diag_mean(aa).A = cell2mat(cellfun(@nanmean,crossday_PV_main_diag(aa).A,'UniformOutput',false));
        %for B
        crossday_PV_main_diag(aa).B = cellfun(@(x) diag(x,0),crossday_PV_mat.exp(aa).B,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_PV_main_diag_mean(aa).B = cell2mat(cellfun(@nanmean,crossday_PV_main_diag(aa).B,'UniformOutput',false));
        
    end
    
    for aa=1:nb_animal
        for d_idx=0:(size(crossday_PV_main_diag_mean(aa).A,1) -1)
            %take each diagnonal starting from main
            %for A trials
            PV_diag_distance(aa).A{d_idx+1} = diag(crossday_PV_main_diag_mean(aa).A,d_idx);
            %for B trials
            PV_diag_distance(aa).B{d_idx+1} = diag(crossday_PV_main_diag_mean(aa).B,d_idx);
        end
        %take mean of each diagonal day distance for A (animal x day
        %separation)
        %A trials
        PV_diag_distance_mean.A(aa,:) = cellfun(@nanmean, PV_diag_distance(aa).A);
        %B trials
        PV_diag_distance_mean.B(aa,:) = cellfun(@nanmean, PV_diag_distance(aa).B);
    end
    
    
end


%% Get diagnonal and mean of PV for each day / mean / sem from individually matched sessions

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

%% Cross-correlation distance data
%run only for short term recall data

if exp_type == 2
    %get mean from each day-distance correlation
    day_PV_distance_mean_mean.A = nanmean(PV_diag_distance_mean.A,1);
    day_PV_distance_mean_mean.B = nanmean(PV_diag_distance_mean.B,1);
    
    %get sem from each day-distance correlation
    day_PV_distance_mean_sem.A = nanstd(PV_diag_distance_mean.A,0,1)./sqrt(sum(~isnan(PV_diag_distance_mean.A),1));
    day_PV_distance_mean_sem.B = nanstd(PV_diag_distance_mean.B,0,1)./sqrt(sum(~isnan(PV_diag_distance_mean.B),1));
end

%% For export
mean_PV.A = day_PV_diag_mean_mean_exp.A;
mean_PV.B = day_PV_diag_mean_mean_exp.B;

sem_PV.A = day_PV_diag_sem_mean_exp.A;
sem_PV.B = day_PV_diag_sem_mean_exp.B;

%data for statistical analysis
raw.day_PV_diag_mean.A = day_PV_diag_mean.exp.A;
raw.day_PV_diag_mean.B = day_PV_diag_mean.exp.B;

%export all cross-day correlation distance data if experiment type is short term
%recall
if exp_type == 2
    %mean
    mean_PV.all_corr.A = day_PV_distance_mean_mean.A;
    mean_PV.all_corr.B = day_PV_distance_mean_mean.B;
    
    %sem
    sem_PV.all_corr.A = day_PV_distance_mean_sem.A;
    sem_PV.all_corr.B = day_PV_distance_mean_sem.B;
    
    %means used to calculate the mean (export data for analysis)
    raw.all_coor.PV_diag_distance_mean.A = PV_diag_distance_mean.A;
    raw.all_coor.PV_diag_distance_mean.B = PV_diag_distance_mean.B;
    
end


end

