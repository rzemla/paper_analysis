function [mean_TC,sem_TC,raw] = filter_convert_day_return_AB_correlation_across_time(exp_struct,excl_day_combined_day_nan,exp_type,day_range,tuning_type)

%% Define parameters

%define number of animals
nb_animal = sum(~cellfun(@isempty,excl_day_combined_day_nan(:,exp_type)));


%% Extract all day 2 day cross correlations for short-term recall animals for Figure 4F plots
%run only for short-term recall since assuming behavioral steady-state
%there

%if short term recall experiment
%if exp_type == 2
    %extract TC correlation matrices for each animal that are day matched
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
%                     A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
%                     crossday_TC_mat.exp(aa).A{dd1,dd2} = A_mat{ses2day_idx_d1,ses2day_idx_d2};
                    
                    %early AB matrix (d1)
                    early_mat = getfield(exp_struct.TC_corr_match{aa}.tc_corr_match,tuning_type);
                    crossday_TC_mat.exp(aa).AB_early{dd1,dd2} = early_mat.TC_corr_all_day{ses2day_idx_d1, ses2day_idx_d2}.AB_AB_early;
                    
                    %late AB matrix (relative day)
                    late_mat = getfield(exp_struct.TC_corr_match{aa}.tc_corr_match,tuning_type);
                    crossday_TC_mat.exp(aa).AB_later{dd1,dd2} = late_mat.TC_corr_all_day{ses2day_idx_d1, ses2day_idx_d2}.AB_AB_later;                    
                     
                    
                    %for B trials
%                     B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
%                     crossday_TC_mat.exp(aa).B{dd1,dd2} = B_mat{ses2day_idx_d1,ses2day_idx_d2};
                    
                else
                    crossday_TC_mat.exp(aa).AB_early{dd1,dd2} = [];
                    crossday_TC_mat.exp(aa).AB_later{dd1,dd2} = [];
                end
                
            end
        end
    end
    
    %set empty cells in cross-day for each animal to nan
    %return cell indices where the cells are empty
    %for each animal
    for aa=1:nb_animal
        
        %for early correlation
        empty_log = cellfun(@isempty,crossday_TC_mat.exp(aa).AB_early);
        %set the empty cells to nan
        crossday_TC_mat.exp(aa).AB_early(empty_log) = {nan};
        
        %for late correlation
        empty_log = cellfun(@isempty,crossday_TC_mat.exp(aa).AB_later);
        %set the empty cells to nan
        crossday_TC_mat.exp(aa).AB_later(empty_log) = {nan};
    end
    
    %if long term recall experiment, compress matrix to reflect neighboring
    %session correlation rather than neighboring days (which does not make
    %sense)
    if exp_type == 3
       %for each animal
       for aa=1:nb_animal
           for dd1 = 1:size(day_range,2)
               for dd2 = 1:size(day_range,2)
                   %for A trials
                   crossday_compress_TC_mat.exp(aa).AB_early{dd1,dd2} = crossday_TC_mat.exp(aa).AB_early{day_range(dd1),day_range(dd2)};
                   %crossday_compress_TC_mat.exp(aa).A{dd1,dd2} = crossday_TC_mat.exp(aa).A{dd1,dd2};
                   %for B trials
                   crossday_compress_TC_mat.exp(aa).AB_later{dd1,dd2} = crossday_TC_mat.exp(aa).AB_later{day_range(dd1),day_range(dd2)};
                   %crossday_compress_TC_mat.exp(aa).B{dd1,dd2} = crossday_TC_mat.exp(aa).B{dd1,dd2};
               end
           end
       end
    
    %overwrite crossday_TC_mat with crossday_compress_TC_mat going forward
    crossday_TC_mat = crossday_compress_TC_mat;
       
    end
    
    
    %extract each diagnonal and upward for each day-by-day match
    %for each animal
    for aa=1:nb_animal
        %for each off diagonal
        %for early
        crossday_TC_main_diag(aa).AB_early = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).AB_early,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_TC_main_diag_mean(aa).AB_early = cell2mat(cellfun(@nanmean,crossday_TC_main_diag(aa).AB_early,'UniformOutput',false));
        %for later
        crossday_TC_main_diag(aa).AB_later = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).AB_later,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_TC_main_diag_mean(aa).AB_later = cell2mat(cellfun(@nanmean,crossday_TC_main_diag(aa).AB_later,'UniformOutput',false));
        
    end
    
    %divide each neuron by the ratio and then take median for EACH ANIMAL
    for aa=1:nb_animal
        %take first row of AB early
        crossday_TC_rel_D1(aa).AB_combined(1,:) = crossday_TC_main_diag(aa).AB_early(1,:);
        %take first row of AB later
        crossday_TC_rel_D1(aa).AB_combined(2,:) = crossday_TC_main_diag(aa).AB_later(1,:);
        
        %ratio of correlations for each neurons in each animal
        for dd=1:size(day_range,2)
            %later divided by early (normalized to early)
            corr_ratio_AB_combined{aa,dd} = crossday_TC_rel_D1(aa).AB_combined{2,dd}./crossday_TC_rel_D1(aa).AB_combined{1,dd};
        end
    end
    
    %get median of ratio of correlations for all neurons for each animal
    %use median b/c mean significantly skews data when high values are
    %generated from the ratio
    AB_corr_ratio_by_neuron_med = cellfun(@nanmedian, corr_ratio_AB_combined);
    
    %pool all ratio'ed individually neurons by day
    %for each day - bundle neurons together from all animals
    for dd=1:size(corr_ratio_AB_combined,2)
        pooled_ratioed_neurons{dd} = cell2mat(corr_ratio_AB_combined(:,dd));
    end
    
    %take median of the pooled ratios
    AB_corr_ratio_by_neuron_med_pooled = cellfun(@nanmedian, pooled_ratioed_neurons);
    
    %extract means along the 1st row (comparisons relative to D1)
    %animal x day
    for aa=1:nb_animal
        %for early
        AB_corr_TC_mean.AB_early(aa,:) = crossday_TC_main_diag_mean(aa).AB_early(1,:);
        %for later
        AB_corr_TC_mean.AB_later(aa,:) = crossday_TC_main_diag_mean(aa).AB_later(1,:);
    end

%block out neuron by neuron analysis
%{    
    %% Neuron-by-neuron analysis - merge all correlation scores from each - skip for now (return later if necessary)
    %neuron together into 1 cell
        %extract each diagnonal and upward for each day-by-day match
    %for each animal
    for aa=1:nb_animal
        %for each off diagonal
        %for A
        crossday_TC_main_diag_neuron(aa).A = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).A,'UniformOutput',false);

        %for B
        crossday_TC_main_diag_neuron(aa).B = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).B,'UniformOutput',false);
        
    end
    
    %place each ii,jj day for each animal into subcells and merge with
    %subcells (easier than dealing with 3D matrix)
    for aa=1:nb_animal
        %size of matrix for A or B trials is the same
        for ii = 1:size(crossday_TC_main_diag_neuron(aa).A,2)
            for jj = 1:size(crossday_TC_main_diag_neuron(aa).A,2)
                %for A trials
                crossday_TC_main_diag_neuron_cell_merge.A{ii,jj}{aa} = crossday_TC_main_diag_neuron(aa).A{ii,jj};
                %for B trials
                crossday_TC_main_diag_neuron_cell_merge.B{ii,jj}{aa} = crossday_TC_main_diag_neuron(aa).B{ii,jj};
            end
        end
        
    end
    
    %combine correlation values from each animal across days
    %size of matrix for A or B trials is the same
    for ii = 1:size(crossday_TC_main_diag_neuron(aa).A,2)
        for jj = 1:size(crossday_TC_main_diag_neuron(aa).A,2)
            %merge into one cell each cell (all TC correlations of neurons by
            %day combination
            %for A trials
            crossday_TC_main_diag_neuron_merge.A{ii,jj} = cell2mat(crossday_TC_main_diag_neuron_cell_merge.A{ii,jj}');
            %for B trials
            crossday_TC_main_diag_neuron_merge.B{ii,jj} = cell2mat(crossday_TC_main_diag_neuron_cell_merge.B{ii,jj}');
        end
    end
    
    %take mean of each merge matrix from all animals
    crossday_TC_main_diag_neuron_merge_mean.A = cellfun(@nanmean,crossday_TC_main_diag_neuron_merge.A);
    crossday_TC_main_diag_neuron_merge_mean.B = cellfun(@nanmean,crossday_TC_main_diag_neuron_merge.B);
    
    
    %extract the distance diagonals
    for d_idx=0:(size(crossday_TC_main_diag_neuron_merge_mean.A,1) -1)
        %take each diagnonal starting from main
        %for A trials
        TC_diag_distance_neuron.A{d_idx+1} = diag(crossday_TC_main_diag_neuron_merge_mean.A,d_idx);
        %for B trials
        TC_diag_distance_neuron.B{d_idx+1} = diag(crossday_TC_main_diag_neuron_merge_mean.B,d_idx);
    end
    
    
    %merge all neurons together from the diagonal days 

    %this approach is used to extract cell diagonals using linear indices
    %(if merging of all neurons from distance-day correlations is desired)
    %take diagonal of each cell and merge into distance matrix
    
    %generate the linear indices for the cell diagonals
    idxs_for_diag_cell_extract = reshape(1:81,9,9);
    
    for d_idx=0:(size(crossday_TC_main_diag_neuron_merge_mean.A,1) -1)
        %get linear indices of the desired diagnonal
        idx_temp = diag(idxs_for_diag_cell_extract,d_idx);
        %for A trials
        TC_diag_distance_all_neuron.A{d_idx+1} = cell2mat(crossday_TC_main_diag_neuron_merge.A(idx_temp));
        %for B trials
        TC_diag_distance_all_neuron.B{d_idx+1} = cell2mat(crossday_TC_main_diag_neuron_merge.B(idx_temp));
    end
    
%}
%end

%% Get ratio as well as individual values

%mean of later correlation / mean of early correlation for each animal -
%mean calculated on later and early correlation (not matched by neurons as
%below)
AB_corr_ratio = AB_corr_TC_mean.AB_later./AB_corr_TC_mean.AB_early;

%median of ratio of correlations for each neuron
AB_corr_ratio_by_neuron_med;

%% QC plot
figure
subplot(1,3,1)
hold on
plot(nanmean(AB_corr_ratio,1))

subplot(1,3,2)
hold on
plot(nanmean(AB_corr_ratio_by_neuron_med,1))

subplot(1,3,3)
hold on
title('Median for pooled neurons');
plot(AB_corr_ratio_by_neuron_med_pooled)

%% Mean and SEM of TCs + 95% CI around the median for pooled neurons

%mean ratio of means
AB_corr_ratio_mean_mean = nanmean(AB_corr_ratio,1);

%get sem across animals using correct count of animals per day
AB_corr_ratio_mean_sem = nanstd(AB_corr_ratio,0,1)./sqrt(sum(~isnan(AB_corr_ratio),1));

%get mean of median of neurons
AB_corr_ratio_by_neuron_med_mean = nanmean(AB_corr_ratio_by_neuron_med,1);

%get sem
AB_corr_ratio_by_neuron_med_sem = nanstd(AB_corr_ratio_by_neuron_med,0,1)./sqrt(sum(~isnan(AB_corr_ratio_by_neuron_med),1));

%get 95% CI around median from pooled analysis
AB_corr_ratio_pooled_95ci = 1.57.*((cellfun(@iqr, pooled_ratioed_neurons))./sqrt(cellfun(@(x) size(x,1), pooled_ratioed_neurons)));

%% QC - Plot median and 95% CI 
% figure
% hold on
% errorbar(1:9,AB_corr_ratio_by_neuron_med_pooled,AB_corr_ratio_pooled_95ci )

%% For export
%by animals stats
%mean
mean_TC.animal.AB_corr_ratio = AB_corr_ratio_mean_mean;
mean_TC.neuron.AB_corr_ratio = AB_corr_ratio_by_neuron_med_mean;
%for pooled neurons, use median
mean_TC.pooled.AB_corr_ratio_median = AB_corr_ratio_by_neuron_med_pooled;

%sem
sem_TC.animal.AB_corr_ratio = AB_corr_ratio_mean_sem;
sem_TC.neuron.AB_corr_ratio = AB_corr_ratio_by_neuron_med_sem;
%for pooled neurons, use 95% CI
sem_TC.pooled.AB_corr_ratio_95ci = AB_corr_ratio_pooled_95ci;

%mean per animal used to make plots (export to Prism)
raw.animal.AB_corr_ratio = AB_corr_ratio;
raw.neuron.AB_corr_ratio = AB_corr_ratio_by_neuron_med;
raw.pooled.AB_corr_ratio = pooled_ratioed_neurons;


%{

%% Cross-correlation distance data - get mean and sem 
%run only for short term recall data

if exp_type == 2
    %get mean from each day-distance correlation
    day_TC_distance_mean_mean.A = nanmean(TC_diag_distance_mean.A,1);
    day_TC_distance_mean_mean.B = nanmean(TC_diag_distance_mean.B,1);
    
    %get sem from each day-distance correlation
    day_TC_distance_mean_sem.A = nanstd(TC_diag_distance_mean.A,0,1)./sqrt(sum(~isnan(TC_diag_distance_mean.A),1));
    day_TC_distance_mean_sem.B = nanstd(TC_diag_distance_mean.B,0,1)./sqrt(sum(~isnan(TC_diag_distance_mean.B),1));
    
    %get mean from each day mean correlation from merging all neurons from
    %all animals
    TC_diag_distance_neuron_mean.A = cellfun(@nanmean,TC_diag_distance_neuron.A);
    TC_diag_distance_neuron_mean.B = cellfun(@nanmean,TC_diag_distance_neuron.B);
    
    %get sem
    %A trials
    TC_diag_distance_neuron_sem.A = cellfun(@nanstd,TC_diag_distance_neuron.A)./...   
    sqrt(cellfun(@sum,(cellfun(@(x) ~isnan(x),TC_diag_distance_neuron.A,'UniformOutput',false))));

    %B trials
    TC_diag_distance_neuron_sem.B = cellfun(@nanstd,TC_diag_distance_neuron.B)./...   
    sqrt(cellfun(@sum,(cellfun(@(x) ~isnan(x),TC_diag_distance_neuron.B,'UniformOutput',false))));    
    
    %get mean from each day from merging all neurons from all animals from
    %all distance day correlations
    TC_diag_distance_all_neuron_mean.A = cellfun(@nanmean,TC_diag_distance_all_neuron.A);
    TC_diag_distance_all_neuron_mean.B = cellfun(@nanmean,TC_diag_distance_all_neuron.B);
    
    %get sem
    %A
    TC_diag_distance_all_neuron_sem.A = cellfun(@nanstd,TC_diag_distance_all_neuron.A)./...
    sqrt(cellfun(@sum,(cellfun(@(x) ~isnan(x),TC_diag_distance_all_neuron.A,'UniformOutput',false))));
    %B
    TC_diag_distance_all_neuron_sem.B = cellfun(@nanstd,TC_diag_distance_all_neuron.B)./...
    sqrt(cellfun(@sum,(cellfun(@(x) ~isnan(x),TC_diag_distance_all_neuron.B,'UniformOutput',false))));    
     
end

%% Combine all neurons together and calculate the TC correlation score on the population

%bundle all the diagonal values (individual corrlations for all neurons
%from each day (to get mean and sem across neurons)

%for each day, merge the tuning curve correlation values
for dd=1:size(day_TC_diag.exp.A,2)
    %for matching A trials
    day_TC_neurons_combined.A{dd} = cat(1,day_TC_diag.exp.A{:,dd});
    %for matching B trials
    day_TC_neurons_combined.B{dd} = cat(1,day_TC_diag.exp.B{:,dd});
end

%get mean for pooled neuron correlation (there shouldn't be any nans)
day_TC_neuron_mean.exp.A = cellfun(@nanmean,day_TC_neurons_combined.A);
day_TC_neuron_mean.exp.B = cellfun(@nanmean,day_TC_neurons_combined.B);

%get sem for pooled neuron correlation (there shouldn't be any nans)
day_TC_neuron_sem.exp.A = cellfun(@nanstd,day_TC_neurons_combined.A)./...
sqrt(cellfun(@(x) size(x,1), day_TC_neurons_combined.A));

day_TC_neuron_sem.exp.B = cellfun(@nanstd,day_TC_neurons_combined.B)./...
sqrt(cellfun(@(x) size(x,1), day_TC_neurons_combined.B));

%QC check that the number of neurons for correlations adds to neurons
%exported previously - should add up since no nan values

%QC check this - checks out
%isequal(cell2mat(day_TC_diag.exp.A(:,5)), day_TC_neurons_combined.A{5})

%% For export from cross-cross and rel D1 analysis

%by neuron stats
%mean
mean_TC.neuron.A = day_TC_neuron_mean.exp.A;
mean_TC.neuron.B = day_TC_neuron_mean.exp.B;

%sem
sem_TC.neuron.A = day_TC_neuron_sem.exp.A;
sem_TC.neuron.B = day_TC_neuron_sem.exp.B;


%correlation of each neuron (export to Prism)
raw.neuron.day_TC_neurons_combined = day_TC_neurons_combined;


%export all cross-day correlation distance data if experiment type is short term
%recall
if exp_type == 2
    %mean
    mean_TC.all_corr.A = day_TC_distance_mean_mean.A;
    mean_TC.all_corr.B = day_TC_distance_mean_mean.B;
    
    %sem
    sem_TC.all_corr.A = day_TC_distance_mean_sem.A;
    sem_TC.all_corr.B = day_TC_distance_mean_sem.B;
    
    %means used to calculate the mean (export data for analysis)
    raw.all_corr.TC_diag_distance_mean.A = TC_diag_distance_mean.A;
    raw.all_corr.TC_diag_distance_mean.B = TC_diag_distance_mean.B;
    
    %export neuron merge mean and sem data
    %cumulative mean from merging all animal data into distance-day
    %analysis
    
    %%% means from all neuron bundled by animals across days
    %each point is a mean taken for that day from the neurons of all
    %animals
    
    %mean
    mean_TC.all_corr.neuron_mean_day.A = TC_diag_distance_neuron_mean.A;
    mean_TC.all_corr.neuron_mean_day.B = TC_diag_distance_neuron_mean.B;
    
    %sem
    sem_TC.all_corr.neuron_mean_day.A = TC_diag_distance_neuron_sem.A;
    sem_TC.all_corr.neuron_mean_day.B = TC_diag_distance_neuron_sem.B;
    
    %raw
    raw.all_corr.neuron_mean_day.TC_diag_distance_neuron.A = TC_diag_distance_neuron.A;
    raw.all_corr.neuron_mean_day.TC_diag_distance_neuron.B = TC_diag_distance_neuron.B;
    
    %%% all neurons taken together
    
    %mean
    mean_TC.all_corr.neuron_ind_all.A = TC_diag_distance_all_neuron_mean.A;
    mean_TC.all_corr.neuron_ind_all.B = TC_diag_distance_all_neuron_mean.B;
    %sem
    sem_TC.all_corr.neuron_ind_all.A = TC_diag_distance_all_neuron_sem.A;
    sem_TC.all_corr.neuron_ind_all.B = TC_diag_distance_all_neuron_sem.B;
    
    %raw
    raw.all_corr.neuron_ind_all.TC_diag_distance_all_neuron.A = TC_diag_distance_all_neuron.A;
    raw.all_corr.neuron_ind_all.TC_diag_distance_all_neuron.B = TC_diag_distance_all_neuron.B;
    
    %DONE - make sure to label variables in Word document
    
end
%}



end

