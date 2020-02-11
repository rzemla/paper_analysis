function [mean_TC,sem_TC,raw] = filter_convert_day_return_mean_sem_TC_global(exp_struct,excl_day_combined_day_nan,exp_type,day_range,tuning_type)

%% Extract TC correlation matrices

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
            A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
            day_TC_mat.exp.A{aa,dd} = A_mat{1,ses2day_idx}; 
            
            %day_TC_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.tuning_type.A{1,ses2day_idx};
            %for B trials
            B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
            day_TC_mat.exp.B{aa,dd} = B_mat{1,ses2day_idx};
            
            %day_TC_mat.exp.B{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.ts.B{1,ses2day_idx};
        else
            day_TC_mat.exp.A{aa,dd} = [];
            day_TC_mat.exp.B{aa,dd} = [];
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
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 3);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 6);
                
            elseif dd==5
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 3);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 7);
                
            end
            %if no assignmenet for that day (make sure both days exist in this case)
            if size(ses2day_idx,2) == 2
                %animal x day TC matrix cell
                %for A trials
                %day_TC_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.ts.A{ses2day_idx(1),ses2day_idx(2)};
                A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
                day_TC_mat.exp.A{aa,dd} = A_mat{ses2day_idx(1),ses2day_idx(2)}; 
                
                %for B trials
                %day_TC_mat.exp.B{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.ts.B{ses2day_idx(1),ses2day_idx(2)};
                B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
                day_TC_mat.exp.B{aa,dd} = B_mat{ses2day_idx(1),ses2day_idx(2)};
                
            else
                day_TC_mat.exp.A{aa,dd} = [];
                day_TC_mat.exp.B{aa,dd} = [];
                
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
                    A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
                    crossday_TC_mat.exp(aa).A{dd1,dd2} = A_mat{ses2day_idx_d1,ses2day_idx_d2};
                    
                    %for B trials
                    B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
                    crossday_TC_mat.exp(aa).B{dd1,dd2} = B_mat{ses2day_idx_d1,ses2day_idx_d2};
                    
                else
                    crossday_TC_mat.exp(aa).A{dd1,dd2} = [];
                    crossday_TC_mat.exp(aa).B{dd1,dd2} = [];
                end
                
            end
        end
    end
    
    %set empty cells in cross-day for each animal to nan
    %return cell indices where the cells are empty
    %for each animal
    for aa=1:nb_animal
        
        %for A trials
        empty_log = cellfun(@isempty,crossday_TC_mat.exp(aa).A);
        %set the empty cells to nan
        crossday_TC_mat.exp(aa).A(empty_log) = {nan};
        
        %for B trials
        empty_log = cellfun(@isempty,crossday_TC_mat.exp(aa).B);
        %set the empty cells to nan
        crossday_TC_mat.exp(aa).B(empty_log) = {nan};
    end
    
    
    %extract each diagnonal and upward for each day-by-day match
    %for each animal
    for aa=1:nb_animal
        %for each off diagonal
        %for A
        crossday_TC_main_diag(aa).A = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).A,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_TC_main_diag_mean(aa).A = cell2mat(cellfun(@nanmean,crossday_TC_main_diag(aa).A,'UniformOutput',false));
        %for B
        crossday_TC_main_diag(aa).B = cellfun(@(x) diag(x,0),crossday_TC_mat.exp(aa).B,'UniformOutput',false);
        %take mean (convert to matrix)
        crossday_TC_main_diag_mean(aa).B = cell2mat(cellfun(@nanmean,crossday_TC_main_diag(aa).B,'UniformOutput',false));
        
    end
    
    for aa=1:nb_animal
        for d_idx=0:(size(crossday_TC_main_diag_mean(aa).A,1) -1)
            %take each diagnonal starting from main
            %for A trials
            TC_diag_distance(aa).A{d_idx+1} = diag(crossday_TC_main_diag_mean(aa).A,d_idx);
            %for B trials
            TC_diag_distance(aa).B{d_idx+1} = diag(crossday_TC_main_diag_mean(aa).B,d_idx);
        end
        %take mean of each diagonal day distance for A (animal x day
        %separation)
        %A trials
        TC_diag_distance_mean.A(aa,:) = cellfun(@nanmean, TC_diag_distance(aa).A);
        %B trials
        TC_diag_distance_mean.B(aa,:) = cellfun(@nanmean, TC_diag_distance(aa).B);
    end
    
    
    %% Neuron-by-neuron analysis - merge all correlation scores from each
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
    
    
    
end

%% Get diagnonal and mean of TC for each day (by animal)

%get diagnonal for A/B
day_TC_diag.exp.A = cellfun(@diag,day_TC_mat.exp.A,'UniformOutput',false);
day_TC_diag.exp.B = cellfun(@diag,day_TC_mat.exp.B,'UniformOutput',false);

%get mean of diagnonal for A
day_TC_diag_mean.exp.A = cell2mat(cellfun(@nanmean,day_TC_diag.exp.A,'UniformOutput',false));
day_TC_diag_mean.exp.B = cell2mat(cellfun(@nanmean,day_TC_diag.exp.B,'UniformOutput',false));

%get mean across animals of the mean of the diagnonals
day_TC_diag_mean_mean_exp.A = nanmean(day_TC_diag_mean.exp.A,1);
day_TC_diag_mean_mean_exp.B = nanmean(day_TC_diag_mean.exp.B,1);

%get sem across animals using correct count of animals per day
day_TC_diag_sem_mean_exp.A = nanstd(day_TC_diag_mean.exp.A,0,1)./sqrt(sum(~isnan(day_TC_diag_mean.exp.A),1));
day_TC_diag_sem_mean_exp.B = nanstd(day_TC_diag_mean.exp.B,0,1)./sqrt(sum(~isnan(day_TC_diag_mean.exp.B),1));

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

%% For export
%by animals stats
%mean
mean_TC.animal.A = day_TC_diag_mean_mean_exp.A;
mean_TC.animal.B = day_TC_diag_mean_mean_exp.B;
%sem
sem_TC.animal.A = day_TC_diag_sem_mean_exp.A;
sem_TC.animal.B = day_TC_diag_sem_mean_exp.B;

%by neuron stats
%mean
mean_TC.neuron.A = day_TC_neuron_mean.exp.A;
mean_TC.neuron.B = day_TC_neuron_mean.exp.B;
%sem
sem_TC.neuron.A = day_TC_neuron_sem.exp.A;
sem_TC.neuron.B = day_TC_neuron_sem.exp.B;

%mean per animal used to make plots (export to Prism)
raw.animal.day_TC_diag_mean = day_TC_diag_mean;
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




end

