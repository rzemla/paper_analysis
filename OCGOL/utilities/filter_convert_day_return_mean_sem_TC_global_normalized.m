function [mean_TC,sem_TC] = filter_convert_day_return_mean_sem_TC_global_normalized(exp_struct,excl_day_combined_day_nan,exp_type,day_range,tuning_type)

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
            %A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
            %day_TC_mat.exp.A{aa,dd} = A_mat{1,ses2day_idx}; 
            
            %get normalized TC correlation matrices
            norm_corr_mat = getfield(exp_struct.TC_corr_match{aa}.tc_corr_match,tuning_type,'TC_corr_all_day_nonsorted');
            
            day_TC_mat.exp.A{aa,dd} = norm_corr_mat{1,ses2day_idx}.A;
            
            %day_TC_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.tuning_type.A{1,ses2day_idx};
            %for B trials
            %B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
            day_TC_mat.exp.B{aa,dd} = norm_corr_mat{1,ses2day_idx}.B;
            
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
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 2);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 6);
                
            elseif dd==5
                ses2day_idx(1) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 2);
                ses2day_idx(2) = find(excl_day_combined_day_nan{aa,exp_type}(2,:) == 7);
                
            end
            %if no assignmenet for that day (make sure both days exist in this case)
            if size(ses2day_idx,2) == 2
                %animal x day TC matrix cell
                %for A trials
                %day_TC_mat.exp.A{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.ts.A{ses2day_idx(1),ses2day_idx(2)};
                %get normalized matrices
                norm_corr_mat = getfield(exp_struct.TC_corr_match{aa}.tc_corr_match,tuning_type,'TC_corr_all_day_nonsorted');
            
                %A_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'A');
                day_TC_mat.exp.A{aa,dd} = norm_corr_mat{ses2day_idx(1),ses2day_idx(2)}.A;
                
                
                %for B trials
                %day_TC_mat.exp.B{aa,dd} = exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses.ts.B{ses2day_idx(1),ses2day_idx(2)};
                %B_mat = getfield(exp_struct.PV_TC_corr(aa).PV_TC_corr.TCcorr_all_ses,tuning_type,'B');
                %day_TC_mat.exp.B{aa,dd} = B_mat{ses2day_idx(1),ses2day_idx(2)};
                
                day_TC_mat.exp.B{aa,dd} = norm_corr_mat{ses2day_idx(1),ses2day_idx(2)}.B;
                
            else
                day_TC_mat.exp.A{aa,dd} = [];
                day_TC_mat.exp.B{aa,dd} = [];
                
            end
        end
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


end

