function [ROI_idx_tuning_class] = return_tuned_classes(tunedLogical,pf_count_filtered_log)


%incoporate place field present and event filter
%SI
tuned_count_filt{1}(1) = size(find((tunedLogical.si.onlyA_tuned & pf_count_filtered_log(1,:)) ==1),2);
tuned_count_filt{1}(2) = size(find((tunedLogical.si.onlyB_tuned & pf_count_filtered_log(2,:)) ==1),2);
tuned_count_filt{1}(3) = size(find((tunedLogical.si.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1),2);
tuned_count_filt{1}(4) = size(pf_count_filtered_log,2) -  sum(tuned_count_filt{1});

%TS
tuned_count_filt{2}(1) = size(find((tunedLogical.ts.onlyA_tuned & pf_count_filtered_log(1,:)) ==1),2);
tuned_count_filt{2}(2) = size(find((tunedLogical.ts.onlyB_tuned & pf_count_filtered_log(2,:)) ==1),2);
tuned_count_filt{2}(3) = size(find((tunedLogical.ts.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1),2);
tuned_count_filt{2}(4) = size(pf_count_filtered_log,2) -  sum(tuned_count_filt{2});

%spatial information (not filtered aside from stat sig outcome)
tuned_count{1}(1) = size(find(tunedLogical.si.onlyA_tuned ==1),2);
tuned_count{1}(2) = size(find(tunedLogical.si.onlyB_tuned ==1),2);
tuned_count{1}(3) = size(find(tunedLogical.si.AandB_tuned ==1),2);
tuned_count{1}(4) = size(find(tunedLogical.si.neither ==1),2); 

%tuning specificity
tuned_count{2}(1) = size(find(tunedLogical.ts.onlyA_tuned ==1),2);
tuned_count{2}(2) = size(find(tunedLogical.ts.onlyB_tuned ==1),2);
tuned_count{2}(3) = size(find(tunedLogical.ts.AandB_tuned ==1),2);
tuned_count{2}(4) = size(find(tunedLogical.ts.neither ==1),2); 

%% Export indices corresponding to each class

%SI
ROI_idx_tuning_class.si.Aonly = find((tunedLogical.si.onlyA_tuned & pf_count_filtered_log(1,:)) ==1);
ROI_idx_tuning_class.si.Bonly = find((tunedLogical.si.onlyB_tuned & pf_count_filtered_log(2,:)) ==1);
ROI_idx_tuning_class.si.AB = find((tunedLogical.si.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1);

%create neither indexes (by means of exclusion)
neitherLog = zeros(1,sum(tuned_count{1}));
neitherLog(ROI_idx_tuning_class.si.Aonly) = 1;
neitherLog(ROI_idx_tuning_class.si.Bonly) = 1;
neitherLog(ROI_idx_tuning_class.si.AB) = 1;

%neither idxs
ROI_idx_tuning_class.si.N = find(neitherLog == 0);

%generate count and check that it sums to total nb of neurons
ROI_nb_class.si(1) = size(ROI_idx_tuning_class.si.Aonly,2);
ROI_nb_class.si(2) = size(ROI_idx_tuning_class.si.Bonly,2);
ROI_nb_class.si(3) = size(ROI_idx_tuning_class.si.AB,2);
ROI_nb_class.si(4) = size(ROI_idx_tuning_class.si.N,2);

%check if total equals nb of neurons as input
isequal(sum(ROI_nb_class.si),size(pf_count_filtered_log,2))

%TS
ROI_idx_tuning_class.ts.Aonly = find((tunedLogical.ts.onlyA_tuned & pf_count_filtered_log(1,:)) ==1);
ROI_idx_tuning_class.ts.Bonly = find((tunedLogical.ts.onlyB_tuned & pf_count_filtered_log(2,:)) ==1);
ROI_idx_tuning_class.ts.AB = find((tunedLogical.ts.AandB_tuned & (pf_count_filtered_log(1,:) & pf_count_filtered_log(2,:))) ==1);

%create neither indexes (by means of exclusion)
neitherLog = zeros(1,sum(tuned_count{2}));
neitherLog(ROI_idx_tuning_class.ts.Aonly) = 1;
neitherLog(ROI_idx_tuning_class.ts.Bonly) = 1;
neitherLog(ROI_idx_tuning_class.ts.AB) = 1;

%neither idxs
ROI_idx_tuning_class.ts.N = find(neitherLog == 0);

%generate count and check that it sums to total nb of neurons
ROI_nb_class.ts(1) = size(ROI_idx_tuning_class.ts.Aonly,2);
ROI_nb_class.ts(2) = size(ROI_idx_tuning_class.ts.Bonly,2);
ROI_nb_class.ts(3) = size(ROI_idx_tuning_class.ts.AB,2);
ROI_nb_class.ts(4) = size(ROI_idx_tuning_class.ts.N,2);

%check if total equals nb of neurons as input
isequal(sum(ROI_nb_class.ts),size(pf_count_filtered_log,2))


%% Generate the equivalent logicals for each class of neurons in SI and TS
%needed as input for some downstream functions

%generate 0 vectors equivalent to all neurons and populate 
%si or ts have same nb of neurons, so doesn't matter
blank_vector = false(1,sum(ROI_nb_class.si));

%SI
%create blanks and preallocate
%Aonly
ROI_idx_tuning_class.si.log.Aonly = blank_vector;
ROI_idx_tuning_class.si.log.Aonly(ROI_idx_tuning_class.si.Aonly) = 1;
%Bonly
ROI_idx_tuning_class.si.log.Bonly = blank_vector;
ROI_idx_tuning_class.si.log.Bonly(ROI_idx_tuning_class.si.Bonly) = 1;
%A&B only
ROI_idx_tuning_class.si.log.AB= blank_vector;
ROI_idx_tuning_class.si.log.AB(ROI_idx_tuning_class.si.AB) = 1;
%neither
ROI_idx_tuning_class.si.log.N = blank_vector;
ROI_idx_tuning_class.si.log.N(ROI_idx_tuning_class.si.N) = 1;

%TS
%create blanks and preallocate
%Aonly
ROI_idx_tuning_class.ts.log.Aonly = blank_vector;
ROI_idx_tuning_class.ts.log.Aonly(ROI_idx_tuning_class.ts.Aonly) = 1;
%Bonly
ROI_idx_tuning_class.ts.log.Bonly = blank_vector;
ROI_idx_tuning_class.ts.log.Bonly(ROI_idx_tuning_class.ts.Bonly) = 1;
%A&B only
ROI_idx_tuning_class.ts.log.AB= blank_vector;
ROI_idx_tuning_class.ts.log.AB(ROI_idx_tuning_class.ts.AB) = 1;
%neither
ROI_idx_tuning_class.ts.log.N = blank_vector;
ROI_idx_tuning_class.ts.log.N(ROI_idx_tuning_class.ts.N) = 1;


end

