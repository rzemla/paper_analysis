function [lme_stats] = two_way_rm_lme(data_in_1,data_in_2)

%number of time points
nb_timepoints = size(data_in_1,2);

%number of total animals
nb_animals = size(data_in_1,1) + size(data_in_2,1);

%correlation score vector
corr_score = [data_in_1; data_in_2];
corr_score = corr_score(:);

%subject data
subject = categorical(repmat(1:nb_animals,nb_timepoints,1))';
subject = subject(:);

%behavior_type vector(1 = learn; 2 = recall)
behavior_type = [ones(size(data_in_1,1),nb_timepoints); 2*ones(size(data_in_2,1),nb_timepoints)];
behavior_type = categorical(behavior_type(:)); 

%time vector
%define ordinal time rankining
for ii=1:nb_timepoints
    %each time point is a cell
   time_value(ii) = {num2str(ii)};
end

time_vec = repmat([1:nb_timepoints]',1,nb_animals)';
time_vec = categorical(categorical(time_vec(:)),time_value,'Ordinal',true);

%organize test data into table
lme_tbl = table(corr_score, subject, time_vec, behavior_type, 'VariableNames',{'corr', 'subject','time','behavior'});

%fit LME with subjects being random factor
lme = fitlme(lme_tbl,'corr ~ 1+ time*behavior + (1|subject)','FitMethod','REML','DummyVarCoding','effects');

%anova on lme
lme_stats = anova(lme,'DFMethod','satterthwaite');


end

