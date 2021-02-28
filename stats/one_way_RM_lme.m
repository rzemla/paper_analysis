function [lme1_stats] = one_way_RM_lme(data_mat,nb_animals, nb_time_points)

%vectorize input data values
data_mat = data_mat(:);

%subject/animal data
subject = categorical(repmat(1:nb_animals,nb_time_points,1))';
subject = subject(:);

%time vector (define time points)
%define ordinal time ranking
for ii=1:nb_time_points
    %each time point is a cell
   time_value(ii) = {num2str(ii)};
end

%define vector
time_vec = repmat([1:7]',1,6)';
time_vec = categorical(categorical(time_vec(:)),time_value,'Ordinal',true);

%create entry table for 
lme1_tbl = table(data_mat, subject, time_vec, 'VariableNames',{'corr', 'subject','time'});

%fit LME with subjects being random factor
lme1 = fitlme(lme1_tbl,'corr ~ 1+ time + (1|subject)','FitMethod','REML','DummyVarCoding','effects');

%anova on lme
lme1_stats = anova(lme1,'DFMethod','satterthwaite');

end

