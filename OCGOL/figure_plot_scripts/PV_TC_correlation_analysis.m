function [outputArg1,outputArg2] = PV_TC_correlation_analysis(short_term_learn, short_term_recall, long_term_recall)

%% Number of animals in each category

nb_st_learn = size(short_term_learn.PV_TC_corr,2);
nb_st_recall = size(short_term_recall.PV_TC_corr,2);
nb_lt_recall = size(long_term_recall.PV_TC_corr,2);


%% Load in non-normalized PV correlation matrices for each class of experiments

for aa=1:nb_st_learn
    %short term learn
    st_learn.PVmat.A = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
    st_learn.PVmat.B = short_term_learn.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;
end

for aa=1:nb_st_recall
    %short term recall
    st_recall.PVmat.A = short_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
    st_recall.PVmat.B = short_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;
end

for aa=1:nb_lt_recall
    %long term recall
    lt_recall.PVmat.A = long_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.A;
    lt_recall.PVmat.B = long_term_recall.PV_TC_corr(aa).PV_TC_corr.PVcorr_all_ses_no_nan.B;    
end


%% Filter out low quality sessions (learning and LT recall)

for aa=1:nb_st_learn
    st_learn.PVmat.A
end



%% 




end

