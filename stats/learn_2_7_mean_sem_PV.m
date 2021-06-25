function [mean_2_7,sem_2_7,nb_samp] = learn_2_7_mean_sem_PV(source_data)

%day x animal rc
d2_7 = source_data([2,7],:) ;

%r - animals matching, c - day (2 vs 7)
d2_7_trim = d2_7(:,sum(~isnan(d2_7),1) ==2)';
nb_samp = size(d2_7_trim,1);

%mean and sem
mean_2_7 = mean(d2_7_trim,1);
sem_2_7 = std(d2_7_trim,0,1)./sqrt(nb_samp);

end

