function [mean_out, sem_out, nb_animals] = neighbor_corr_mean_sem(nei_mat)

%trim matrix
nei_mat_trim = nei_mat(sum(~isnan(nei_mat(:, [1,6])),2) ==2,[1,6]);

%nb_animals
nb_animals = size(nei_mat_trim,1);
%mean
mean_out = mean(nei_mat_trim,1);
%sem
sem_out = std(nei_mat_trim,0,1)./sqrt(nb_animals);


end

