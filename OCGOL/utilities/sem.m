function [sem_output] = sem(data_matrix)
%calculate standard error of mean (no nans) and assumes each treatment has
%same amount of subjects
%input: data matrix -
        %rows = animals/subjects
        %columns = treatment groups 
%output: vector of sem at each points


sem_output = std(data_matrix,0,1)./sqrt(size(data_matrix,1));


end

