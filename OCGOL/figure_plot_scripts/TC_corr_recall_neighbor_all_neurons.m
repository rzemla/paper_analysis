function [mean_TC_recall_neighbor_all,sem_TC_recall_neighbor_all] = TC_corr_recall_neighbor_all_neurons(neighbor_STC_recall_days_all)

%% Calculate TC correlation neighbor recall

%correlation coefcient for each neuron
for dd=1:8
    %check if empty
    if ~isempty(neighbor_STC_recall_days_all.A{dd})
        %A
        TC_neighbor_STC_recall_days_all.A{dd} = diag(corr(neighbor_STC_recall_days_all.A{dd}(:,1:100)',neighbor_STC_recall_days_all.A{dd}(:,101:200)'));
    else
        TC_neighbor_STC_recall_days_all.A{dd} = [];
    end
    
    if ~isempty(neighbor_STC_recall_days_all.B{dd})
        %B
        TC_neighbor_STC_recall_days_all.B{dd} = diag(corr(neighbor_STC_recall_days_all.B{dd}(:,1:100)',neighbor_STC_recall_days_all.B{dd}(:,101:200)'));
    else
        TC_neighbor_STC_recall_days_all.B{dd} = [];
    end
end

%get mean for each day (across all animals)
%A
mean_TC_recall_neighbor_all.A = cellfun(@nanmean,TC_neighbor_STC_recall_days_all.A);
%B
mean_TC_recall_neighbor_all.B = cellfun(@nanmean,TC_neighbor_STC_recall_days_all.B);

%% Get sem - neighbor RECALL 
%A
%can't calculate in one shot b/c 0 empty vectors returned by std func
std_TC_recall_neighbor_all.A = cellfun(@(x) std(x,0,1), TC_neighbor_STC_recall_days_all.A,'UniformOutput',false);
%days missing (to fill with nan for consistency of numbering
nb_missing = sum(cellfun(@isempty,std_TC_recall_neighbor_all.A));
%set empty cells to nan
std_TC_recall_neighbor_all.A(cellfun(@isempty,std_TC_recall_neighbor_all.A))= num2cell(nan(1,nb_missing));
%calculate sem
sem_TC_recall_neighbor_all.A = cell2mat(std_TC_recall_neighbor_all.A)./sqrt(cellfun(@(x) size(x,1),TC_neighbor_STC_recall_days_all.A,'UniformOutput',true));

%B
%can't calculate in one shot b/c 0 empty vectors returned by std func
std_TC_recall_neighbor_all.B = cellfun(@(x) std(x,0,1), TC_neighbor_STC_recall_days_all.B,'UniformOutput',false);
%days missing (to fill with nan for consistency of numbering
nb_missing = sum(cellfun(@isempty,std_TC_recall_neighbor_all.B));
%set empty cells to nan
std_TC_recall_neighbor_all.B(cellfun(@isempty,std_TC_recall_neighbor_all.B))= num2cell(nan(1,nb_missing));
%calculate sem
sem_TC_recall_neighbor_all.B = cell2mat(std_TC_recall_neighbor_all.B)./sqrt(cellfun(@(x) size(x,1),TC_neighbor_STC_recall_days_all.B,'UniformOutput',true));


end

