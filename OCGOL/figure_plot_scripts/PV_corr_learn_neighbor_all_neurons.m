function [mean_mean_PV_recall_neighbor_all,sem_PV_recall_neighbor_all] = PV_corr_learn_neighbor_all_neurons(neighbor_STC_recall_days_all)

%% Calculate TC correlation neighbor recall

%correlation coefcient for each neuron
for dd=1:8
    %for every animal
    for aa=1:size(neighbor_STC_recall_days_all.A,1)
        %check if empty
        if ~isempty(neighbor_STC_recall_days_all.A{aa,dd})
            %A
            PV_neighbor_STC_recall_days_all.A{aa,dd} = diag(corr(neighbor_STC_recall_days_all.A{aa,dd}(:,1:100),neighbor_STC_recall_days_all.A{aa,dd}(:,101:200)));
        else
            PV_neighbor_STC_recall_days_all.A{aa,dd} = [];
        end
        
        if ~isempty(neighbor_STC_recall_days_all.B{aa,dd})
            %B
            PV_neighbor_STC_recall_days_all.B{aa,dd} = diag(corr(neighbor_STC_recall_days_all.B{aa,dd}(:,1:100),neighbor_STC_recall_days_all.B{aa,dd}(:,101:200)));
        else
            PV_neighbor_STC_recall_days_all.B{aa,dd} = [];
        end
    end
end

%get mean for each day - for each animal
%A
mean_PV_recall_neighbor_each.A = cellfun(@nanmean,PV_neighbor_STC_recall_days_all.A);
%B
mean_PV_recall_neighbor_each.B = cellfun(@nanmean,PV_neighbor_STC_recall_days_all.B);

%get mean of means (across animals)
mean_mean_PV_recall_neighbor_all.A  = nanmean(mean_PV_recall_neighbor_each.A,1);
mean_mean_PV_recall_neighbor_all.B  = nanmean(mean_PV_recall_neighbor_each.B,1);

%animal in each session
nb_animals_by_ses = sum(~isnan(mean_PV_recall_neighbor_each.A),1);

%% Get sem - neighbor RECALL 
%A
%can't calculate in one shot b/c 0 empty vectors returned by std func
sem_PV_recall_neighbor_all.A = nanstd(mean_PV_recall_neighbor_each.A,0,1)./sqrt(nb_animals_by_ses);

%B
sem_PV_recall_neighbor_all.B = nanstd(mean_PV_recall_neighbor_each.B,0,1)./sqrt(nb_animals_by_ses);


end

