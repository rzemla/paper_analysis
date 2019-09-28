function [outputArg1,outputArg2] = detect_SCE_assembly(sce_activity_matrix,options)
%% Detects SCEs associated with reccurent cell assemblies (Fig. 2 Malvache 2016)

%% Matrix ROI participation vs. SCE (row x column)
%all active ROISs (input ROIs - manually selected across entire session)
assembly_mat = sce_activity_matrix;

%% Normalized covariance matrix and squared Euclidean distance between columns

%normalized covariance matrix on SCEs
norm_cov_SCE = corr(assembly_mat,'Rows','complete');
%covariance matrix on SCEs
cov_SCE = cov(assembly_mat);

%% Remove NaN rows from normalized covariance matrix

%create copy
norm_cov_SCE_nonan = norm_cov_SCE;
%remove rows and columns that are nan
remove_norm_cov_idx = find(isnan(norm_cov_SCE_nonan(1,:)) ==1);
norm_cov_SCE_nonan(remove_norm_cov_idx,:) =[];
norm_cov_SCE_nonan(:,remove_norm_cov_idx) =[];

%replace no-nan'd version as the input
norm_cov_SCE = norm_cov_SCE_nonan;

figure
imagesc(cov_SCE)
hold on
colorbar

%calculated squared euclidean distance between columns of norm cov matrix
sqE = pdist2(norm_cov_SCE, norm_cov_SCE,'squaredeuclidean');

%% Run k-mean clustering on columns of square Euclidean metric
%100 iterations
%run for 2-19 clusters
%idx - SCE cluster assignments
%parfor ii=2:19
clust_max = options.clust_max;

parfor ii=2:clust_max
    disp(['# clusters: ', num2str(ii)])
    for iter_nb=1:100
        disp(['Iteration nb: ', num2str(iter_nb)])
        [cluster_assign{ii-1}{iter_nb}, ~] = kmeans(sqE,ii,'MaxIter',iter_nb);
    end
end

%% Get sihouette value
%calculate the silhoutte value for each cluster, take max
parfor cc=1:(clust_max-1)
    disp(cc);
    for iter_nb=1:100
        %mean silhouette value for each iteration
        sil_vals{cc}{iter_nb} = silhouette(norm_cov_SCE,cluster_assign{cc}{iter_nb});
    end
end

%get average silhoutte value of each cluster and each iteration (for
%sorting)
%for each cluster
for cc=1:(clust_max-1)
    %take cluster assignment for that iteration of k-means
    for iter_nb=1:100
        %for each set of cluster assign indices
        for ii=1:max(cluster_assign{cc}{iter_nb})
            %take mean sil value for that cluster
            mean_clust_score{cc}{iter_nb}(ii) = mean(sil_vals{cc}{iter_nb}(find(cluster_assign{cc}{iter_nb} == ii)));
        end
    end
end

%mean of silouette values
for cc=1:(clust_max-1)
    mean_sil_vals{cc} = cellfun(@mean,sil_vals{cc},'UniformOutput',true);
end

figure
boxplot(cell2mat(mean_sil_vals')')

%figure out cluster number to use cluster assignment based in iter of
%k-means
[~, nb_clust] = max(cellfun(@max,mean_sil_vals));
%get the itertion of k-means where max occurred 
[~,select_k_iter] = max(mean_sil_vals{nb_clust})

%get cluster sort order by mean sihouette value - descending order
[B,I_sort] = sort(mean_clust_score{nb_clust}{select_k_iter},'descend') 

%generate cluster sort order
%numerical order
idx_order = size(cov_SCE,1);
%create sorted cluster
idx_cluster_order = [];
%nb_clust = 25;
%select_k_iter = 11;

%make cluster order index vector
linear_idx = 1;
for ii=I_sort
    cluster_SCEs_nb{linear_idx} = find(cluster_assign{nb_clust}{select_k_iter} ==ii);
    idx_cluster_order = [idx_cluster_order;find(cluster_assign{nb_clust}{select_k_iter} ==ii)];
    %update linear idx
    linear_idx = linear_idx + 1;
end

%2 cluster reordered matrix
norm_cov_SCE_clust= norm_cov_SCE(idx_cluster_order, idx_cluster_order);

cov_SCE_clust =  cov_SCE(idx_cluster_order, idx_cluster_order);

%% Plot norm covariance and covariance
figure('Position',[2058 546 1029 420])
subplot(1,2,1)
imagesc(norm_cov_SCE_clust)
hold on
title('Normalized covariance')
axis square
colormap('jet')
colorbar
subplot(1,2,2)
imagesc(cov_SCE_clust)
hold on
title('Covariance')
colorbar
axis square

%% Extract SCE's within cluster

% %sort assembly matrix
% assembly_mat_sort = assembly_mat(:,idx_cluster_order);
% %take out subset of each all SCE assemblies
% assembly_subset =  assembly_mat(:,cluster_SCEs_nb{4});
% 
% %sort by number of matches
% [~,I_subset_sort] = sort(sum(assembly_subset,2),'descend');
% assembly_subset_sort = assembly_subset(I_subset_sort,:);
% 
% figure;
% imagesc(assembly_subset_sort)
% 
% %plot cluster sort assembly matrix
% figure
% imagesc(assembly_mat_sort)

end

