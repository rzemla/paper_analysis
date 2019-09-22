%% Detects SCEs associated with reccurent cell assemblies (Fig. 2 Malvache 2016)

%% Rough SCE selection

select_SCE = find(seq_overlap>5);


%% Create matrix ROI vs. SCE (row x column)
%all active ROISs (input ROIs - manually selected across entire session)
assembly_mat = zeros(size(pca_input,2),size(select_SCE,2));

%fill in the cells that participated in each sce
for ss=1:size(select_SCE,2)
    assembly_mat(SCE_ROIs{select_SCE(ss)},ss) = 1;
end

%% Normalized covariance matrix and squared Euclidean distance between columns
assembly_mat = sce_activity_matrix;

%normalized covariance matrix
%norm_cov_SCE = corrcoef(assembly_mat);
norm_cov_SCE = corr(assembly_mat);
cov_SCE = cov(assembly_mat);

%calculated squared euclidean distance between columns of norm cov matrix
sqE = pdist2(norm_cov_SCE, norm_cov_SCE,'squaredeuclidean');

%% Run k-mean clustering on columns of square Euclidean metric
%100 iterations
%run for 2-19 clusters
%idx - SCE cluster assignments
%parfor ii=2:19
parfor ii=2:19
    disp(['# clusters: ', num2str(ii)])
    for iter_nb=1:100
        disp(['Iteration nb: ', num2str(iter_nb)])
        [cluster_assign{ii-1}{iter_nb}, ~] = kmeans(sqE,ii,'MaxIter',iter_nb);
    end
end

%% Get sihouette value
%calculate the silhoutte value for each cluster, take max
parfor cc=1:18
    disp(cc);
    for iter_nb=1:100
        %mean silhouette value for each iteration
        sil_vals{cc}{iter_nb} = silhouette(norm_cov_SCE,cluster_assign{cc}{iter_nb});
    end
end

%mean of silouette values
for cc=1:18
    mean_sil_vals{cc} = cellfun(@mean,sil_vals{cc},'UniformOutput',true);
end

figure
boxplot(cell2mat(mean_sil_vals')')

%generate cluster sort order
%numerical order
idx_order = size(cov_SCE,1);
%create sorted cluster
idx_cluster_order = [];
nb_clust = 10;
%make cluster order index vector
for ii=1:nb_clust
    idx_cluster_order = [idx_cluster_order;find(cluster_assign{nb_clust}{44} ==ii)];
end

%2 cluster reordered matrix
norm_cov_SCE_clust= norm_cov_SCE(idx_cluster_order, idx_cluster_order);

cov_SCE_clust =  cov_SCE(idx_cluster_order, idx_cluster_order);

figure
%imagesc(norm_cov_SCE_clust)
imagesc(cov_SCE_clust)
hold on
colormap('jet')
colorbar


%% Try sorting by assembly matrix by cluster numbers
%sort SCEs according to clusters

%sub-cluster order
cluster_sizes = histcounts(cluster_assign{10},0.5:1:11.5);
%sort in descending order
[~,cluster_sort_idx] = sort(cluster_sizes,'descend');

%preallocate
SCE_cluster_sort = [];

%for each cluster
for cc=cluster_sort_idx
    SCE_cluster_add{cc} = find(cluster_assign{10} == cc)
    SCE_cluster_sort = [SCE_cluster_sort; SCE_cluster_add{cc}];
end
figure;
imagesc(assembly_mat(:,SCE_cluster_sort))

%sort ROIs in each cluster in descending order
%sum events and sort is descending order
[~,ROI_sort] = sort(sum(assembly_mat(:,SCE_cluster_add{1, 3}),2),'descend');

assembly_mat(ROI_sort,SCE_cluster_add{1, 3})
%% Plot Euclidean distance matrix

figure;
imagesc(sqE)
hold on
colormap('jet')

%% Plot
figure;
imagesc(assembly_mat)
hold on
ylabel('Cell #');
xlabel('SCE #');
