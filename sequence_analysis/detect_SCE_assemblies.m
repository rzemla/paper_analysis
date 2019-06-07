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

%normalized covariance matrix
norm_cov_SCE = corrcoef(assembly_mat);
%cov_SCE = cov(assembly_mat);

%calculated squared euclidean distance between columns of norm cov matrix
sqE = pdist2(norm_cov_SCE, norm_cov_SCE,'squaredeuclidean');

%% Run k-mean clustering on columns of square Euclidean metric
%100 iterations
%run for 2-19 clusters
%idx - SCE cluster assignments
for ii=2:19
    [cluster_assign{ii-1}, ~] = kmeans(sqE,ii,'MaxIter',100);
end

%calculate the silhoutte value for each cluster, take max

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
