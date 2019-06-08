%% Detects SCEs associated with reccurent cell assemblies (Fig. 2 Malvache 2016)

%% Rough SCE selection

select_SCE_thres = find(seq_overlap > sce_threshold);
select_SCE_run_seq_min = SCE_run_seqs;

select_SCE = intersect(select_SCE_thres, select_SCE_run_seq_min);

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
    %for each interation
    for iter=1:100
        [cluster_assign{ii-1}{iter}, ~] = kmeans(sqE,ii,'MaxIter',iter);
    end
end

%calculate the silhoutte value for each cluster, take max
%% Calculate silhoutte of each member of each cluster

%calculate silhoutte values for each set of clusters
for cc=1:size(cluster_assign,2)
    for iter=1:100
        sil_value{cc}{iter} = silhouette(sqE, cluster_assign{cc}{iter});
    end
end

%average silhouette values for each cluster
for cc=1:size(cluster_assign,2)
    avg_sil_values(:,cc)  = mean(cell2mat(sil_value{cc}));
end

%plot box plots of silhoutte values for each cluster set
figure
hold on;
title('Average silhoutte values of each k-means cluster');
boxplot(avg_sil_values)

%max sil value
[~,I_sil] = max(max(avg_sil_values));

%get average of each set of silhoutte values) and index of max average (old

%[~,I_sil] = max(mean(cell2mat(sil_value)));

%% Sort by clister size by assembly matrix by cluster numbers
%sort SCEs according to clusters

%sub-cluster order
cluster_sizes = histcounts(cluster_assign{I_sil}{100},0.5:1:(I_sil+1+0.5));
%sort in descending order
[~,cluster_sort_idx] = sort(cluster_sizes,'descend');

%preallocate
SCE_cluster_sort = [];

%for each cluster
for cc=cluster_sort_idx
    SCE_cluster_add{cc} = find(cluster_assign{I_sil}{100} == cc);
    SCE_cluster_sort = [SCE_cluster_sort; SCE_cluster_add{cc}];
end
figure;
imagesc(assembly_mat(:,SCE_cluster_sort))

%% Covariance matrix of cluster-sorted SCE

cov_asm_clusSort = cov(assembly_mat(:,SCE_cluster_sort));
cov_asm_clusSort = corrcoef(assembly_mat(:,SCE_cluster_sort));
figure;
imagesc(cov_asm_clusSort);
hold on;
xlabel('SCE');
ylabel('SCE');
colormap('jet')
colorbar


%% Now 

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
