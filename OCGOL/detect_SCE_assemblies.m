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

%calculated squared euclidean distance between columns of norm cov matrix
sqE = pdist2(norm_cov_SCE, norm_cov_SCE,'squaredeuclidean');

%% Plot
figure;
imagesc(assembly_mat)
hold on
ylabel('Cell #');
xlabel('SCE #');
