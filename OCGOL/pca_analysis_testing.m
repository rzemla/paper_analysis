%% RUN sequence detection assembly


%% Input data

%small set for testing
%columns - ROIS
%rows - frame time
%pca_input =traces(st_idx:end_idx,input_neuron_idxs_sorted);

%all active restricted neurons
pca_input =traces;

%A traces

pca_input = session_vars{1, 1}.Imaging_split{1, 4}.trace_restricted;

%make raw copy
pca_input_raw = pca_input;
%all restricted time and traces
%pca_input =traces;

%% Smooth calcium traces with gaussian with sigma = 5

%5 sigma gaussian kernel
%5seconds -  150 frames
options.sigma_filter =  150;
gaussFilter = define_Gaussian_kernel(options);

for rr=1:size(pca_input,2)
    pca_input(:,rr) =conv(pca_input(:,rr),gaussFilter, 'same');
end


%% Run basic PCA

%calculate mean of each ROI - get from PCA as mu

%pca
[coef, score, latent,~,explained,mu] = pca(pca_input);

%extract 1st component (highest explained variance)
%principal component * loading
figure;
for ii=1:5
compNb = ii;
first_comp_recon = score(:,compNb)*coef(:,compNb)';
% add mean for each ROI
%first_comp_recon = bsxfun(@plus,first_comp_recon,mu);
first_comp_recon = first_comp_recon + mu;

%plot as colormap - display first 5 comps
subplot(1,5,ii)
imagesc(first_comp_recon')
hold on;
title(['PCA - component: ', num2str(compNb)]);
colorbar;
caxis([0 1])
end

%% Derivative of principal componenet for use as template

%derivate of first PCA
diff_comp = diff(score(:,2));

%% Display PCA

%% Correlate to each neuron (smoothed? or raw) - try both
%pad with 
%same result if raw or smoothed
corr_values = corr(pca_input, [0;diff_comp]);

%plot histogram
figure;
hold on
histogram(corr_values, 10)


%% Otsu's method for correlation threshold calculation (done on entire population of active cells)

[corr_counts,~] = histcounts(corr_values);
%calculate Otsu' threshold
T = otsuthresh(corr_counts);

find(corr_values >  T)

%% For offset PCA - similar to PCA test - 
%construct covariance matrix from original dF/F and 1 timeframes shifted
%matrix (circshift) - 2 separate matrices


%use cov to construct covariance matrix
%use pcacov to run pca on cov matrix

%try this approach first to get same output as with regular PCA
%need to stack the two matrixs together and select correct set do PCA; if
%loaded separate will run as two long vectors

%original matrix
orig_mat = pca_input;
%remove last row
orig_mat = orig_mat(1:end-1,:);

%shifted matrix
shift_mat = circshift(pca_input,1);
%remove first row
shift_mat = shift_mat(2:end,:);

%get mean of each ROI (variable)
data_mean = repmat(mean(pca_input,1),size(pca_input,1)-1,1);


%center the inputs as well
cov_pca =  cov([orig_mat-data_mean,shift_mat-data_mean]);
%cov_pca1 =  cov([pca_input,pca_input]);
%extract correct sub-matrix that corresponds to the cross-covariance
%cov_cross = cov_pca1(1:16,17:32);
cov_cross = cov_pca(569:end,1:568);

%cov_pca2 = cov([pca_input,circshift(pca_input,3)]);

%pca on covariance matrix (same explained variance as reg pca)
[coef_pcacov, latened_pcacov, explained_pca_cov] = pcacov(cov_cross);

%scores - representation of X (data) in principal component space
%need to generate this to reconstruct the data

scores_cov = orig_mat*coef_pcacov;
%reconstruct using chosen component
compNb = 2;
first_comp_recon_pca_cov = scores_cov(:,compNb)*coef_pcacov(:,compNb)'  + data_mean;

%plot as colormap
figure;
imagesc(first_comp_recon_pca_cov')
hold on;
title(['PCA - component: ', num2str(compNb)]);
colorbar;


%% Derivative of principal componenet for use as template

%derivate of first PCA
diff_comp_cov = diff(scores_cov(:,5));

%% Display PCA

%% Correlate to each neuron (smoothed? or raw) - try both
%pad with 
%same result if raw or smoothed
corr_values_cov = corr(orig_mat, [0;diff_comp_cov]);

%plot histogram
figure;
hold on
histogram(corr_values_cov)


%% Otsu's method for correlation threshold calculation (done on entire population of active cells)

[corr_counts_cov,~] = histcounts(corr_values_cov,100);
%calculate Otsu' threshold
T = otsuthresh(corr_counts_cov);

find(corr_values_cov >  T)

%% OLD BELOW

%% PCA using covariance matrix (no double matrix for offset test) - as good as above
data_mean = repmat(mean(pca_input,1),size(pca_input,1),1);
%center the inputs as well
%cov_pca1 =  cov(pca_input-data_mean);
cov_pca1 =  cov(pca_input);

%pca on covariance matrix (same explained variance as reg pca)
[coef_pcacov, latened_pcacov, explained_pca_cov] = pcacov(cov_pca1);

%scores - representation of X (data) in principal component space
%need to generate this to reconstruct the data

scores_cov = pca_input*coef_pcacov;
%reconstruct using chosen component
compNb = 1;
first_comp_recon_pca_cov = scores_cov(:,compNb)*coef_pcacov(:,compNb)'  + data_mean;

%plot as colormap
figure;
imagesc(first_comp_recon_pca_cov')
hold on;
title(['PCA - component: ', num2str(compNb)]);
colorbar;

diff_mat = first_comp_recon_pca_cov' - first_comp_recon';


