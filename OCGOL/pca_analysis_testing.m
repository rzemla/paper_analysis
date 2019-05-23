

%small set for testing
%columns - ROIS
%rows - frame time
pca_input =traces(1:5000,input_neuron_idxs_sorted);

%all restricted time and traces
pca_input =traces;

%% Smooth calcium traces with gaussian with sigma = 5

%% Run basic PCA

%pca
[coef, score, latent,~,explained,mu] = pca(pca_input);

%extract 1st component (highest explained variance)
%principal component * loading
compNb = 1;
first_comp_recon = score(:,compNb)*coef(:,compNb)';

%plot as colormap
figure;
imagesc(first_comp_recon')
hold on;
title(['PCA - component: ', num2str(compNb)]);
colorbar;

