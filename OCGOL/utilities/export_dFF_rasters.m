function [rate_mean_dFF] = export_dFF_rasters(session_vars,remapping_corr_idx)

%remapping neuron idx
remap_neuron_idx = remapping_corr_idx.final.rate_remap_all;

for ii=1:size(remap_neuron_idx,2)
    %correct A trials
    rate_mean_dFF(ii,:) = [nanmean(session_vars{1}.Place_cell{1}.dF_lap_map_ROI{remap_neuron_idx(ii)},1),...
         nanmean(session_vars{1}.Place_cell{2}.dF_lap_map_ROI{remap_neuron_idx(ii)},1)];
    
    %correct B trials
   
end

%% 
figure
imagesc(rate_mean_dFF)
hold on
colormap('jet')
colorbar
caxis([0 2])
end

