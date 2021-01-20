function [ROI_sig] = whichSignificant2(Place_cell)

    if(~iscell(Place_cell))
        Place_cell = {Place_cell};
    end

    minEvents = 3;

    pval_tuning_specificity = sum(bsxfun(@lt,Place_cell{1}.Tuning_Specificity.tuning_specificity,Place_cell{1}.Tuning_Specificity.Shuffle.tuning_specificity),1)./sum(~isnan(bsxfun(@lt,Place_cell{1}.Tuning_Specificity.tuning_specificity,Place_cell{1}.Tuning_Specificity.Shuffle.tuning_specificity)),1);
    null95 = quantile(Place_cell{1}.Tuning_Specificity.Shuffle.tuning_specificity,0.95);
%     sig_Tuning_Specificity = Place_cell{1}.Tuning_Specificity.tuning_specificity>null95;
%     

    
    
    allShuffles = cell2mat(reshape(Place_cell{1}.Spatial_Info.Shuffle.spatial_info,1,1,[]));
    shuffleMean = squeeze(mean(allShuffles,3));
    corrected_shuffles = bsxfun(@minus, allShuffles, shuffleMean);
    
    original_spatial_info = Place_cell{1}.Spatial_Info.Spatial_Info;
    corrected_spatial_info = original_spatial_info-shuffleMean;

    maxInfo_original_shuffle = squeeze(max(allShuffles,[],1));
    maxInfo_corrected_shuffle = squeeze(max(corrected_shuffles,[],1));

    maxInfo_original_info = squeeze(max(original_spatial_info,[],1));
    maxInfo_corrected_info = squeeze(max(corrected_spatial_info,[],1));
    
    pval_corrected_v_corrected = sum(bsxfun(@lt,maxInfo_corrected_info',maxInfo_corrected_shuffle),2)./sum(~isnan(maxInfo_corrected_shuffle),2);
    pval_original_v_original = sum(bsxfun(@lt,maxInfo_original_info',maxInfo_original_shuffle),2)./sum(~isnan(maxInfo_original_shuffle),2);
    pval_corrected_v_original = sum(bsxfun(@lt,maxInfo_corrected_info',maxInfo_original_shuffle),2)./sum(~isnan(maxInfo_original_shuffle),2);
    
    
    
    
    
    
    allShuffles2 = cell2mat(reshape(Place_cell{1}.Spatial_Info.Shuffle.spatial_info_Skaggs,1,1,[]));
    shuffleMean2 = squeeze(mean(allShuffles2,3));
    corrected_shuffles2 = bsxfun(@minus, allShuffles2, shuffleMean2);
    
    original_spatial_info2 = Place_cell{1}.Spatial_Info.Spatial_Info_Skaggs;
    corrected_spatial_info2 = original_spatial_info2-shuffleMean2;

    maxInfo_original_shuffle2 = squeeze(max(allShuffles2,[],1));
    maxInfo_corrected_shuffle2 = squeeze(max(corrected_shuffles2,[],1));

    maxInfo_original_info2 = squeeze(max(original_spatial_info2,[],1));
    maxInfo_corrected_info2 = squeeze(max(corrected_spatial_info2,[],1));
    
    pval_corrected_v_corrected2 = sum(bsxfun(@lt,maxInfo_corrected_info2',maxInfo_corrected_shuffle2),2)./sum(~isnan(maxInfo_corrected_shuffle2),2);
    pval_original_v_original2 = sum(bsxfun(@lt,maxInfo_original_info2',maxInfo_original_shuffle2),2)./sum(~isnan(maxInfo_original_shuffle2),2);
    pval_corrected_v_original2 = sum(bsxfun(@lt,maxInfo_corrected_info2',maxInfo_original_shuffle2),2)./sum(~isnan(maxInfo_original_shuffle2),2);
    
    
    
    nEvents = sum(Place_cell{1}.Spatial_Info.Run_onset_bin,1);
    
    ROI_sig.Tuning_Specificity = find((pval_tuning_specificity<0.05));
    
    ROI_sig.Info_Content.cvc = find(pval_corrected_v_corrected2'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Content.ovo = find(pval_original_v_original2'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Content.cvo = find(pval_corrected_v_original2'<0.05 & (nEvents>=minEvents));
    
    ROI_sig.Info_Content.cvc = find(pval_corrected_v_corrected2'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Content.ovo = find(pval_original_v_original2'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Content.cvo = find(pval_corrected_v_original2'<0.05 & (nEvents>=minEvents));
    
    ROI_sig.Info_Rate.cvc = find(pval_corrected_v_corrected'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Rate.ovo = find(pval_original_v_original'<0.05 & (nEvents>=minEvents));
    ROI_sig.Info_Rate.cvo = find(pval_corrected_v_original'<0.05 & (nEvents>=minEvents));
    
    ROI_sig.Overall.ic_cvc = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Content.cvc);
    ROI_sig.Overall.ic_ovo = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Content.ovo);
    ROI_sig.Overall.ic_cvo = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Content.cvo);
    
    ROI_sig.Overall.ir_cvc = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Rate.cvc);
    ROI_sig.Overall.ir_ovo = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Rate.ovo);
    ROI_sig.Overall.ir_cvo = intersect(ROI_sig.Tuning_Specificity,ROI_sig.Info_Rate.cvo);
    
end