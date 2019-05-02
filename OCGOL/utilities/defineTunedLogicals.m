function [tunedLogical] = defineTunedLogicals(animal_data)
%defines logical vectors for different classes of tuned neurons

%% Get logicals of neurons tuned using each metric SI or TS

for aa=1:size(animal_data)
    
    %using spatial information metric for each class
    tunedLogical(aa).si.Atuned = animal_data{aa}.Place_cell{1}.Spatial_Info.significant_ROI == 1;
    tunedLogical(aa).si.Btuned = animal_data{aa}.Place_cell{2}.Spatial_Info.significant_ROI == 1;
    
    %using tuning specificity metric for each class
    tunedLogical(aa).ts.Atuned = animal_data{aa}.Place_cell{1}.Tuning_Specificity.significant_ROI == 1;
    tunedLogical(aa).ts.Btuned = animal_data{aa}.Place_cell{2}.Tuning_Specificity.significant_ROI == 1;
    
    %si combininations
    tunedLogical(aa).si.AandB_tuned = tunedLogical(aa).si.Atuned & tunedLogical(aa).si.Btuned;
    tunedLogical(aa).si.AorB_tuned = tunedLogical(aa).si.Atuned | tunedLogical(aa).si.Btuned;
    
    tunedLogical(aa).si.onlyA_tuned = tunedLogical(aa).si.Atuned & ~tunedLogical(aa).si.Btuned;
    tunedLogical(aa).si.onlyB_tuned = ~tunedLogical(aa).si.Atuned & tunedLogical(aa).si.Btuned;
    
    tunedLogical(aa).si.neither = ~tunedLogical(aa).si.Atuned & ~tunedLogical(aa).si.Btuned;
    
    %ts combinations
    tunedLogical(aa).ts.AandB_tuned = tunedLogical(aa).ts.Atuned & tunedLogical(aa).ts.Btuned;
    tunedLogical(aa).ts.AorB_tuned = tunedLogical(aa).ts.Atuned | tunedLogical(aa).ts.Btuned;
    
    tunedLogical(aa).ts.onlyA_tuned = tunedLogical(aa).ts.Atuned & ~tunedLogical(aa).ts.Btuned;
    tunedLogical(aa).ts.onlyB_tuned = ~tunedLogical(aa).ts.Atuned & tunedLogical(aa).ts.Btuned;
    
    tunedLogical(aa).ts.neither = ~tunedLogical(aa).ts.Atuned & ~tunedLogical(aa).ts.Btuned;
    
end

end
