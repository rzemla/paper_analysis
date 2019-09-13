function [tunedLogical] = defineTunedLogicals(animal_data,options)
%defines logical vectors for different classes of tuned neurons

%% Get logicals of neurons tuned using each metric SI or TS


%for each animal
for aa=options.sessionSelect%1:size(animal_data,2)
    %only use correct A or B vs all A or B
    if options.allCorrect == 1
        
        %using spatial information metric for each class
        tunedLogical(aa).si.Atuned = animal_data{aa}.Place_cell{1}.Spatial_Info.significant_ROI == 1;
        tunedLogical(aa).si.Btuned = animal_data{aa}.Place_cell{2}.Spatial_Info.significant_ROI == 1;
        
        %using tuning specificity metric for each class
        tunedLogical(aa).ts.Atuned = animal_data{aa}.Place_cell{1}.Tuning_Specificity.significant_ROI == 1;
        tunedLogical(aa).ts.Btuned = animal_data{aa}.Place_cell{2}.Tuning_Specificity.significant_ROI == 1;
    else
        %using spatial information metric for each class
        tunedLogical(aa).si.Atuned = animal_data{aa}.Place_cell{4}.Spatial_Info.significant_ROI == 1;
        tunedLogical(aa).si.Btuned = animal_data{aa}.Place_cell{5}.Spatial_Info.significant_ROI == 1;
        
        %using tuning specificity metric for each class
        tunedLogical(aa).ts.Atuned = animal_data{aa}.Place_cell{4}.Tuning_Specificity.significant_ROI == 1;
        tunedLogical(aa).ts.Btuned = animal_data{aa}.Place_cell{5}.Tuning_Specificity.significant_ROI == 1;
    end
        
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

    %A and B tuned by either SI or TS criteria
%find A&B tuned by either criterion (input option below)
%unclassified neurons get moved into mixed category
tuned_A_si_ts = tunedLogical(1).si.Atuned  | tunedLogical(1).ts.Atuned;
tuned_B_si_ts = tunedLogical(1).si.Btuned  | tunedLogical(1).ts.Btuned;
tunedLogical(aa).tuned_AB_si_ts = tuned_A_si_ts & tuned_B_si_ts;
    
end

end

