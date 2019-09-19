function [tunedLogical] = defineTunedLogicals(animal_data,options)
%defines logical vectors for different classes of tuned neurons

%% Get logicals of neurons tuned using each metric SI or TS


%for each animal
for ss=options.sessionSelect%1:size(animal_data,2)
    %only use correct A or B vs all A or B
    if options.allCorrect == 1
        
        %using spatial information metric for each class
        tunedLogical(ss).si.Atuned = animal_data{ss}.Place_cell{1}.Spatial_Info.significant_ROI == 1;
        tunedLogical(ss).si.Btuned = animal_data{ss}.Place_cell{2}.Spatial_Info.significant_ROI == 1;
        
        %using tuning specificity metric for each class
        tunedLogical(ss).ts.Atuned = animal_data{ss}.Place_cell{1}.Tuning_Specificity.significant_ROI == 1;
        tunedLogical(ss).ts.Btuned = animal_data{ss}.Place_cell{2}.Tuning_Specificity.significant_ROI == 1;
    else
        %using spatial information metric for each class
        tunedLogical(ss).si.Atuned = animal_data{ss}.Place_cell{4}.Spatial_Info.significant_ROI == 1;
        tunedLogical(ss).si.Btuned = animal_data{ss}.Place_cell{5}.Spatial_Info.significant_ROI == 1;
        
        %using tuning specificity metric for each class
        tunedLogical(ss).ts.Atuned = animal_data{ss}.Place_cell{4}.Tuning_Specificity.significant_ROI == 1;
        tunedLogical(ss).ts.Btuned = animal_data{ss}.Place_cell{5}.Tuning_Specificity.significant_ROI == 1;
    end
        
    %si combininations
    tunedLogical(ss).si.AandB_tuned = tunedLogical(ss).si.Atuned & tunedLogical(ss).si.Btuned;
    tunedLogical(ss).si.AorB_tuned = tunedLogical(ss).si.Atuned | tunedLogical(ss).si.Btuned;
    
    tunedLogical(ss).si.onlyA_tuned = tunedLogical(ss).si.Atuned & ~tunedLogical(ss).si.Btuned;
    tunedLogical(ss).si.onlyB_tuned = ~tunedLogical(ss).si.Atuned & tunedLogical(ss).si.Btuned;
    
    tunedLogical(ss).si.neither = ~tunedLogical(ss).si.Atuned & ~tunedLogical(ss).si.Btuned;
    
    %ts combinations
    tunedLogical(ss).ts.AandB_tuned = tunedLogical(ss).ts.Atuned & tunedLogical(ss).ts.Btuned;
    tunedLogical(ss).ts.AorB_tuned = tunedLogical(ss).ts.Atuned | tunedLogical(ss).ts.Btuned;
    
    tunedLogical(ss).ts.onlyA_tuned = tunedLogical(ss).ts.Atuned & ~tunedLogical(ss).ts.Btuned;
    tunedLogical(ss).ts.onlyB_tuned = ~tunedLogical(ss).ts.Atuned & tunedLogical(ss).ts.Btuned;
    
    tunedLogical(ss).ts.neither = ~tunedLogical(ss).ts.Atuned & ~tunedLogical(ss).ts.Btuned;

    %A and B tuned by either SI or TS criteria
%find A&B tuned by either criterion (input option below)
%unclassified neurons get moved into mixed category
tuned_A_si_ts{ss} = tunedLogical(ss).si.Atuned  | tunedLogical(ss).ts.Atuned;
tuned_B_si_ts{ss} = tunedLogical(ss).si.Btuned  | tunedLogical(ss).ts.Btuned;
tunedLogical(ss).tuned_AB_si_ts = tuned_A_si_ts{ss} & tuned_B_si_ts{ss};
    
end

end

