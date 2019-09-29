function [SCE] = SCE_participation(session_vars,SCE,PV_TC_corr,options)

%% Define parameters

%which sessions to analyze
sessionSelect = options.sessionSelect;
selectTrial = options.selectTrial;

%A vs.B TC scores for each session
for ss=sessionSelect
    TC_scores{ss} = PV_TC_corr.TC_same_day_diag{ss};
end

%% Total ROIs in each session
%row - neuron on/off
%column - each SCE
for ss=1:size(sessionSelect,2)
    ses_nbROI(ss) = size(session_vars{ss}.Place_cell{selectTrial(1)}.Spatial_tuning_curve,2);
end

%% Create SCE activity matrix for each sessions


%fill in each matrix
%for each session
for ss =sessionSelect
    %create blank SCE matrix
    SCE{ss}.sce_activity_matrix = zeros(ses_nbROI(ss),SCE{ss}.nbSCE);
    
    %fill cell participation for each SCE
    for cc=1:SCE{ss}.nbSCE
        SCE{ss}.sce_activity_matrix(SCE{ss}.SCE_unique_ROIs{cc},cc) = 1;
    end
end

%% Assign A vs. B TC score to each ROI in each SCE

%get TC score of A vs. B trials for each
for ss=sessionSelect
    for cc=1:SCE{ss}.nbSCE
        SCE{ss}.TC_score_SCE_sorted{cc} = TC_scores{ss}(SCE{ss}.SCE_unique_ROIs_sorted{cc});
        
    end
end

%mean TC score of neurons participating in each SCE
for ss=sessionSelect
    SCE{ss}.meanTC = cellfun(@nanmean,SCE{ss}.TC_score_SCE_sorted);
end

%boxplot wrapper for different size vectors
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%figure
%boxplot2(SCE{2}.TC_score_SCE_sorted)

%% Get SCE position of onset
%get normalized postion vector for each session
for ss=sessionSelect
    pos_norm{ss} = session_vars{ss}.Behavior.resampled.normalizedposition;
end

for ss=sessionSelect
    %frames onset of each SCE use as input to norm position
    SCE{ss}.SCE_pos_start = pos_norm{ss}(SCE{ss}.sync_idx(SCE{ss}.sync_range(:,1)));
end


%% Split SCEs in categories 
%compare side-by-side

for ss=sessionSelect
    %check if in learning session or recall session
    if selectTrial(1) == 1
        %for recall take in only correct and incorrect trials
        SCE{ss}.sce_activity.A = SCE{ss}.sce_activity_matrix(:,(SCE{ss}.sce_assign_log ==2));
        SCE{ss}.sce_activity.B = SCE{ss}.sce_activity_matrix(:,(SCE{ss}.sce_assign_log ==3));
    %for learning take both correct and incorrect A/B trials
    elseif selectTrial(1) == 4
        SCE{ss}.sce_activity.A = SCE{ss}.sce_activity_matrix(:,(SCE{ss}.sce_assign_log ==2 | SCE{ss}.sce_assign_log ==20));
        SCE{ss}.sce_activity.B = SCE{ss}.sce_activity_matrix(:,(SCE{ss}.sce_assign_log ==3 | SCE{ss}.sce_assign_log ==30));        
        
    end
end



%counts of how many SCEs each ROI is involved in
for ss=sessionSelect
    %first col - A; second col - B; last - all
    SCE{ss}.ROI_SCE_count(:,1) = sum(SCE{ss}.sce_activity.A,2);
    SCE{ss}.ROI_SCE_count(:,2) = sum(SCE{ss}.sce_activity.B,2);
    SCE{ss}.ROI_SCE_count(:,3) = sum(SCE{ss}.sce_activity_matrix,2);
end

% figure
% imagesc(cov(SCE{ss}.sce_activity.A))
% hold on
% colormap('jet')
% %caxis([0 1])
% colorbar

% figure
% imagesc([SCE{ss}.sce_activity.A, SCE{ss}.sce_activity.A])
% 
% corr(SCE{ss}.sce_activity.A,SCE{ss}.sce_activity.B)

% %find neurons that are repeatedly recruited by SCE
% roi_sce_part.A = sum(sce_activity{ss}.A,2);
% roi_sce_part.B = sum(sce_activity{ss}.B,2);
% 
% %combine in one mattrix
% roi_sce_split = [roi_sce_part.A,roi_sce_part.B ];


% 



end

