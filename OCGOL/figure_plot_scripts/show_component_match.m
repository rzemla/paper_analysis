function [outputArg1,outputArg2] = show_component_match(CNMF_learn,reg_learn)

%% Define animal and sessions to use
%which animal
learn_animal = 2;
% which two sessions to display
ses_comp = [1 3];

%% Select the shared filtered matches between sessions
%remove soma parsed components from Coor_kp - coordinate outlines from CNMF
%ses_1
soma_keep{1} = CNMF_learn.removeROI_learn{learn_animal}{ses_comp(1)}.compSelect;
%ses 2
soma_keep{2} = CNMF_learn.removeROI_learn{learn_animal}{ses_comp(2)}.compSelect;
%coordinate outlines from both sessions
coor_keep{1} = CNMF_learn.CNMF_vars_learn{learn_animal}{ses_comp(1)}.Coor_kp(soma_keep{1});
coor_keep{2} = CNMF_learn.CNMF_vars_learn{learn_animal}{ses_comp(2)}.Coor_kp(soma_keep{2});

%extract filtered matching matrix for the two days
extract_match_idx = find(sum(~isnan(reg_learn{learn_animal}.registered.multi.assigned_filtered(:,ses_comp)),2) ==2);
%select session match matrix
vis_match_matrix = reg_learn{learn_animal}.registered.multi.assigned_filtered(extract_match_idx,ses_comp);

%% Extract all matches selected by the program between sessions
extract_match_idx_all = find(sum(~isnan(reg_learn{learn_animal}.registered.multi.assigned(:,ses_comp)),2) ==2);
vis_match_matrix_all = reg_learn{learn_animal}.registered.multi.assigned(extract_match_idx_all,ses_comp);

%find overlapping componenets with filted out matches and remove
[~,remove_overlap_idx,~] = intersect(vis_match_matrix_all(:,1),vis_match_matrix(:,1));

%componenet match with exclusions
vis_match_matrix_excl = vis_match_matrix_all;
vis_match_matrix_excl(remove_overlap_idx,:) = [];

%% Plot match filtered componenets between the selected sessions - included
%plot BW outline of the component
f= figure('Position',[2100 150 1460 540]);
%set figure background color to white
set(f,'color','w');

subaxis(1,3,1, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
%subplot(1,3,1)
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
title('Day 1')
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix(:,1)'
    %plot componenet outline
    
    plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1);
    %f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1,'EdgeColor','none');
    %alpha(f1,0.5)
end

%subplot(1,3,2)
subaxis(1,3,2,'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(2)}.template);
hold on
title('Day 3')
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)
%plot all selected ROIs as green
for ROI = vis_match_matrix(:,2)'
    %plot componenet outline
    plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m', 'LineWidth',1);
    %    f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    %alpha(f2,0.5)
end

%combined
subaxis(1,3,3, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);

imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
title('Merge')
axes(gca);
axis square
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix(:,1)'
    %plot componenet outline
    
    plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1);
     %f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    %alpha(f1,0.5)
end

for ROI = vis_match_matrix(:,2)'
    %plot componenet outline
    plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m', 'LineWidth',1);
%f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m','EdgeColor','none')
 %alpha(f2,0.5)
end

%save figure 4a component match example learning
%disp('Saving match ROIs STC ')
%export_fig(f ,fullfile('G:\Google_drive\task_selective_place_paper\input_figures_to_illustrator\Figure_4_figures',...
 %   'componenent_matching.png'),'-r300')


%% Plot match filtered componenets between the selected sessions - EXCLUDED
%plot BW outline of the component
f= figure('Position',[2100 150 1460 540]);
%set figure background color to white
set(f,'color','w');

subaxis(1,3,1, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
%subplot(1,3,1)
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
axes(gca);
axis square
title('Day 1')
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix_excl(:,1)'
    %plot componenet outline
    
    plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1);
    %f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c', 'LineWidth',1,'EdgeColor','none');
    %alpha(f1,0.5)
end

%subplot(1,3,2)
subaxis(1,3,2,'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(2)}.template);
hold on
axes(gca);
axis square
title('Day 3')
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)
%plot all selected ROIs as green
for ROI = vis_match_matrix_excl(:,2)'
    %plot componenet outline
    plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m', 'LineWidth',1);
    %    f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    %alpha(f2,0.5)
end

%combined
subaxis(1,3,3, 'SpacingHorizontal', 0.0015,...
    'SpacingVertical',0.001,...
    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);

imagesc(CNMF_learn.templates_learn{learn_animal}{ses_comp(1)}.template);
hold on
axes(gca);
axis square
title('Merge')
xticks(gca,[])
yticks(gca,[])
grayMap = brighten(gray,0.2);
colormap(gca,grayMap)

%plot all selected ROIs as green
for ROI = vis_match_matrix_excl(:,1)'
    %plot componenet outline
    
    plot(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'c-', 'LineWidth',1);
    % f1= fill(coor_keep{1}{ROI}(1,:),coor_keep{1}{ROI}(2,:),'y', 'LineWidth',1,'EdgeColor','none');
    %alpha(f1,0.5)
end

for ROI = vis_match_matrix_excl(:,2)'
    %plot componenet outline
    plot(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m-', 'LineWidth',1);
    %f2= fill(coor_keep{2}{ROI}(1,:),coor_keep{2}{ROI}(2,:),'m','EdgeColor','none')
    %alpha(f2,0.5)
end

%save figure 4a component match example learning
disp('Saving match excluded ROIs STC ')
export_fig(f ,fullfile('G:\Google_drive\task_selective_place_paper\input_figures_to_illustrator\Figure_4_figures',...
    'componenent_matching_ex.png'),'-r300')

end

