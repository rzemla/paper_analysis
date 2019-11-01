function [outputArg1,outputArg2] = cumulative_performance_plot(short_term_learn,short_term_recall)

%function of sessions and not days

%% Extract peformance variables

perf_learn = short_term_learn.perf;
perf_recall = short_term_recall.perf;


%% Sessions for each animal

%number of animals for learn and recall
nb_learn = size(perf_learn,2); 
nb_recall = size(perf_recall,2); 

%learn
for aa=1:nb_learn
    learn_ses_nb(aa) = size(perf_learn{aa}.ses_perf,2);
end

%recall
for aa=1:nb_recall
    recall_ses_nb(aa) = size(perf_recall{aa}.ses_perf,2);
end

%% Make one 3-D matrix with fractional performance from all animals

%recall (sessions collected from all animals - same)
for aa = 1:size(perf_recall,2)
    perf_recall_comb(:,:,aa) = perf_recall{aa}.ses_perf;
end

%learning
%make blank matrix (rows - combined,A,B; col - consecutive ses, 3rd dim -
%animal number)
perf_learn_comb = nan(3, max(learn_ses_nb),nb_learn);

for aa = 1:size(perf_learn,2)
    perf_learn_comb(:,1:learn_ses_nb(aa),aa) = perf_learn{aa}.ses_perf;
end

%get means
mean_recall_perf = nanmean(perf_recall_comb,3);
mean_learn_perf = nanmean(perf_learn_comb,3);

%get stds
std_recall_perf = nanstd(perf_recall_comb,0,3);
std_learn_perf = nanstd(perf_learn_comb,0,3);

%get number of animals per session (learning)
nb_learn_per_ses = sum(~isnan(perf_learn_comb),3);
nb_recall_per_ses = sum(~isnan(perf_recall_comb),3);

%get sems
sem_recall_perf = std_recall_perf./sqrt(nb_recall_per_ses);
sem_learn_perf = std_learn_perf./sqrt(nb_learn_per_ses);

%% Plot
%recall - dash
%learning - solid
f = figure('Position',[2210 350 510 470]);
hold on
set(f,'color','w');
axis square
%combined - change color for each type
%subplot(1,3,tt)
tt=1
xlim([0.5 9.5])
ylim([0 1.1])
set(gca,'linewidth',2)
set(gca,'FontSize',21)
yticks(0:0.2:1)
ylabel('Fraction of correct trials')
xticks((1:9))
xticklabels({'5A5B','5A5B','3A3B','3A3B','Random','Random','Random','Random','Random'})
xtickangle(45)
%plot performance line (85%)
%plot([0.5 6.5],[0.85 0.85],'k--','LineWidth',2);
%plot([1 2 3 6 7 8 9],mean_recall_perf(1,:),'m-')

if tt ==1
    %all
    color_vec = [139, 0, 139]/255;
elseif tt ==2
    %A
    color_vec = [65,105,225]/255;
elseif tt==3
    %B
    color_vec = [ 220,20,60]/255;
end
%errorbar + mean  - learning
errorbar(1:9,mean_learn_perf(tt,:),sem_learn_perf(tt,:),'Color', color_vec, 'LineStyle', '-','LineWidth',1.5)
%errorbar + mean - recall
%errorbar([1 2 3 6 7 8 9],mean_recall_perf(tt,:),sem_recall_perf(tt,:),'Color', color_vec, 'LineStyle', '-','LineWidth',1.5)

%save performance figure
disp('Saving performance figure ')
export_fig(f ,fullfile('G:\Google_drive\task_selective_place_paper\input_figures_to_illustrator\Figure_4_figures',...
    'learning_performance.svg'),'-r600')

%plot([1:6],mean_learning_perf(1,:),'LineStyle','--','LineWidth',2,'Color', [139, 0, 139]/255)

%shaded plot parameters
% s = shadedErrorBar(1:6,squeeze(perf_learning_comb(1,:,:))',{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
% set(s.edge,'LineWidth',1,'LineStyle','-','Color',[[139, 0, 139]/255, 0.2]) %last value add transparency value
% s.mainLine.LineWidth = 2;
% s.mainLine.Color = [139, 0, 139]/255;
% s.patch.FaceColor = [139, 0, 139]/255;


end

