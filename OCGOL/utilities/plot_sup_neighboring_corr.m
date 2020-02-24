function [outputArg1,outputArg2] = plot_sup_neighboring_corr(PV_mean_sem,TC_mean_sem)


%% Short term data only

%% Plot data ordered for iterative input - 1-6 index subplotting

input_data.A{1} = PV_mean_sem.neighbor.st_learn.A(:,1:6);
input_data.B{1} = PV_mean_sem.neighbor.st_learn.B(:,1:6);

input_data.A{2} = PV_mean_sem.neighbor.st_recall.A;
input_data.B{2} = PV_mean_sem.neighbor.st_recall.B;

input_data.A{3} = TC_mean_sem.neighbor.ts.st_learn.A;
input_data.B{3} = TC_mean_sem.neighbor.ts.st_learn.B;

input_data.A{4} = TC_mean_sem.neighbor.ts.st_recall.A;
input_data.B{4} = TC_mean_sem.neighbor.ts.st_recall.B;

input_data.A{5} = TC_mean_sem.neighbor.si.st_learn.A;
input_data.B{5} = TC_mean_sem.neighbor.si.st_learn.B;

input_data.A{6} = TC_mean_sem.neighbor.si.st_recall.A;
input_data.B{6} = TC_mean_sem.neighbor.si.st_recall.B;

title_labels{1} = 'PV correlation - Learn';
title_labels{2} = 'PV correlation - Recall';
title_labels{3} = 'TC correlation - TS - Learn';
title_labels{4} = 'TC correlation - TS - Recall';
title_labels{5} = 'TC correlation - SI - Learn';
title_labels{6} = 'TC correlation - SI - Recall';


%%
figure('Position',[2567 70 657 926])
%learn vs raw PV
for ii=1:6
    subplot(3,2,ii)
    hold on
    axis square
    ylim([0 1])
    xlabel('Day comparison')
    ylabel('Correlation score')
    xticks([1:8])
    xlim([0 9])
    xticklabels({'1 vs. 2','2 vs. 3','3 vs. 4','4 vs. 5','5 vs. 6','6 vs. 7','7 vs. 8','8 vs. 9'})
    xtickangle(45)
    title(title_labels{ii})
    %learn
    %if left panel
    if intersect(ii, 1:2:5)
        lA = plot_error_line(input_data.A{ii},'--',2,[65,105,225]/255);
        lB = plot_error_line(input_data.B{ii},'--',2,[220,20,60]/255);
    else %right panel
        lA = plot_error_line(input_data.A{ii},'-',2,[65,105,225]/255);
        lB = plot_error_line(input_data.B{ii},'-',2,[220,20,60]/255);
    end
    
    set(gca,'FontSize',12)
    set(gca,'Linewidth',1.5)
    
    legend([lA,lB],{'A','B'},'location','southwest')
end



end

