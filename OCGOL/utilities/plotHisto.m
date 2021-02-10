function [outputArg1,outputArg2] = plotHisto(lickCombined,edges,barColor)

Aarea = zeros(1,50);
%Aarea(28:30) =1; 
Aarea(35:38) =1; 

Barea = zeros(1,50);
%Barea(12:14) =1; 
Barea(15:18) =1; 

%zone colors
red_z = [220,20,60]./255;
blue_z = [65,105,225]./255;

bar(Aarea,1,'FaceColor',blue_z,'FaceAlpha', 0.2,'EdgeColor','none')
bar(Barea,1,'FaceColor',red_z,'FaceAlpha', 0.2,'EdgeColor','none')

[N,edges,bin] = histcounts(lickCombined,edges,'Normalization','probability');
bar(N, 1, 'EdgeColor', 'none','FaceColor', barColor);
xticks([1 25 50])
xticklabels({'0.0','0.5','1.0'});
xlabel('Normalized position');

ylim([0 0.5])
yticks([0 0.5]);
yticklabels({'0.0','0.5'});
ylabel('Fraction of licks')

set(gca,'FontSize',12)
set(gca,'LineWidth',2)

 %rectangle('Position', [12 0 2 0.7 ], 'FaceColor','g', 'EdgeColor','none');
% rectangle('Position', [28 0 2 0.7 ], 'FaceColor','g', 'EdgeColor','none');

end

