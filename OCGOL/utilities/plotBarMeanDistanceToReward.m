function [outputArg1,outputArg2] = plotBarMeanDistanceToReward(mean_Distance_TS_sig)

%input:
%mean distance to reward vector for each session (0-6 days)

figure;
hold on
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
title('Distance of place fields to reward');
ylabel('Distance to reward');
xticks([1 2 3 4 5 6 7]);
xticklabels({'RF Day 0','GOL 1 Day 1','GOL 1 Day 2','GOL 1 Day 3',...
    'GOL 2 Day 4','GOL 2 Day 5','GOL 2 Day 6'})
xtickangle(45)
yticks([0.7854,1.5708])
ylim([0.7854,2])
yticklabels({'\pi/4','\pi/2'})
bar(mean_Distance_TS_sig,'b')

%plot chance level difference (pi/2)
refL = refline(0,1.5708);
refL.Color = [0.5 0.5 0.5];
refL.LineStyle = '--';


end

