function [outputArg1,outputArg2] = distanceToRew_vs_performance_scatter(mean_Distance_TS_sig,frac_rew)
%plot distance to reward vs performance (fraction of licks in reward zone)
%scatterplot
%inputs:
%1:vector of mean distances to reward
%2:vector of fraction of licks in reward zone values

%make scatter plot
figure;
hold on
ylim([0 0.9]);
xlim([0.7854 2.25]);
ylabel('Fraction of licks in reward zone')
xticks([0.7854, 1.1781,1.5708, 1.9635])
xticklabels({'\pi/4','3\pi/8','\pi/2','5\pi/8'})
xlabel('Mean distance to reward');
scatter(mean_Distance_TS_sig,frac_rew,'b');

end

