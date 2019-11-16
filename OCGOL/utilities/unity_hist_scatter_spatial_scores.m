function [outputArg1,outputArg2] = unity_hist_scatter_spatial_scores(si,options)

%split into 2 categories - below and above line for each class
%make positive vs. negative based on distance
si_cumul_above.Aonly = si.Aonly_cumul(si.Aonly_cumul(:,2) > si.Aonly_cumul(:,1),:);
si_cumul_below.Aonly = si.Aonly_cumul(si.Aonly_cumul(:,2) <= si.Aonly_cumul(:,1),:);

si_cumul_above.Bonly = si.Bonly_cumul(si.Bonly_cumul(:,2) > si.Bonly_cumul(:,1),:);
si_cumul_below.Bonly = si.Bonly_cumul(si.Bonly_cumul(:,2) <= si.Bonly_cumul(:,1),:);

si_cumul_above.AB = si.AB_cumul(si.AB_cumul(:,2) > si.AB_cumul(:,1),:);
si_cumul_below.AB = si.AB_cumul(si.AB_cumul(:,2) <= si.AB_cumul(:,1),:);

%si distance
%make negative
dist_unity.Aonly_above = (-1)*point_to_line_distance(si_cumul_above.Aonly,[0 0],[1 1]);
dist_unity.Bonly_above = (-1)*point_to_line_distance(si_cumul_above.Bonly,[0 0],[1 1]);
dist_unity.AB_above = (-1)*point_to_line_distance(si_cumul_above.AB,[0 0],[1 1]);

%make positive
dist_unity.Aonly_below = point_to_line_distance(si_cumul_below.Aonly,[0 0],[1 1]);
dist_unity.Bonly_below = point_to_line_distance(si_cumul_below.Bonly,[0 0],[1 1]);
dist_unity.AB_below = point_to_line_distance(si_cumul_below.AB,[0 0],[1 1]);

%merge
dist_unity.A_cumul = [dist_unity.Aonly_above; dist_unity.Aonly_below];
dist_unity.B_cumul = [dist_unity.Bonly_above; dist_unity.Bonly_below];
dist_unity.AB_cumul = [dist_unity.AB_above; dist_unity.AB_below];

%get max height of centerline based on bin height of A&B neurons
line_height = max(histcounts(dist_unity.AB_cumul,15,'Normalization','probability'));

%histogram
figure('Position',[2108 546 560 195])
hold on
set(gca,'TickLength',[0, 0])
axis off
xlim(options.xlims)
h1 = histogram(dist_unity.AB_cumul,15,'Normalization','probability','FaceColor', [139, 0, 139]./255)
h1.EdgeColor =  [1 1 1];

h2= histogram(dist_unity.A_cumul,15,'Normalization','probability','FaceColor',[65,105,225]./255)
h2.EdgeColor =  [1 1 1];

h3 =histogram(dist_unity.B_cumul,15,'Normalization','probability','FaceColor',[220,20,60]./255)
h3.EdgeColor =  [1 1 1];

%plot centerline
plot([0 0],[0 line_height],'k','LineWidth',2)


end

