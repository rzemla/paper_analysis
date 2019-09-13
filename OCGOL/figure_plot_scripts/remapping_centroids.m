function [outputArg1,outputArg2] = remapping_centroids(path_dir,options)

%% Load in data from each session directory
for ee=1:size(path_dir,2)
    load_data_path{ee} = fullfile(path_dir{ee},'cumul_analysis','place_field_centers_remap.mat');
    binCenter_data{ee} = load(string(load_data_path{ee}));
end

%load(fullfile(path_dir{rr},'cumul_analysis','place_field_centers_remap.mat'),'bin_center');

%% Simple scatter plot for global remapper far
figure
hold on
title('far')
xlabel('Center of Place field A')
ylabel('Center of Place field B')
for ee=1:size(path_dir,2)
    scatter(binCenter_data{ee}.bin_center.global_far(1,:),binCenter_data{ee}.bin_center.global_far(2,:),'m')
end
plot([0 100],[0 100],'k--')

% 
merge_far_binCenters = [];
figure
hold on
title('near + far ')
xlabel('Center of Place field A')
ylabel('Center of Place field B')
for ee=1:size(path_dir,2)
    scatter(binCenter_data{ee}.bin_center.global_near(1,:),binCenter_data{ee}.bin_center.global_near(2,:),'g')
    %scatter(binCenter_data{ee}.bin_center.global_far(1,:),binCenter_data{ee}.bin_center.global_far(2,:),'r')
%merge far centroids
merge_far_binCenters = [merge_far_binCenters,binCenter_data{ee}.bin_center.global_far,binCenter_data{ee}.bin_center.global_near];
end

plot([0 100],[0 100],'k--')

%tile histogram
figure;
hold on
nbins = 15;
xedges = 0:5:100;
xedges = [0,33,66, 100];
yedges = xedges;
xlim([0 100]);
ylim([0 100])
h = histogram2(merge_far_binCenters(1,:),merge_far_binCenters(2,:),xedges,yedges,'DisplayStyle','tile','ShowEmptyBins','on');
%h = histogram2(merge_far_binCenters(1,:),merge_far_binCenters(2,:),nbins,'DisplayStyle','tile','ShowEmptyBins','on');
colormap('hot')
colorbar

%find all y's below x values
less_than_unity = merge_far_binCenters(1,:)>merge_far_binCenters(2,:);
figure
hold on
histogram2(merge_far_binCenters(1,less_than_unity),merge_far_binCenters(2,less_than_unity),xedges,yedges,'DisplayStyle','tile','ShowEmptyBins','on')

%only low track PV correlated animals
figure
hold on
title('Low PV correlated animal - far global')
xlabel('Center of Place field A')
ylabel('Center of Place field B')
for ee=options.lowPVcorr
    scatter(binCenter_data{ee}.bin_center.global_far(1,:),binCenter_data{ee}.bin_center.global_far(2,:),'m')
end
plot([0 100],[0 100],'k--')

%plot partial remappers
figure
hold on
for ee=1:size(path_dir,2)
for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
    if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
        scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
    else
        scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
    end
end
end

%% Plot common scatter
figure; hold on; scatter(binCenter_data{1}.bin_center.common(1,:),binCenter_data{1}.bin_center.common(2,:) )
%combine common centers into one matrix
common_centers = [];
for ee=1:size(path_dir,2)
    common_centers = [common_centers, binCenter_data{ee}.bin_center.common];
end


figure
histogram(mean(common_centers,1),0:5:100)

%% 2 subplots - store common vs. partial center for each A or B find
storeIdx =0;
figure
subplot(1,2,1)
hold on
for ee=1:size(path_dir,2)
    for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
        if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
            storeIdx = storeIdx + 1;
            scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
            %store B mean common vs. B partial in matrix
%all animals
            partial.com_parB(storeIdx,:) = [mean(binCenter_data{ee}.bin_center.partial_com(:,rr)), binCenter_data{ee}.bin_center.partial_far(2,rr)];
        else
            %scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
        end
    end
end
storeIdx =0;
subplot(1,2,2)
hold on
for ee=1:size(path_dir,2)
for rr=1:size(binCenter_data{ee}.bin_center.partial_com,2)
    if isnan(binCenter_data{ee}.bin_center.partial_far(1,rr))
        %scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(2,rr),'r')
    else
storeIdx = storeIdx + 1;
        scatter(mean(binCenter_data{ee}.bin_center.partial_com(:,rr)),binCenter_data{ee}.bin_center.partial_far(1,rr),'b')
            %store A mean common vs. A partial in matrix
%all animals            
partial.com_parA(storeIdx,:) = [mean(binCenter_data{ee}.bin_center.partial_com(:,rr)), binCenter_data{ee}.bin_center.partial_far(1,rr)];
    end
end
end

%% Convert partial mappings to negative of difference and plot histogram
diff_par.A = (-1)*(partial.com_parA(:,1) -  partial.com_parA(:,2));
diff_par.B = (-1)*(partial.com_parB(:,1) -  partial.com_parB(:,2));

A_fwd = find(diff_par.A > 0);
A_rew = find(diff_par.A < 0);
[h,p] = kstest2((-1)*diff_par.A(A_fwd), diff_par.A(A_rew))

B_fwd = find(diff_par.B > 0);
B_rew = find(diff_par.B < 0);
[h,p] = kstest2((-1)*diff_par.B(B_fwd), diff_par.B(B_rew))


%test difference of number of neurons ahead of centroid
kstest2(diff_par.A,diff_par.A);

%Ks test of difference between the 2 distribtuions
[h,p] = kstest2(diff_par.A,diff_par.B);

%ks test on distribution of partial field centers
[h,p] =kstest2(partial.com_parA(:,2),partial.com_parB(:,2))

%% Bin as a function of center and get boxplots of distrubution
%get bin indices
%10 even spaced bins
partial_bin_idx.A = discretize(partial.com_parA(:,1),[0:10:100]);
partial_bin_idx.B = discretize(partial.com_parB(:,1),[0:10:100]);

%3 behaviorally spaced bins (1-35, 35-75, 75-100)
partial_bin_idx_3.A = discretize(partial.com_parA(:,1),[0,35,75,100]);
partial_bin_idx_3.B = discretize(partial.com_parB(:,1),[0,35,75,100]);

%
%for each bin (10 even)
for bb=1:10
    count.A{bb} = partial.com_parA(find(partial_bin_idx.A == bb),2)
    count.B{bb} = partial.com_parB(find(partial_bin_idx.B == bb),2)
end

%for each bin (3 behavior)
for bb=1:3
    count_3.A{bb} = partial.com_parA(find(partial_bin_idx_3.A == bb),2)
    count_3.B{bb} = partial.com_parB(find(partial_bin_idx_3.B == bb),2)
end
%boxplot 2 wrapper(plot distributions of diff sizes
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%merge both cell in alternate manner
merge_hist_counts(1:2:20) = count.A;
merge_hist_counts(2:2:21) = count.B;

merge_hist_counts_3(1:2:5) = count_3.A;
merge_hist_counts_3(2:2:6) = count_3.B;

figure;
%Distribution of partial indexes  
hold on
boxplot2(merge_hist_counts_3)

%use matlab boxplot to generate groups
grouping_box{1} = repmat(1,size(count_3.A{1},1),1)
grouping_box{2} = repmat(2,size(count_3.B{1},1),1)
grouping_box{3} = repmat(3,size(count_3.A{2},1),1)
grouping_box{4} = repmat(4,size(count_3.B{2},1),1)
grouping_box{5} = repmat(5,size(count_3.A{3},1),1)
grouping_box{6} = repmat(6,size(count_3.B{3},1),1)
%catx=repmat({'I' 'II','III'},1,100);
c=cell2mat(grouping_box')';
c(find(c==1 | c==3 | c==5)) = 1;
c(find(c==2 |c==4 | c==6)) = 2

g=gramm('x',cell2mat(grouping_box'),'y',cell2mat(merge_hist_counts_3'),'color',c);
g.set_color_options('map',[65,105,225; 220,20,60; 255,0,255; 0 0 0]/255);

%%%
% Plot raw data as points
g.geom_point();
g.stat_boxplot('width',1,'dodge',2);
g.geom_vline('xintercept',[2.5,4.5])
% Set appropriate names for legends
g.set_names('column','Origin','x','Zone','y','Place field center','color','Trial type');
%%%
% Set figure title
%g.set_title('Fuel economy of new cars between 1970 and 1982');
%%%
% Do the actual drawing
figure('Position',[100 100 800 400]);
g.draw();

% 
% y = iosr.statistics.tab2box(cell2mat(grouping_box'),cell2mat(count_3.A'));
% 
% figure;
% bp = iosr.statistics.boxPlot(1:3,y,'Notch',false);
% bp.lineColor{2} = 'b';
% 
% figure;
% iosr.statistics.boxPlot([count_3.A(1);count_3.B(1)])

figure
boxplot(count_3.A{1},ones(size(count_3.A{1},1),1))

figure;
%Distribution of partial indexes  
hold on
subplot(1,2,1)
boxplot2(count_3.A)
hold on
subplot(1,2,2)
boxplot2(count_3.B)

[h,p] = ranksum(count_3.A{1},count_3.B{1})
[h,p] = ranksum(count_3.A{2},count_3.B{2})
[h,p] = ranksum(count_3.A{3},count_3.B{3})

figure
boxplot2([[count_3.A(1), count_3.B(1)];[count_3.A(2), count_3.B(2)]])

%% Boxplot of distributions relative to mean of centroid in each zone

merge_new = [merge_hist_counts_3([1,3,5]);merge_hist_counts_3([2,4,6])]

xlab={'1','2','3'} 
col=[220,20,60, 255;
65,105,225, 255;
]; 
%0, 0, 255, 200]; 
col=col/255;

f=figure 
hold on
%plot zone separator lines
[f2,x,group,positions,labelpos] =  multiple_boxplot(merge_new',xlab,{'A','B'},col') 
%overlay boxplot to add median line 
z= boxplot(x,group, 'positions', positions);
lines = findobj(z, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines, 'LineWidth',2)
xticks([1.3750, 2.1250, 2.875])
xticklabels({'Zone I', 'Zone II', 'Zone III'})
ylabel('Place field center')
plot([1.7500 ,1.7500],[-10 110],'k--','LineWidth',1.5)
plot([2.500 ,2.500],[-10 110],'k--','LineWidth',1.5)
set(gca,'FontSize',16)

%% Distribution of partially remapping fields for A and B plot and KS test for significanee (Figure 2F)
figure('Position',[2200 317 508 433])
hold on
title('Distribution of partially remapping fields')
[fA,xA] = ecdf(partial.com_parA(:,2))
[fB,xB] = ecdf(partial.com_parB(:,2))
ap = stairs(xA,fA,'LineWidth',2,'Color',[65,105,225]/255);
bp = stairs(xB,fB,'LineWidth',2,'Color',[220,20,60]/255);
xlabel(gca,'Track positon')
ylabel('Cumulative probability')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)
legend([ap bp],{'A','B'},'Location','northwest')

[h,p] =kstest2(partial.com_parB(:,2), partial.com_parA(:,2))

%normalized histogram for each set of partial remapping neurons (inset for
%cdf
figure('Position',[1199 400 208 420])
subplot(2,1,1)
hold on;
ylim([0 0.2])
yticks([0 0.1 0.2])
ylabel('Norm. density')
set(gca,'FontSize',16)
histogram(partial.com_parA(:,2),[0:10:100],'Normalization','probability','FaceColor', [65,105,225]/255)

subplot(2,1,2)
hold on
ylim([0 0.2])
yticks([0 0.1 0.2])
ylabel('Norm. density')
xlabel('Track position')
set(gca,'FontSize',16)
histogram(partial.com_parB(:,2),[0:10:100],'Normalization','probability','FaceColor',[220,20,60]/255)


figure;
hold on
title('')
xlabel('Distance of field from mean of place field center')
xlim([-100 100])
histogram(diff_par.A,-100:10:100)
histogram(diff_par.B,-100:10:100)
plot([0 0],[0 35],'k--')

figure;
hold on
title('')
xlabel('')
xlim([0 100])
ylim([0 40])
histogram(partial.com_parA(:,2),0:5:100)
%histogram(partial.com_parB(:,2),0:10:100)
%histogram(diff_par.B,-100:10:100)
%plot([0 0],[0 35],'k--')

%% For partial, scatter plot of centroids - sort my mean of common place field center
%play with later
%{
[~,Isort_com] = sort(mean(bin_center.partial_com),'descend');
%sort columns by common mean center
bin_center.partial_com_sort = bin_center.partial_com(:,Isort_com);
bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

figure
%plot common
hold on
for rr=1:size(bin_center.partial_com,2)
    scatter(bin_center.partial_com_sort(1,rr),rr,'g')
    scatter(bin_center.partial_com_sort(2,rr),rr,'g')
end
%plot partial
for rr=1:size(bin_center.partial_com,2)
if ~isnan(bin_center.partial_far_sort(1,rr))
    scatter(bin_center.partial_far_sort(1,rr),rr,'b')
else
    scatter(bin_center.partial_far_sort(2,rr),rr,'r')
end
end

%sort global far by B trials
[~,Isort_glo_far_B] = sort(bin_center.global_far(2,:),'descend');
[~,Isort_glo_far_A] = sort(bin_center.global_far(1,:),'descend');
%sort columns by common mean center
bin_center.global_far_sortB = bin_center.global_far(:,Isort_glo_far_B);
bin_center.global_far_sortA = bin_center.global_far(:,Isort_glo_far_A);
%bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

%Plot global far
figure;
subplot(1,2,1)
hold on;
for rr=1:size(bin_center.global_far,2)
    scatter(bin_center.global_far_sortA(1,rr),rr,'b')
    scatter(bin_center.global_far_sortA(2,rr),rr,'r')
end
subplot(1,2,2)
hold on;
for rr=1:size(bin_center.global_far,2)
    scatter(bin_center.global_far_sortB(1,rr),rr,'b')
    scatter(bin_center.global_far_sortB(2,rr),rr,'r')
end

figure;
hold on
histogram(bin_center.global_far(1,:),0:10:100)

figure;
hold on
histogram(bin_center.global_far(2,:),0:10:100)

figure
hold on
title('Far global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.global_far_sortA(1,:),bin_center.global_far_sortA(2,:))

figure
hold on
title('Near global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.global_near(1,:),bin_center.global_near(2,:))
plot([0 100],[0 100],'k--')

figure
hold on
title('Common - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(bin_center.common(1,:),bin_center.common(2,:))
plot([0 100],[0 100],'k--')

%}
end

