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
%mean of common field first and then partial field
%mean of common bin goes here
%group into 10 bins
partial_bin_idx.A = discretize(partial.com_parA(:,1),[0:10:100]);
partial_bin_idx.B = discretize(partial.com_parB(:,1),[0:10:100]);

%3 behaviorally spaced bins (1-35, 35-75, 75-100)
%split mean of common bin here by zones - %group into 3 bins
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

% g=gramm('x',cell2mat(grouping_box'),'y',cell2mat(merge_hist_counts_3'),'color',c);
% g.set_color_options('map',[65,105,225; 220,20,60; 255,0,255; 0 0 0]/255);

%%%
% Plot raw data as points
% g.geom_point();
% g.stat_boxplot('width',1,'dodge',2);
% g.geom_vline('xintercept',[2.5,4.5])
% % Set appropriate names for legends
% g.set_names('column','Origin','x','Zone','y','Place field center','color','Trial type');
% %%%
% % Set figure title
% %g.set_title('Fuel economy of new cars between 1970 and 1982');
% %%%
% % Do the actual drawing
% figure('Position',[100 100 800 400]);
% g.draw();


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

%% Distance relative to PF center position
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
histogram(partial.com_parB(:,2),0:5:100)
%histogram(diff_par.B,-100:10:100)
%plot([0 0],[0 35],'k--')

%%
%% For partial, scatter plot of centroids - sort my mean of common place field center

figure
hold on
for ss=1:11
    
    [~,Isort_com] = sort(mean(binCenter_data{ss}.bin_center.partial_com),'descend');
    %sort columns by common mean center
    bin_center.partial_com_sort = binCenter_data{ss}.bin_center.partial_com(:,Isort_com);
    bin_center.partial_far_sort = binCenter_data{ss}.bin_center.partial_far(:,Isort_com);
    
    
    %plot common
    hold on
    for rr=1:size(binCenter_data{ss}.bin_center.partial_com,2)
        scatter(bin_center.partial_com_sort(1,rr),rr,'g')
        scatter(bin_center.partial_com_sort(2,rr),rr,'g')
    end
    
    %plot partial
    for rr=1:size(binCenter_data{ss}.bin_center.partial_com,2)
        if ~isnan(bin_center.partial_far_sort(1,rr))
            scatter(bin_center.partial_far_sort(1,rr),rr,'b')
        else
            scatter(bin_center.partial_far_sort(2,rr),rr,'r')
        end
    end
end

figure;
for ss=1:11
%sort global far by B trials
[~,Isort_glo_far_B] = sort(binCenter_data{1}.bin_center.global_far(2,:),'descend');
[~,Isort_glo_far_A] = sort(binCenter_data{1}.bin_center.global_far(1,:),'descend');
%sort columns by common mean center
bin_center.global_far_sortB = binCenter_data{1}.bin_center.global_far(:,Isort_glo_far_B);
bin_center.global_far_sortA = binCenter_data{1}.bin_center.global_far(:,Isort_glo_far_A);
%bin_center.partial_far_sort =bin_center.partial_far(:,Isort_com);

%Plot global far

subplot(1,2,1)
hold on;
for rr=1:size(binCenter_data{1}.bin_center.global_far,2)
    scatter(bin_center.global_far_sortA(1,rr),rr,'b')
    scatter(bin_center.global_far_sortA(2,rr),rr,'r')
end
subplot(1,2,2)
hold on;
for rr=1:size(binCenter_data{1}.bin_center.global_far,2)
    scatter(bin_center.global_far_sortB(1,rr),rr,'b')
    scatter(bin_center.global_far_sortB(2,rr),rr,'r')
end
end

%% Make one common, global far, global near partial matrix
global_far.A = [];
global_far.B = [];
global_near.A = [];
global_near.B = [];
common = [];
partial_common = [];
rate = [];

for ee=1:11
    global_far.A = [global_far.A, binCenter_data{ee}.bin_center.global_far(1,:)];
    global_far.B = [global_far.B, binCenter_data{ee}.bin_center.global_far(2,:)];
    common = [common, binCenter_data{ee}.bin_center.common];
    partial_common = [partial_common,mean(binCenter_data{ee}.bin_center.partial_com,1) ];
    rate = [rate, binCenter_data{ee}.bin_center.rate];
    global_near.A = [global_near.A, binCenter_data{ee}.bin_center.global_near(1,:)]
    global_near.B = [global_near.B, binCenter_data{ee}.bin_center.global_near(2,:)]
end

%split by animal into separate cell
for ee=1:11
    global_near_each{ee} = binCenter_data{ee}.bin_center.global_near;
end


%combine global fars
global_far_comb = [global_far.A; global_far.B];

bin_space = 10;
ylim_def = 0.3;

figure
hold on
plot(global_far_comb)

figure
hold on
%scatter(global_near.A,global_near.B)
%unity line
axis square
plot([0 100], [0 100],'k-');
less_than_unity = global_near.A >= global_near.B;
greater_than_unity = global_near.A < global_near.B;
scatter(global_near.A(less_than_unity),global_near.B(less_than_unity),14,'MarkerFaceColor','b','MarkerEdgeColor','b') 
scatter(global_near.A(greater_than_unity),global_near.B(greater_than_unity),14,'MarkerFaceColor','r','MarkerEdgeColor','r')

%fit line
x = global_near.A';                                        % Create Data
y = global_near.B';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

[B,BINT,R,RINT,STATS] = regress( global_near.B',[global_near.A',ones(size(global_near.B,2),1)] )

plot(Xnew, YCI, '--k')
plot(x, Ypred,'-r')
hold off
%grid

%% Binomial test for shift (edges excluded)

%ecdf this
zone1_idx = find(global_near.A > 5 & global_near.A <=35);
zone2_idx = find(global_near.A > 35 & global_near.A <= 75);
zone3_idx = find(global_near.A > 75 & global_near.A <= 95);

%field difference B-A
field_pos_diff = global_near.B - global_near.A;

zone1_shifts = field_pos_diff(zone1_idx);
zone2_shifts = field_pos_diff(zone2_idx);
zone3_shifts = field_pos_diff(zone3_idx);

ratio_pos_neg.I = length(find(zone1_shifts > 0))/length(zone1_idx)
ratio_pos_neg.II = length(find(zone2_shifts > 0))/length(zone2_idx)
ratio_pos_neg.III = length(find(zone3_shifts > 0))/length(zone3_idx)

pOut = myBinomTest(length(find(zone1_shifts > 0)),length(zone1_idx),0.5)
pOut = myBinomTest(length(find(zone2_shifts > 0)),length(zone2_idx),0.5)
pOut = myBinomTest(length(find(zone3_shifts > 0)),length(zone3_idx),0.5)


%% Plot a stacked columns
figure
hold on
ba = bar([ratio_pos_neg.I, 1-ratio_pos_neg.I; ratio_pos_neg.II, 1-ratio_pos_neg.II; ratio_pos_neg.III, 1-ratio_pos_neg.III],...
        'stacked','FaceColor','flat')
xticks([1 2 3])
xticklabels({'Early','Mid','Late'})
%plot 50% chance line
plot([-0.2 4.5],[0.5 0.5],'k--','LineWidth',1)

ylim([0 1.2])
xlim([-0.2 4.2])
ylabel('Fraction of neurons');
yticks([0 0.2 0.4 0.6 0.8 1])
%set bar colors
ba(1).CData = [0,100,0]/255;
ba(2).CData = [1,1,1]*0.8;
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%% 6 zones - not more informative than 3 zones
if 0
    %ecdf this - break into 6 zones
    zone1_idx = find(global_near.A > 5 & global_near.A <=20);
    zone2_idx = find(global_near.A > 20 & global_near.A <=35);
    zone3_idx = find(global_near.A > 35 & global_near.A <= 60);
    zone4_idx = find(global_near.A > 60 & global_near.A <= 75);
    zone5_idx = find(global_near.A > 75 & global_near.A <= 85);
    zone6_idx = find(global_near.A > 85 & global_near.A <= 95);
    
    %field difference B-A
    field_pos_diff = global_near.B - global_near.A;
    
    zone1_shifts = field_pos_diff(zone1_idx);
    zone2_shifts = field_pos_diff(zone2_idx);
    zone3_shifts = field_pos_diff(zone3_idx);
    zone4_shifts = field_pos_diff(zone4_idx);
    zone5_shifts = field_pos_diff(zone5_idx);
    zone6_shifts = field_pos_diff(zone6_idx);
    
    ratio_pos_neg.I = length(find(zone1_shifts > 0))/length(zone1_idx)
    ratio_pos_neg.II = length(find(zone2_shifts > 0))/length(zone2_idx)
    ratio_pos_neg.III = length(find(zone3_shifts > 0))/length(zone3_idx)
    ratio_pos_neg.IV = length(find(zone4_shifts > 0))/length(zone4_idx)
    ratio_pos_neg.V = length(find(zone5_shifts > 0))/length(zone5_idx)
    ratio_pos_neg.VI = length(find(zone6_shifts > 0))/length(zone6_idx)
    
    pOut = myBinomTest(length(find(zone1_shifts > 0)),length(zone1_idx),0.5)
    pOut = myBinomTest(length(find(zone2_shifts > 0)),length(zone2_idx),0.5)
    pOut = myBinomTest(length(find(zone3_shifts > 0)),length(zone3_idx),0.5)
    pOut = myBinomTest(length(find(zone4_shifts > 0)),length(zone4_idx),0.5)
    pOut = myBinomTest(length(find(zone5_shifts > 0)),length(zone5_idx),0.5)
    pOut = myBinomTest(length(find(zone6_shifts > 0)),length(zone6_idx),0.5)
end

%% Classify global remappers into respective zones
%I-II
%II-III
%I-III
%ecdf this
zone_global_idx(1).A = find(global_far.A > 0 & global_far.A <= 35);
zone_global_idx(2).A = find(global_far.A > 35 & global_far.A <= 75);
zone_global_idx(3).A = find(global_far.A > 75 & global_far.A <= 100);

zone_global_idx(1).B = find(global_far.B > 0 & global_far.B <= 35);
zone_global_idx(2).B = find(global_far.B > 35 & global_far.B <= 75);
zone_global_idx(3).B = find(global_far.B > 75 & global_far.B <= 100);

%merge into one matrix and assign zones
merge_zone_matrix = zeros(2,size(global_far.A,2));

merge_zone_matrix(1,zone_global_idx(1).A) = 1;
merge_zone_matrix(1,zone_global_idx(2).A) = 2;
merge_zone_matrix(1,zone_global_idx(3).A) = 3;

merge_zone_matrix(2,zone_global_idx(1).B) = 1;
merge_zone_matrix(2,zone_global_idx(2).B) = 2;
merge_zone_matrix(2,zone_global_idx(3).B) = 3;

%count each type of transition
%I-II, II-I
one_two = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 2);
two_one = (merge_zone_matrix(2,:) == 1) & (merge_zone_matrix(1,:) == 2);
one_two_all = one_two | two_one;

%1-III, III-I
one_three = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 3);
three_one = (merge_zone_matrix(2,:) == 1) & (merge_zone_matrix(1,:) == 3);
one_three_all = one_three | three_one;

%II-III,III-II
two_three = (merge_zone_matrix(1,:) == 2) & (merge_zone_matrix(2,:) == 3);
three_two = (merge_zone_matrix(2,:) == 2) & (merge_zone_matrix(1,:) == 3);
two_three_all = two_three | three_two;

%self I-I, II-II, III-III
self_one = (merge_zone_matrix(1,:) == 1) & (merge_zone_matrix(2,:) == 1);
self_two = (merge_zone_matrix(1,:) == 2) & (merge_zone_matrix(2,:) == 2);
self_three = (merge_zone_matrix(1,:) == 3) & (merge_zone_matrix(2,:) == 3);

self_all = ((self_one | self_two) | self_three);

zone_count = [length(find(self_all == 1)), length(find(one_two_all == 1)), length(find(two_three_all == 1)),length(find(one_three_all == 1))]

%set up count for each transition (I-II, II-III, I-III)
transition_stacked = [length(find(one_two ==1)),length(find(two_one ==1));...
                    length(find(two_three ==1)),length(find(three_two ==1));...
                    length(find(one_three ==1)),length(find(three_one ==1))];

group_count = sum(transition_stacked,2);


figure;
hold on
ba2 = bar(transition_stacked./sum(transition_stacked,2),'stacked','FaceColor','flat');
xticks([1 2 3])
xticklabels({'Early <-> Mid','Mid <-> Late','Early <-> Late'})
ylim([0 1.2])
xlim([-0.2 4.2])
ylabel('Fraction of neurons');
yticks([0 0.2 0.4 0.6 0.8 1])
xtickangle(45);
%set bar colors
ba2(1).CData = [65,105,225]/255;
ba2(2).CData = [220,20,60]/255;
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%binomial test for transitions
pOut = myBinomTest(transition_stacked(1,1),group_count(1),0.5)
pOut = myBinomTest(transition_stacked(2,1),group_count(2),0.5)
pOut = myBinomTest(transition_stacked(3,1),group_count(3),0.5)

%% Common distribution 

%account for edge when calculating the mean of the common centers 
%find neurons with difference greater than 
split_common_idxs = find(diff(common,1) > 20 | diff(common,1) < -20);

%find min idx
%add 100 to min
%take mean
%subtract 100
for ii=1:size(split_common_idxs,2)
    [val_temp, idx_temp] = min(common(:,split_common_idxs(ii)));
    if idx_temp == 2
        upd_common_val(ii) = ceil(((val_temp+100) + common(1,split_common_idxs(ii)))/2 - 100);
    elseif idx_temp == 1
        upd_common_val(ii) = ceil(((val_temp+100) + common(2,split_common_idxs(ii)))/2 - 100);
    end
end
%change 0 outputs to 100 bin
upd_common_val((upd_common_val == 0)) = 100;

%take the mean and update values
mean_common = mean(common,1);
mean_common(split_common_idxs) = upd_common_val;


figure
hold on
ylabel('Normalized density')
xlabel('Track position')
histogram(mean_common,[0:10:100],'Normalization','probability','FaceColor',[139, 0, 139]/255)
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

%test against uniform distribition (count of common neurons in each bin)
[h,p]= chi2gof(mean_common, 'Edges', [0:10:100],'Expected',ones(1,10)*(size(mean_common,2)/10))

%% Plot scatter as a function of place field center difference - account for edge effects

%field difference B-A
field_pos_diff = global_near.B - global_near.A;

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation
%shifting behind track edge
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;

%function of mean position 
%mean_pos = mean([global_near.A; global_near.B],1); 

%convert to normalized space, assign as input variables
x = (global_near.A./100)';
y = (field_pos_diff./100)';

figure;
hold on
scatter(x,  y',14,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)


%% Scatter - Try this for global remapping neurons


%field difference B-A
field_pos_diff = global_far.B - global_far.A;

%merge global positions
global_far_merge = [global_far.A; global_far.B];

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation - no account 
%shifting behind track edge
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%function of mean position 
mean_pos = mean([global_far.A; global_far.B],1); 

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;


%convert to normalized space, assign as input variables
x = (global_far.A./100)';
y = (field_pos_diff./100)';

figure;
hold on
scatter(x,  y',14,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
%max position difference line for global
plot([0 1],[0.15 0.15],'r--')
plot([0 1],[-0.15 -0.15],'r--')
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)

%% Combine all 3 classes

%merge global positions
all_far_merge = [[global_far.A; global_far.B], [global_near.A; global_near.B], common, rate];

%field difference B-A
field_pos_diff = all_far_merge(2,:) -  all_far_merge(1,:);

%global_near.B(track_fwd_idx(1)) - global_near.A(track_fwd_idx(1))

%account for edge effect on difference calculation - no account 
%shifting behind track edge
if 0 
track_rev_idx = find(field_pos_diff >50);
%shifting forward of track edge
track_fwd_idx = find(field_pos_diff <-50);

%account for edge effect
field_pos_diff(track_rev_idx) = field_pos_diff(track_rev_idx)-100;
field_pos_diff(track_fwd_idx) = field_pos_diff(track_fwd_idx)+100;
end

%convert to normalized space, assign as input variables
x = (all_far_merge(1,:)./100)';
y = (field_pos_diff./100)';



figure;
hold on
scatter(x,  y',9,'MarkerFaceColor',[139,0,139]/255,'MarkerEdgeColor',[139,0,139]/255)
axis square

%fit line
%x = global_near.A';                                        % Create Data
%y = field_pos_diff';               % Create Data
mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals
[Ypred,YCI] = predict(mdl, x);                      % Fitted Regression Line & Confidence Intervals

Xnew = linspace(min(x), max(x), 1000);
[~,YCI] = predict(mdl, Xnew'); 

%[B,BINT,R,RINT,STATS] = regress(y,[x,ones(size(x,1),1)])

%plot no correlation 0 line
plot([0 1],[0 0],'LineStyle',':','Color',[34, 139, 34]/255,'LineWidth',2)
plot(Xnew, YCI,'LineStyle', '--','Color',[0.5 0.5 0.5],'LineWidth',1)

plot(x, Ypred,'-k','LineWidth',2)
%max position difference line for global
%plot([0 1],[0.15 0.15],'r--')
%plot([0 1],[-0.15 -0.15],'r--')
xlabel('Normalized position relative to A')
ylabel('B - A field position difference')
set(gca,'FontSize',18)
set(gca,'LineWidth',2)
hold off

%pearson R correlation
[r,p] = corrcoef(x,  y)




%% plot histograms
figure;
hold on
histogram(global_near.A(less_than_unity),[0:20:100])
histogram(global_near.A(greater_than_unity),[0:20:100])



[h,p] =kstest2(global_near.A(less_than_unity),global_near.B(less_than_unity))

scatter(global_far_comb(1,:),global_far_comb(2,:))

figure
subplot(5,1,1)
hold on
ylim([0 ylim_def])
histogram(mean(common,1),[0:bin_space:100],'Normalization','probability')
subplot(5,1,2)
hold on
ylim([0 ylim_def])
histogram(global_far_comb(1,:),[0:bin_space:100],'Normalization','probability')
subplot(5,1,3)
hold on
ylim([0 ylim_def])
histogram(global_far_comb(2,:),[0:bin_space:100],'Normalization','probability')
subplot(5,1,4)
hold on
ylim([0 ylim_def])
histogram(partial_common,[0:bin_space:100],'Normalization','probability')
subplot(5,1,5)
hold on
ylim([0 ylim_def])
histogram(partial.com_parA(:,2),[0:bin_space:100],'Normalization','probability')
histogram(partial.com_parB(:,2),[0:bin_space:100],'Normalization','probability')

%% Ecdf of values above
figure;
hold on
ecdf(mean(common,1))
ecdf(mean(rate,1))
ecdf(global_far_comb(1,:))
ecdf(global_far_comb(2,:))
% ecdf(partial_common)
% ecdf(partial.com_parA(:,2))
% ecdf(partial.com_parB(:,2))

figure;
hold on
histogram(binCenter_data{1}.bin_center.global_far(1,:),0:10:100)

figure;
hold on
histogram(binCenter_data{1}.bin_center.global_far(2,:),0:10:100)

figure
hold on
title('Far global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.global_far_sortA(1,:),binCenter_data{1}.bin_center.global_far_sortA(2,:))

figure
hold on
title('Near global remappers - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.global_near(1,:),binCenter_data{1}.bin_center.global_near(2,:))
plot([0 100],[0 100],'k--')

figure
hold on
title('Common - center of A vs center of B')
xlabel('A field center')
ylabel('B field center')
scatter(binCenter_data{1}.bin_center.common(1,:),binCenter_data{1}.bin_center.common(2,:))
plot([0 100],[0 100],'k--')

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

