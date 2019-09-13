function [outputArg1,outputArg2] = centroid_difference(path_dir)


%read relevant data - A&B tuned (TS)
for ee=1:size(path_dir,2)
    load_data_path_cent_diff{ee} = fullfile(path_dir{ee},'cumul_analysis','centroid_diff.mat');
    centroid_diff_data{ee} = load(string(load_data_path_cent_diff{ee}));
end

%read relevant data - all
for ee=1:size(path_dir,2)
    load_data_path_cent_diff_all{ee} = fullfile(path_dir{ee},'cumul_analysis','centroid_diff_all.mat');
    centroid_diff_data_all{ee} = load(string(load_data_path_cent_diff_all{ee}));
end

%load radian difference for each neuron into 1 vector
%load into cell and convert to vector
for ee=1:size(path_dir,2)
    centroid_diffs{ee} = centroid_diff_data{ee}.cent_diff_AandB.angle_diff;
    centroid_bins{ee} = centroid_diff_data{ee}.cent_diff_AandB.max_bin;
end

%centroid bins and diffs for all neurons
for ee=1:size(path_dir,2)
    centroid_diffs_all{ee} = centroid_diff_data_all{ee}.cent_diff.angle_diff;
    centroid_bins_all{ee} = centroid_diff_data_all{ee}.cent_diff.max_bin;
end

%read relevant data - A&B tuned (TS)
for ee=1:size(path_dir,2)
    load_data_path_select_ROI{ee} = fullfile(path_dir{ee},'cumul_analysis','select_ROI_criteria.mat');
    select_ROI_data{ee} = load(string(load_data_path_select_ROI{ee}));
end

%% Generate logical vector that checks for min 5 events and screen for mutual TS/SI tuning
%for all sessions
for ss=1:size(path_dir,2)
    %make sure at least 1 field is significant in each trial
    for tt=1:2
        sig_fields{ss}{tt} = cellfun(@(x) sum(x,2),select_ROI_data{ss}.select_fields{1}{tt},'UniformOutput',false);
        sig_fields_emp{ss}{tt} = cellfun(@isempty,sig_fields{ss}{tt});
        %convert empty cells to 0
        sig_fields{ss}{tt}(sig_fields_emp{ss}{tt}) = num2cell(zeros(1,length(find(sig_fields_emp{ss}{tt} == 1))));
        %convert to vector
        sig_fields_vec{ss}{tt} = cell2mat(sig_fields{ss}{tt});
    end
    %neuron with significant fields in both trial types
    sig_fields_log{ss} = sig_fields_vec{ss}{1} & sig_fields_vec{ss}{2};
end

%and the sig_fields_vector with the A&B either criteria tuning vector for
%each session
%for all sessions
for ss=1:size(path_dir,2)
    select_ROI_for_centroid{ss} = select_ROI_data{ss}.tunedLogical.tuned_AB_si_ts & sig_fields_log{ss};
end

%% Extract centroid diff for SI+TS A&B tuned neurons 
%empty vector to hold all the centroid diffs
centroid_diff_all_cumul = [];
centroid_diff_all_bins_cumul.A = [];
centroid_diff_all_bins_cumul.B = [];

for ss=1:size(path_dir,2)
    centroid_diff_all_cumul = [centroid_diff_all_cumul centroid_diffs_all{ss}(select_ROI_for_centroid{ss})];
    centroid_diff_all_bins_cumul.A = [centroid_diff_all_bins_cumul.A centroid_bins_all{ss}(1,select_ROI_for_centroid{ss})];
    centroid_diff_all_bins_cumul.B = [centroid_diff_all_bins_cumul.B centroid_bins_all{ss}(2,select_ROI_for_centroid{ss})];
end

%% Downsample into 10 bins
bin_idx.A = discretize(centroid_diff_all_bins_cumul.A,[0:10:100]);
bin_idx.B = discretize(centroid_diff_all_bins_cumul.B,[0:10:100]);

%for each bin get the centroid diff values
for bb=1:10
    cent_diff_binned.A{bb} = centroid_diff_all_cumul(find(bin_idx.A == bb));
    cent_diff_binned.B{bb} = centroid_diff_all_cumul(find(bin_idx.B == bb));
end

%get median value for each bin
figure;
hold on
plot(cellfun(@median,cent_diff_binned.A),'b')
plot(cellfun(@median,cent_diff_binned.B),'r')

for bb=1:10
    [p(bb),h] = ranksum(cent_diff_binned.A{bb},cent_diff_binned.B{bb});
end

%% Display cumulative histogram for SI+TS A&B tuned neurons 
figure
hold on
histogram(centroid_diff_all_cumul,(0:2*pi/32:pi),'FaceColor',[139, 0, 139]/255);
ylabel('Neuron count');
xlabel('Angular difference [rad]');
xticks([pi/4 pi/2 3*pi/4,pi])
xticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
title('Centroid difference')
set(gca,'FontSize',16)

%% Scatterplot of max centroid diff for SI+TS A&B tuned neurons 

% %find minimum bin for each neuron and combine with corresponding ang diff 
% for ee=1:size(path_dir,2)
%     %min_bins_cent_diff{ee} = min(centroid_bins{ee},[],1);
%     min_bins_cent_diff{ee} = centroid_bins{ee}(2,:);
%     %relative to A trials
%     
%     min_bins_cent_diff{ee}(2,:) = centroid_diffs{ee};
% end
% 
% %convert coupled min with cent diff into 1 matrix (from cells)
% combined_min_bins_cent_diff = cell2mat(min_bins_cent_diff);

%plot (cumlative) scatter as a fxn of minimum bin location of max transient
figure('Position', [2060 380 770 580]);
subplot(1,2,1)
hold on
set(gca,'FontSize',16);
title('Positional centroid difference')
scatter(centroid_diff_all_bins_cumul.A,centroid_diff_all_cumul,'filled','SizeData',10,'MarkerFaceColor',[139,0,139]/255)
xlabel('Spatial bin')
ylabel('Centroid difference [rad]');
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})

subplot(1,2,2)
hold on
set(gca,'FontSize',16);
title('Positional centroid difference')
scatter(centroid_diff_all_bins_cumul.B,centroid_diff_all_cumul,'filled','SizeData',10,'MarkerFaceColor',[139,0,139]/255)
xlabel('Spatial bin')
ylabel('Centroid difference [rad]');
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})

%% For boxplots as function of spatial bin
%segregate centroid diffs into 100 bins for EACH animal
%for each animal
for ss=1:size(path_dir,2)
    %for each bin
    for nbin=1:100
        idx_log_nbin = min_bins_cent_diff{ss}(1,:) == nbin;
        boxplot_counts_each{ss}{nbin} =  min_bins_cent_diff{ss}(2,idx_log_nbin);
    end
end

%get median at each bin
for ss=1:size(path_dir,2)
    %for each bin
    %for nbin=1:100
        med_bin(ss,:) = cellfun(@nanmedian,boxplot_counts_each{ss});
    %end
end

%get std fof medians
std_med_bin = nanstd(med_bin,0,1);

%make anonymous functions with parameters as input to shaded error bar func
med_func = @(x) nanmedian(x,1);
mad_func = @(x) mad(x,1,1); %using median along rows

%line plot with std at each spatial bin
figure('Position', [2220 270 730 420]);
plot(nanmedian(med_bin,1))
%shaded std
s = shadedErrorBar(1:100,med_bin,{med_func,mad_func},'lineprops','-','transparent',true,'patchSaturation',0.20);
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
ylim([0 pi])
title('Positional centroid difference');
xlabel('Spatial bin');
ylabel('Centroid difference [rad]');

%group centroid diff by bin (across all)
for nbin=1:100
    idx_log_nbin = combined_min_bins_cent_diff(1,:) == nbin;
    boxplot_counts{nbin} = combined_min_bins_cent_diff(2,idx_log_nbin);
end

%warpper for boxplot
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%get median

%plot the boxplots
f= figure('Position',[1940 400 1250 420]);
boxplot2(boxplot_counts);
hold on
set(gca,'TickLabelInterpreter', 'tex');
ylim([0 pi])
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
xticks([1,20:20:100])
xticklabels({'1','20','40','60','80','100'})
xlabel('Spatial bin')
ylabel('Centroid difference [rad]');
title('Positional centroid distbutions (cumulative)')


%% Max centroid diff only for TS tuned
%convert to vector
centroid_diff_mat = cell2mat(centroid_diffs);

%display as cumulative histogram
figure
hold on
histogram(centroid_diff_mat,(0:2*pi/32:pi),'FaceColor',[1 0 1]);
ylabel('Count');
xlabel('Angular difference [rad]');
xticks([pi/4 pi/2 3*pi/4,pi])
xticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
title('Centroid difference')

%find minimum bin for each neuron and combine with corresponding ang diff 
for ee=1:size(path_dir,2)
    %min_bins_cent_diff{ee} = min(centroid_bins{ee},[],1);
    min_bins_cent_diff{ee} = centroid_bins{ee}(2,:);
    %relative to A trials
    
    min_bins_cent_diff{ee}(2,:) = centroid_diffs{ee};
end

% figure
% hold on
% for ee=1:11
% scatter3(centroid_bins{ee}(1,:), centroid_bins{ee}(2,:),abs(centroid_bins{ee}(1,:)-centroid_bins{ee}(2,:)))
% end

%convert coupled min with cent diff into 1 matrix (from cells)
combined_min_bins_cent_diff = cell2mat(min_bins_cent_diff);

%segregate centroid diffs into 100 bins for EACH animal
%for each animal
for ss=1:size(path_dir,2)
    %for each bin
    for nbin=1:100
        idx_log_nbin = min_bins_cent_diff{ss}(1,:) == nbin;
        boxplot_counts_each{ss}{nbin} =  min_bins_cent_diff{ss}(2,idx_log_nbin);
    end
end

%get median at each bin
for ss=1:size(path_dir,2)
    %for each bin
    %for nbin=1:100
        med_bin(ss,:) = cellfun(@nanmedian,boxplot_counts_each{ss});
    %end
end

%get std fof medians
std_med_bin = nanstd(med_bin,0,1);

%make anonymous functions with parameters as input to shaded error bar func
med_func = @(x) nanmedian(x,1);
mad_func = @(x) mad(x,1,1); %using median along rows

%line plot with std at each spatial bin
figure('Position', [2220 270 730 420]);
plot(nanmedian(med_bin,1))
%shaded std
s = shadedErrorBar(1:100,med_bin,{med_func,mad_func},'lineprops','-','transparent',true,'patchSaturation',0.20);
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
ylim([0 pi])
title('Positional centroid difference');
xlabel('Spatial bin');
ylabel('Centroid difference [rad]');


%construct symmetric matrix for connectogram plot

%plot (cumlative) scatter as a fxn of minimum bin location of max transient
figure('Position', [2060 380 770 580]);
hold on
set(gca,'FontSize',16);
title('Positional centroid difference')
scatter(combined_min_bins_cent_diff(1,:),combined_min_bins_cent_diff(2,:),'filled','SizeData',10,'MarkerFaceColor',[139,0,139]/255)
xlabel('Spatial bin')
ylabel('Centroid difference [rad]');
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})

%group centroid diff by bin (across all)
for nbin=1:100
    idx_log_nbin = combined_min_bins_cent_diff(1,:) == nbin;
    boxplot_counts{nbin} = combined_min_bins_cent_diff(2,idx_log_nbin);
end

%warpper for boxplot
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%get median

%plot the boxplots
f= figure('Position',[1940 400 1250 420]);
boxplot2(boxplot_counts);
hold on
set(gca,'TickLabelInterpreter', 'tex');
ylim([0 pi])
yticks([pi/4 pi/2 3*pi/4,pi])
yticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
xticks([1,20:20:100])
xticklabels({'1','20','40','60','80','100'})
xlabel('Spatial bin')
ylabel('Centroid difference [rad]');
title('Positional centroid distbutions (cumulative)')

end

