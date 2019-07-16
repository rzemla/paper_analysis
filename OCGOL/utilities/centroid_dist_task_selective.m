function [centroid_ct] = centroid_dist_task_selective(tunedLogical, task_selective_ROIs,max_bin_rate,options)
%distribution of max centroids across length of track (bin space)

%was originally session cell with relevant variables - import from learning
%script
animal_data = 1;

switch options.tuning_criterion
    case 'si' %spatial information
        %for each session
        for ss =1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).si.Atuned;
            Btuned{ss} = tunedLogical(ss).si.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).si.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).si.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).si.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).si.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            
            all_neurons{ss} = true(size(Atuned{ss}));
        end
    case 'ts' %spatial information
        for ss =1:size(animal_data,2)
            %spatial information criterion
            Atuned{ss} = tunedLogical(ss).ts.Atuned;
            Btuned{ss} = tunedLogical(ss).ts.Btuned;
            
            AandB_tuned{ss} =  tunedLogical(ss).ts.AandB_tuned;
            AorB_tuned{ss} = tunedLogical(ss).ts.AorB_tuned;
            onlyA_tuned{ss} = tunedLogical(ss).ts.onlyA_tuned;
            onlyB_tuned{ss} = tunedLogical(ss).ts.onlyB_tuned;
            AxorB_tuned{ss} =  AorB_tuned{ss} & ~AandB_tuned{ss};
            
            all_neurons{ss} = true(size(Atuned{ss}));
        end
        
    case 'selective_filtered'
        %neurons idxs associated with selective filtering for
        %task-selectivity
        select_filt_ROIs.A = task_selective_ROIs.A.idx;
        select_filt_ROIs.B = task_selective_ROIs.B.idx;
end


%% Generate vectors with only max rate from A or B tuned neurons by TS tuning criterion
%only A; only B
bin_rate.A = max_bin_rate{1}{1}(select_filt_ROIs.A);
bin_rate.B = max_bin_rate{1}{2}(select_filt_ROIs.B);

%A and B
%bin_rate.AB{1} = max_bin_rate{1}{1}(AandB_tuned{1, 1});
%bin_rate.AB{2} = max_bin_rate{1}{2}(AandB_tuned{1, 1});

%% Generate count of centroids at each bin (1-100) (all neurons)

%for each trial
for tt=1:2
    [centroid_count(tt,:),edges] = histcounts(max_bin_rate{1}{tt}, 0.5:1:100.5);
end

%% Generate count of centroids at each bin (1-100) (A only and B only neurons)
[centroid_ct.A,edges] = histcounts(bin_rate.A, 0.5:1:100.5);
[centroid_ct.B,edges] = histcounts(bin_rate.B, 0.5:1:100.5);

%AandB
%[centroid_ct.AB{1},edges] = histcounts(bin_rate.AB{1}, 0.5:1:100.5);
%[centroid_ct.AB{2},edges] = histcounts(bin_rate.AB{2}, 0.5:1:100.5);

%% Plot histogram of max rate for each trial
% for all neurons
if 0
figure;
subplot(2,1,1)
hold on
title('A trials')
histogram('BinEdges', edges,'BinCounts', centroid_count(1,:))
subplot(2,1,2)
hold on
title('B trials')
histogram('BinEdges', edges,'BinCounts', centroid_count(2,:))
end

%% Plot histogram selective for A or B only
% top - A only; bottom B onyly
figure;
subplot(2,1,1)
hold on
title('A trials - selective')
ylim([0 10])
histogram('BinEdges', edges,'BinCounts', centroid_ct.A)
subplot(2,1,2)
hold on
title('B trials - selective')
ylim([0 10])
histogram('BinEdges', edges,'BinCounts', centroid_ct.B)

%% Plot histogram selective for A or B only
% top - A only; bottom B onyly
if 0
figure;
subplot(2,1,1)
hold on
title('A trials')
ylim([0 10])
histogram('BinEdges', edges,'BinCounts', centroid_ct.AB{1})
subplot(2,1,2)
hold on
title('B trials')
ylim([0 10])
histogram('BinEdges', edges,'BinCounts', centroid_ct.AB{2})
end
%% Export centroid counts for A, B, and A&B selective counts



end

