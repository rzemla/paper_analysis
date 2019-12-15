function [mean_bin_speed,lap_bin_split] = task_sel_speed(tunedLogical,task_selective_ROIs,session_vars,ROI_idx_tuning_class,options)


sessionSelect = options.sessionSelect;

%% Select trial tuned classes of neurons (use logicals as input)
%prefiltered for min 5 events on distinct laps

%S.I.
%for each session
for ss =1:size(session_vars,2)
    %spatial information criterion - regardless if tuned in other session
    Atuned.si{ss} = ROI_idx_tuning_class.si.log.Aonly | ROI_idx_tuning_class.si.log.AB;
    Btuned.si{ss} = ROI_idx_tuning_class.si.log.Bonly | ROI_idx_tuning_class.si.log.AB;

    onlyA_tuned.si{ss} = ROI_idx_tuning_class.si.log.Aonly;
    onlyB_tuned.si{ss} = ROI_idx_tuning_class.si.log.Bonly;    
    AandB_tuned.si{ss} = ROI_idx_tuning_class.si.log.AB;
    neither_tuned.si{ss} = ROI_idx_tuning_class.si.log.N;
    
end

%T.S.
for ss =1:size(session_vars,2)
    %regardless if tuned in other session
    Atuned.ts{ss} = ROI_idx_tuning_class.ts.log.Aonly | ROI_idx_tuning_class.ts.log.AB;
    Btuned.ts{ss} = ROI_idx_tuning_class.ts.log.Bonly | ROI_idx_tuning_class.ts.log.AB;

    onlyA_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.Aonly;
    onlyB_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.Bonly;    
    AandB_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.AB;
    neither_tuned.ts{ss} = ROI_idx_tuning_class.ts.log.N;
end


for ss =1:size(session_vars,2)
    %blank all neuron idx logical vector
    all_neurons{ss} = true(size(Atuned.si{ss}));
    %blank logical equal to size of number of neurons
    blank_log = false(1,size(all_neurons{ss},2));
end

%A/B selective logicals; all A, all B, A&B by either criterion
for ss=1:size(session_vars,2)
    
    %A selective and B selective neurons
    A_sel{ss} = blank_log;
    A_sel{ss}(task_selective_ROIs.A.idx) = 1;
    
    B_sel{ss} = blank_log;
    B_sel{ss}(task_selective_ROIs.B.idx) = 1;
    
    %either criterion - A all,B all, A&B all
    Atuned.si_ts{ss} = Atuned.si{ss} | Atuned.ts{ss};
    Btuned.si_ts{ss} = Btuned.si{ss} | Btuned.ts{ss};
    AandB_tuned.si_ts{ss} = Atuned.si_ts{ss} & Btuned.si_ts{ss}; 
end

%% Convert to indices from logicals for input

%selective
A_sel_idx = find(A_sel{1} ==1);
B_sel_idx = find(B_sel{1} ==1);

%all A or B by either criterion or both
%S.I.
Atuned_idx.si = find(Atuned.si{1} ==1);
Btuned_idx.si = find(Btuned.si{1} ==1);
%T.S.
Atuned_idx.ts = find(Atuned.ts{1} ==1);
Btuned_idx.ts = find(Btuned.ts{1} ==1);
%S.I. or T.S
Atuned_idx.si_ts = find(Atuned.si_ts{1} ==1);
Btuned_idx.si_ts = find(Btuned.si_ts{1} ==1);

%AandB - either criterion or both
%S.I.
AandB_tuned_idx.si = find(AandB_tuned.si{1} ==1);
%T.S.
AandB_tuned_idx.ts = find(AandB_tuned.ts{1} ==1);
%S.I. or T.S.
AandB_tuned_idx.si_ts = find(AandB_tuned.si_ts{1} ==1);

%only A or B by either criterion
%S.I.
onlyA_tuned_idx.si = find(onlyA_tuned.si{1} ==1);
onlyB_tuned_idx.si = find(onlyB_tuned.si{1} ==1);
%T.S.
onlyA_tuned_idx.ts = find(onlyA_tuned.ts{1} ==1);
onlyB_tuned_idx.ts = find(onlyB_tuned.ts{1} ==1);


%% Load in run epoch and speed variables

%across entire restricted session
run_epoch_all_laps = session_vars{1}.Behavior.run_ones;

%across every lap
run_epoch_each_lap = session_vars{1}.Behavior_split_lap.run_ones;

%speed across all complete laps (downsampled to match frames)
speed_all_laps = session_vars{1}.Behavior.speed;

%normalized position across all laps
norm_pos_all_laps = session_vars{1}.Behavior.resampled.normalizedposition;

%trial order of each lap
%2 - A lap
%3 - B lap
%added 0 - (i.e. - 20 or 30 is wrong lap is wrong A or B lap,respectively)
trialOrder = session_vars{1}.Behavior.performance.trialOrder;

%lap number of each frame index
lapNb_fr = session_vars{1}.Behavior.resampled.lapNb;

%% Split by lap
max_lap = max(lapNb_fr);

%get start and end lap idxs across (in frame space) for each lap
for ll=1:max_lap
    lap_idxs(ll,1) = find(lapNb_fr == ll,1,'first');
    lap_idxs(ll,2) = find(lapNb_fr == ll,1,'last');
end

%normalized position and speed on each lap
for ll=1:max_lap
   norm_pos_each_lap{ll} = norm_pos_all_laps(lap_idxs(ll,1):lap_idxs(ll,2));  
   speed_each_lap{ll} = speed_all_laps(lap_idxs(ll,1):lap_idxs(ll,2));
end

%convert speed on each lap to vector
cell2mat(speed_each_lap')

%% Extract bin (run) frame on each lap

%correct A lap indices
corr_lap_idx.A = find(trialOrder == 2);
corr_lap_idx.B = find(trialOrder == 3);

%extract bins associated with each frame
trialType = 1;
[mean_bin_speed.A,lap_bin_split.A] = extract_bin_speed(corr_lap_idx.A,session_vars,speed_each_lap,run_epoch_each_lap,trialType);

trialType = 2;
[mean_bin_speed.B,lap_bin_split.B] = extract_bin_speed(corr_lap_idx.B,session_vars,speed_each_lap,run_epoch_each_lap,trialType);

figure
hold on
plot(mean_bin_speed.A','b')
plot(mean_bin_speed.B','r')

%get sem and mean (on each bin) across laps 
mean_speed_A_laps = mean(mean_bin_speed.A,1);
mean_speed_B_laps = mean(mean_bin_speed.B,1);

%get sem
sem_speed_A_laps = std(mean_bin_speed.A,0,1)./sqrt(size(mean_bin_speed.A,1));
sem_speed_B_laps = std(mean_bin_speed.B,0,1)./sqrt(size(mean_bin_speed.B,1));

cell2mat(lap_bin_split.A')

%% Plot with shades (sem)
figure
hold on
title('Bin speed in run epochs - shades = STD')
xlabel('Bin position');
ylabel('Speed [cm/s]');
%plot(mean_speed_A_laps,'b')
%plot(mean_speed_B_laps,'r')
%errorbar(mean_speed_A_laps,sem_speed_A_laps,'b')
%errorbar(mean_speed_B_laps,sem_speed_B_laps,'r')

%create input function to calculate standard error of mean
sem = @(x) std(x,0,1)./sqrt(size(x,1));

s1 = shadedErrorBar(1:100,mean_bin_speed.A,{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s1.edge,'LineWidth',1.5,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
s1.mainLine.LineWidth = 2;
s1.mainLine.Color = [65,105,225]/255;
s1.patch.FaceColor = [65,105,225]/255;

s2 = shadedErrorBar(1:100,mean_bin_speed.B,{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s2.edge,'LineWidth',1.5,'LineStyle','-','Color',[[220,20,60]/255, 0.2]) %last value add transparency value
s2.mainLine.LineWidth = 2;
s2.mainLine.Color = [220,20,60]/255;
s2.patch.FaceColor = [220,20,60]/255;

legend('A laps','B laps','location','southeast');


%% All lap speed and norm position (QC)
figure
hold on
%norm speed
plot(speed_all_laps', 'g')
%norm position
plot(norm_pos_all_laps-1, 'k')


%% Extract bin code below - modified to be generic for extract_bin_speed fxn
%{
% correct A laps (100 bins)
corr_lap_bins.A = session_vars{1}.Place_cell{1}.Bin{8};

%speed for each frame on corr A laps
speed_corr.A = speed_each_lap(corr_lap_idx.A);

%break run bins into cells (by each corr A lap)
%get run epoch on/off first
run_epoch_corr.A = run_epoch_each_lap(corr_lap_idx.A);

%extract 1 only indices (RUN EPOCH ON)
for ll=1:size(run_epoch_corr.A,2)
    run_only_corr.A{ll} = ll*ones(size(find(run_epoch_corr.A{ll}==1),1),1);
end
%convert to vector
run_only_corr_vec.A = cell2mat(run_only_corr.A');

%get indices for extracting bins as laps
for ll=1:size(run_only_corr.A,2)
    lap_bin_idxs.A(ll,1) = find(run_only_corr_vec.A == ll,1,'first');
    lap_bin_idxs.A(ll,2) = find(run_only_corr_vec.A == ll,1,'last');
end    

%split bins into laps - correct A
for ll=1:size(lap_bin_idxs.A,1)
    corr_lap_bin_split.A{ll} = corr_lap_bins.A(lap_bin_idxs.A(ll,1):lap_bin_idxs.A(ll,2));
end

%generate matching speed values for bins for each lap
for ll=1:size(lap_bin_idxs.A,1)
    speed_match_bin.A{ll} = speed_corr.A{ll}(logical(run_epoch_corr.A{ll}));
end

%get mean speed in each bin on each lap
%for each lap
for ll=1:size(lap_bin_idxs.A,1)
    %for each bin
    for bb=1:100
        speed_bin_lap.A{ll}{bb} = speed_match_bin.A{ll}(find(corr_lap_bin_split.A{ll} == bb));
    end 
end

%take mean speed in each bin on each lap
for ll=1:size(lap_bin_idxs.A,1)
    mean_bin_speed.A(ll,:) = cellfun(@mean,speed_bin_lap.A{ll})
end
%}


end

