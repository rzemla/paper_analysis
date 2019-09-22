function [SCE] = assign_SCE_trials(session_vars,SCE,options)

for ss=options.sessionSelect
    
    %get frames time of first idx of onset
    abs_SCE_idx{ss} = SCE{ss}.sync_idx(SCE{ss}.sync_range(:,1));
    
    %assignmnet log 1- corr A; 2 - corr B; 10 - wrong A; 20 - wrong B
    sce_assign_log{ss} = zeros(size(abs_SCE_idx{ss},1),1);
    
    SCE_start_time{ss} = session_vars{ss}.Imaging_split{3}.time_restricted(abs_SCE_idx{ss});
    
    %all A and B times
    [~,Aall_idx{ss},~] = intersect(SCE_start_time{ss},session_vars{ss}.Imaging_split{4}.time_restricted);
    [~,Ball_idx{ss},~] = intersect(SCE_start_time{ss},session_vars{ss}.Imaging_split{5}.time_restricted);
    %correct A and B times
    [~,Acorr_idx{ss},~] = intersect(SCE_start_time{ss},session_vars{ss}.Imaging_split{1}.time_restricted);
    [~,Bcorr_idx{ss},~] = intersect(SCE_start_time{ss},session_vars{ss}.Imaging_split{2}.time_restricted);
    
    %% Assign to common matrix
    %assign all first as wrong, the correct wrong with correct find
    sce_assign_log{ss}(Aall_idx{ss}) = 10;
    sce_assign_log{ss}(Ball_idx{ss}) = 20;
    
    sce_assign_log{ss}(Acorr_idx{ss}) = 1;
    sce_assign_log{ss}(Bcorr_idx{ss}) = 2;

%% Assign to output SCE struct

SCE{ss}.sce_assign_log = sce_assign_log{ss};

end

%% Find on which lap each SCE occurred

ss=1
%find which lap each SCE belongs to
SCE_start_time{ss}(16)

%% Plot sample SCE dF/F

figure;
for ii=4:79
    select_SCE = ii;
    
    %time range to show for dF/F
    time_range = [SCE{ss}.sync_idx(SCE{ss}.sync_range(select_SCE,1),1)-1000:SCE{ss}.sync_idx(SCE{ss}.sync_range(select_SCE,1),1)+1500]
    
    imagesc(session_vars{ss}.Imaging.trace_restricted(time_range,SCE{ss}.SCE_unique_ROIs{select_SCE})')
    hold on
    colormap('jet')
    caxis([0 2])
    pause;
    clf
end

%line plot here and SCE sorter

figure;

for ii=10:247
    h1 =subplot(2,1,1)
    hold on
    select_SCE = ii;
    
    %time range to show for dF/F
    time_range = [SCE{ss}.sync_idx(SCE{ss}.sync_range(select_SCE,1),1)-1000:SCE{ss}.sync_idx(SCE{ss}.sync_range(select_SCE,1),1)+1500];
    
    hold on
    stepSize = 2;
    step = 0;
    for rr=1:size(SCE{ss}.SCE_unique_ROIs{select_SCE},2)
        plot(session_vars{ss}.Imaging.trace_restricted(time_range,SCE{ss}.SCE_unique_ROIs{select_SCE}(rr))'-step, 'k', 'LineWidth', 1.5);
        step = step - stepSize;
    end
    %red line to mark start of SCE
    plot([1000,1000],[0,2*rr+8],'r--')

    h2 =  subplot(2,1,2)
    hold on
    plot(session_vars{ss}.Behavior.speed(time_range))
    %red line to mark start of SCE
    plot([1000,1000],[-5 20],'r--')
    pause
    cla(h1)
    cla(h2)
end

%%
end

