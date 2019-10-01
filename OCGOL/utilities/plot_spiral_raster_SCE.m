function [outputArg1,outputArg2] = plot_spiral_raster_SCE(plot_raster_vars,session_vars,registered,cat_registered_cell,SCE, options)


%% Load variables

idxMin = plot_raster_vars.idxMin;
r_scaled = plot_raster_vars.r_scaled;
posVectorApprox = plot_raster_vars.posVectorApprox;
x = plot_raster_vars.x;

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% 

%session to look
ss=2


%norm_position for given session
norm_pos = session_vars{ss}.Behavior.resampled.normalizedposition;

%% Collect all run events for each neuron on correct A or B laps
%get normalized position of events, %calclate non zone segregated median

%onset frame (global - 3 related)
for rr=1:size(session_vars{ss}.Events_split{1}.Run.run_onset_offset,2)
    A_onset_frames{rr} = session_vars{ss}.Events_split{selectTrial(1)}.Run.run_onset_offset{rr}(:,1);
    B_onset_frames{rr} = session_vars{ss}.Events_split{selectTrial(2)}.Run.run_onset_offset{rr}(:,1);
end

%get normalized postion for each event (frame of onset)
for rr=1:size(session_vars{ss}.Events_split{1}.Run.run_onset_offset,2)
    A_all_pos{rr} = norm_pos(A_onset_frames{rr});
    B_all_pos{rr} = norm_pos(B_onset_frames{rr});
end

%get median of onsets (not split into zones)
A_pos_med = cellfun(@nanmedian, A_all_pos);
B_pos_med = cellfun(@nanmedian, B_all_pos);


%% Organize to give a number that tells whether frame is correct A or B

%Trial Order
trialOrder = session_vars{ss}.Behavior.performance.trialOrder;

%convert trial order to letter order
for ii=1:size(trialOrder,1)
    if trialOrder(ii) == 2
        trialOrder_text{ii} = 'A';
    elseif trialOrder(ii) == 3
        trialOrder_text{ii} = 'B';
    elseif trialOrder(ii) == 30
         trialOrder_text{ii} = 'Wrong B';
    elseif trialOrder(ii) == 20
        trialOrder_text{ii} = 'Wrong A';
    end
end

%calcium traces - input into get_SCE_order (all trials)
traces = session_vars{ss}.Imaging.trace_restricted;

%total number of laps
lap_count = size(session_vars{ss}.Behavior.performance.trialOrder,1);
%what lap each frame corresponds to
lap_frame = session_vars{ss}.Behavior.resampled.lapNb;

%total number of frames
frame_count = size(traces,1);
%preallocate
frame_trial_assign = zeros(frame_count,1);

%lap number to which each frame belongs
lap_nb_frame = session_vars{ss}.Behavior.resampled.lapNb;

%get binary lap onset vector
lap_start_idx = find(diff(lap_nb_frame) ==1)+1;

%create empty
lap_start_binary = zeros(1,size(session_vars{ss}.Behavior.resampled.lapNb,1));
%populate
lap_start_binary(lap_start_idx) = 1;

%find norm position intervals or reward A and reward B zones
pos_A_start = mean(session_vars{ss}.Behavior.rewards{2}.position_norm);
pos_B_start = mean(session_vars{ss}.Behavior.rewards{1}.position_norm);

%logical vector or A and B rew zone
rew_A_zone_log = norm_pos>=pos_A_start & norm_pos<= (pos_A_start+0.05);
rew_B_zone_log = norm_pos>=pos_B_start & norm_pos<= (pos_B_start+0.05);

%assign trial type
for ll=1:lap_count
    %if A trial
    if session_vars{ss}.Behavior.performance.trialOrder(ll) == 2
        frame_trial_assign(find(lap_frame ==ll)) = 2;
    elseif session_vars{ss}.Behavior.performance.trialOrder(ll) == 3
        frame_trial_assign(find(lap_frame ==ll)) = 3;
    elseif session_vars{ss}.Behavior.performance.trialOrder(ll) == 30
        frame_trial_assign(find(lap_frame ==ll)) = 6;
    elseif session_vars{ss}.Behavior.performance.trialOrder(ll) == 20
        frame_trial_assign(find(lap_frame ==ll)) = 4;
    end
end

%amend if correct or incorrect trial (20 - wrong A; 30 - wrong B)
for ll=1:lap_count
    %if A trial
    if ~session_vars{ss}.Behavior.performance.trialCorrect(ll) == 1
        if session_vars{ss}.Behavior.performance.trialOrder(ll) == 2
            frame_trial_assign(find(lap_frame ==ll)) = 20;
            
        elseif session_vars{ss}.Behavior.performance.trialOrder(ll) == 3
            frame_trial_assign(find(lap_frame ==ll)) = 30;
        end
    end
end 

%find events for each neuron with the 3 lap frame range
%for each ROI in SCE find frame event onset
for cc=1:SCE{ss}.nbSCE  
 sce_nb = cc;
        start_SCE_frame = SCE{ss}.sync_idx(SCE{ss}.sync_range(sce_nb,1));
    
    %get lap # on which SCE occurred
    sce_lap_nb = lap_frame(start_SCE_frame);
    %extract frame range of previous and lap ahead
    surr_lap_frame_idx_range = find(session_vars{ss}.Behavior.resampled.lapNb >=(sce_lap_nb-1) &...
        session_vars{ss}.Behavior.resampled.lapNb <=(sce_lap_nb+1));
    %start and end frame idx
    surr_lap_start_end = [surr_lap_frame_idx_range(1), surr_lap_frame_idx_range(end)];

    st_idx =surr_lap_start_end(1);
    end_idx = surr_lap_start_end(2);
    
    %extract run state fragment
    run_ones{cc}  =session_vars{ss}.Behavior.run_ones(surr_lap_start_end(1):surr_lap_start_end(2));

    %vector with lap starts
    lap_start{cc} = lap_start_binary(surr_lap_start_end(1):surr_lap_start_end(2));

    %reward zones
    rew_A_pos{cc} = rew_A_zone_log(surr_lap_start_end(1):surr_lap_start_end(2));
    rew_B_pos{cc} = rew_B_zone_log(surr_lap_start_end(1):surr_lap_start_end(2));
    
    %extract the lap indices involved in SCE
    laps_involved{cc} = unique(lap_nb_frame(st_idx:end_idx));
    
    %which trial surround the SCE
    trials_text{cc} = trialOrder_text(laps_involved{cc});

    %frame range of each lap
    %preallocate
    abs_lap_idx_st_end{cc} =[];
    for ll=1:size(laps_involved{cc},1)
        abs_lap_idx_st_end{cc}(ll,1) = find(lap_nb_frame == laps_involved{cc}(ll),1,'first');
        abs_lap_idx_st_end{cc}(ll,2) = find(lap_nb_frame == laps_involved{cc}(ll),1,'last');
    end
    
    %GLOBAL EXTRACTION
    %for each ROI within the SCE
    for rc = 1:size(SCE{ss}.SCE_unique_ROIs_sorted{cc},2)
        %get event onsets within the ~3 laps surround (-1/+1) lap around the SCE onset 
        extract_onset_idx_temp = find(session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1)  > st_idx &...
            session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1) < end_idx);
        
        sce_event_onsets{cc}{rc} = session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(extract_onset_idx_temp,1);
    end
    
    %PERI-LAP EXTRACTION
        %for each ROI within the SCE
        for rc = 1:size(SCE{ss}.SCE_unique_ROIs_sorted{cc},2)
            for ll=1:size(laps_involved{cc},1)
                %get event onsets within the ~3 laps surround (-1/+1) lap around the SCE onset
                extract_onset_idx_temp = find(session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1)  > abs_lap_idx_st_end{cc}(ll,1) &...
                    session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1) < abs_lap_idx_st_end{cc}(ll,2));
                
                %absolute frame onsets by lap
                sce_event_onsets_lap{cc}{rc}{ll} = session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(extract_onset_idx_temp,1);
                %get normalized position of onset
                sce_event_onsets_lap_pos{cc}{rc}{ll} = norm_pos(sce_event_onsets_lap{cc}{rc}{ll});
                
            end
        end
end


%get frame onset of every neuron in SCE by lap (nearby 3 laps)
%session_vars{ss}.Events_split{1}.Run.run_onset_offset{1, 1}

%% Check which SCEs exhibit replay (median of all positions on each lap)

%if p less than 0.05 and positive --> forward replay; 
%if spearman correlation negative --> reverse replay;

%get median onset position on each lap in surrounding laps (all event
%positions
tic;
%preallocate for speed
rho = cell(1, SCE{ss}.nbSCE);
p = cell(1, SCE{ss}.nbSCE);

for cc=1:SCE{ss}.nbSCE
    %for each neuron
    for rc=1:size(sce_event_onsets_lap_pos{cc},2)
        %for each lap
        median_onset_pos{cc}(rc,:) = cellfun(@nanmedian,sce_event_onsets_lap_pos{cc}{rc});
        
    end
    %correlate SCE onsets with median position of firing (check replay)
    [rho{cc},p{cc}] =  corr(SCE{ss}.sorted_onset_times{cc}', median_onset_pos{cc} ,'Type','Spearman','rows','complete');
end
toc;


%% Find replay SCEs

p_sig_SCE = cellfun(@(x) sum(x <0.05),p);
sig_SCEs = find(p_sig_SCE ~= 0)
sigRho = rho(sig_SCEs)
sig_p = p(sig_SCEs)


%% Check which SCEs exhibit replay (median of split positions on each lap)
%start - first reward zone
%first reward zone to second reward zone
%second reward zone to end
%create matrix with norm positions for each animal

reward_B_pos = nanmedian(session_vars{ss}.Behavior.rewards{1}.position_norm);
reward_A_pos= nanmedian(session_vars{ss}.Behavior.rewards{2}.position_norm);

%position ranges
pos_range = [0 reward_B_pos reward_A_pos 1];

%if p less than 0.05 and positive --> forward replay; 
%if spearman correlation negative --> reverse replay;

%get median onset position on each lap in surrounding laps (all event
%positions
tic;
%preallocate for speed
rho = cell(1, SCE{ss}.nbSCE);
p = cell(1, SCE{ss}.nbSCE);

zones_idx = cell(1, SCE{ss}.nbSCE);
zone_pos = cell(1, SCE{ss}.nbSCE);

tic;
for cc=1:SCE{ss}.nbSCE
    %for each neuron
    for rc=1:size(sce_event_onsets_lap_pos{cc},2)
        %for each lap
        %median_onset_pos{cc}(rc,:) = sce_event_onsets_lap_pos{cc}{rc};
        %zone I
        zones_idx{cc}{rc}{1} = cellfun(@(x) find(x >= 0 & x<reward_B_pos), sce_event_onsets_lap_pos{cc}{rc},'UniformOutput',false);
        %zone II
        zones_idx{cc}{rc}{2} = cellfun(@(x) find(x >= reward_B_pos & x<reward_A_pos), sce_event_onsets_lap_pos{cc}{rc},'UniformOutput',false);
        %zone III
        zones_idx{cc}{rc}{3} = cellfun(@(x) find(x >= reward_A_pos & x <1), sce_event_onsets_lap_pos{cc}{rc},'UniformOutput',false);
    end
end
toc;

tic;
%get all position onsets in each zone of track
for cc=1:SCE{ss}.nbSCE
    %for each neuron
    for rc=1:size(sce_event_onsets_lap_pos{cc},2)
        %for each lap
        for ll=1:size(sce_event_onsets_lap_pos{cc}{rc},2)
            zones_pos{cc}{rc}{1}{ll} = sce_event_onsets_lap_pos{cc}{rc}{ll}(zones_idx{cc}{rc}{1}{ll});
            zones_pos{cc}{rc}{2}{ll} = sce_event_onsets_lap_pos{cc}{rc}{ll}(zones_idx{cc}{rc}{2}{ll});
            zones_pos{cc}{rc}{3}{ll} = sce_event_onsets_lap_pos{cc}{rc}{ll}(zones_idx{cc}{rc}{3}{ll});
        end
    end
end
toc;

%get median of onset in each zone
for cc=1:SCE{ss}.nbSCE
    %for each neuron
    for rc=1:size(sce_event_onsets_lap_pos{cc},2)
        %zone I
        med_pos_zone{cc}{rc}{1} = cellfun(@nanmedian, zones_pos{cc}{rc}{1});
        %zone II
        med_pos_zone{cc}{rc}{2} =cellfun(@nanmedian, zones_pos{cc}{rc}{2});
        %zone III
        med_pos_zone{cc}{rc}{3} =cellfun(@nanmedian, zones_pos{cc}{rc}{3});
    end
end

%combine into 1 matrix SCE (all ROI with med position in each zone)
for cc=1:SCE{ss}.nbSCE
    for rc=1:size(sce_event_onsets_lap_pos{cc},2)
        %zone I
        med_pos_comb{cc}{1}(rc,:) = med_pos_zone{cc}{rc}{1};
        %zone II
        med_pos_comb{cc}{2}(rc,:) = med_pos_zone{cc}{rc}{2};
        %zone III
        med_pos_comb{cc}{3}(rc,:) = med_pos_zone{cc}{rc}{3};
    end
end

%[rho_zone{cc},p_zone{cc}] =  corr(SCE{ss}.sorted_onset_times{cc}', med_pos_comb{cc}{2},'Type','Spearman','rows','complete')

 %[rho_test, p_test] =corr(SCE{ss}.sorted_onset_times{cc}', med_pos_comb{cc}{2},'Type','Spearman','rows','complete')

 
 %% Plot individual SCE as scatter on lines across all laps (A vs. B trials)
%number of distinct colors equal to the number of ROIs in SCE
cmap_distinct = distinguishable_colors(70);
figure
for cc=369%sig_SCEs%1:SCE{ss}.nbSCE
    %absolute start frame of sync event
    start_SCE_frame = SCE{ss}.sync_idx(SCE{ss}.sync_range(cc,1));
    
    start_SCE_pos = norm_pos(start_SCE_frame);
    
    
    for tt= 1:2
        subplot(1,2,tt)
        hold on
        %title(trials_text{cc}{ll})
        xlim([-0.2 1.2])
        ylim([(-1)*size(SCE{ss}.SCE_unique_ROIs_sorted{cc},2) 0])
        %for each ROI in SCE
        step = -1;
        for rc=1:size(SCE{ss}.SCE_unique_ROIs_sorted{cc},2)
            %plot flat line
            %plot([0 1],[step step],'k-','LineWidth',1)
            %plot scatter of points along line
            scatter(A_all_pos{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)},ones(1,size(A_all_pos{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)},1)).*step,...
                21,'filled','MarkerFaceColor',cmap_distinct(rc,:))
            step = step -1;
        end
        %add SCE position
        plot([start_SCE_pos, start_SCE_pos],[0 step+1],'k--')
    end
    %pause
    %clf
end
 
%% Plot individual SCE as scatter on lines across neighboring laps
%number of distinct colors equal to the number of ROIs in SCE
cmap_distinct = distinguishable_colors(70);
figure
for cc=369%sig_SCEs%1:SCE{ss}.nbSCE
    %absolute start frame of sync event
    start_SCE_frame = SCE{ss}.sync_idx(SCE{ss}.sync_range(cc,1));
    
    start_SCE_pos = norm_pos(start_SCE_frame);
    
    
    for ll= 1:size(sce_event_onsets_lap_pos{cc}{1},2)
        subplot(1,3,ll)
        hold on
        title(trials_text{cc}{ll})
        xlim([-0.2 1.2])
        ylim([(-1)*size(sce_event_onsets_lap_pos{cc},2) 0])
        %for each ROI in SCE
        step = -1;
        for rc=1:size(sce_event_onsets_lap_pos{cc},2)
            %plot flat line
            %plot([0 1],[step step],'k-','LineWidth',1)
            %plot scatter of points along line
            scatter(sce_event_onsets_lap_pos{cc}{rc}{ll},ones(1,size(sce_event_onsets_lap_pos{cc}{rc}{ll},1)).*step,...
                21,'filled','MarkerFaceColor',cmap_distinct(rc,:))
            step = step -1;
        end
        %add SCE position
        plot([start_SCE_pos, start_SCE_pos],[0 step+1],'k--')
    end
    %pause
    %clf
end

%% Figure out lap that event occurred on
figure;
for cc=14
    sce_nb = cc;
    
    %absolute time frame of SCE start (all trials)
    start_SCE_frame = SCE{ss}.sync_idx(SCE{ss}.sync_range(sce_nb,1));
    
    %get lap # on which SCE occurred
    sce_lap_nb = lap_frame(start_SCE_frame);
    %extract frame range of previous and lap ahead
    surr_lap_frame_idx_range = find(session_vars{ss}.Behavior.resampled.lapNb >=(sce_lap_nb-1) &...
        session_vars{ss}.Behavior.resampled.lapNb <=(sce_lap_nb+1));
    %start and end frame idx
    surr_lap_start_end = [surr_lap_frame_idx_range(1), surr_lap_frame_idx_range(end)];
    
    %frame range for plotting
    frame_range = diff(surr_lap_start_end);
    
    %start and end laps idx
    st_idx =surr_lap_start_end(1);
    end_idx = surr_lap_start_end(2);

    %get lap start indices
    lap_start_idxs =find(diff(frame_trial_assign(st_idx:end_idx)) ==1 | diff(frame_trial_assign(st_idx:end_idx)) == -1)+1;

    lap_start_idx2 = find(diff(lap_start_binary(st_idx:end_idx))==1)+1;

    % Plot example SCE - move to separate script
    
    %temporary generate times and speed for select input inverval
    time = session_vars{ss}.Imaging.time_restricted;
    speed = session_vars{ss}.Behavior.speed;
    pos_norm = session_vars{ss}.Behavior.resampled.normalizedposition;
    events = session_vars{ss}.Events.Run.run_onset_binary;
    
    %# of frames before and after onset of SCE
    %plot_range = [1000, 1000];
    
    %start and end points of sorted (10-15 frame range) - 500 ms = 15
    st_evt_sort = start_SCE_frame;
    % st_evt_sort = 11511;
    end_evt_sort = st_evt_sort+15;
    
    
    st_idx =surr_lap_start_end(1);
    end_idx = surr_lap_start_end(2);
    
    subplot(3,2,1)
    imagesc(traces(st_idx:end_idx,SCE{ss}.SCE_unique_ROIs_sorted{sce_nb})')
    hold on;
    xlim([0 frame_range])
    title(['SCE ROIs traces sorted by onset time ','SCE: ',num2str(cc)])
    colormap('hot')
    caxis([0 1])
    
    subplot(3,2,3)
    hold on
    ylabel('Normalized position')
    plot(pos_norm(st_idx:end_idx),'k-','LineWidth',2)
    %plot trial type
    plot(frame_trial_assign(st_idx:end_idx)*0.1)

    xlim([0 frame_range])
    %plot frames where sync event occurred
    %plot([plot_range(1),plot_range(1)],[0 1],'Color',[0.5 0.5 0.5],'LineStyle', '--')
    
    %color blue A trials
    %fr_assign_cut_A  = pos_norm(st_idx:end_idx);
    %plot(find(frame_trial_assign(st_idx:end_idx)==2),fr_assign_cut_A(find(frame_trial_assign(st_idx:end_idx)==2)),'b')
    %color red B trials
    %fr_assign_cut_B  = pos_norm(st_idx:end_idx);
    %plot(find(frame_trial_assign(st_idx:end_idx)==3),fr_assign_cut_B(find(frame_trial_assign(st_idx:end_idx)==3)),'r')
    
    subplot(3,2,5)
    hold on
    plot(speed(st_idx:end_idx),'k','LineWidth',1)
    xlim([0 frame_range])
    ylabel('Speed [cm/s]');
    %plot frames where sync event occurred
    %plot([plot_range(1),plot_range(1)],[-5 25],'Color',[0.5 0.5 0.5],'LineStyle', '--')
    
    subplot(3,2,[2,4,6])
    hold on
    %set(gca,'Color','k')
    xlim([0 frame_range])
    stepSize = 2;
    step = 0;
    
    %number of distinct colors equal to the number of ROIs in SCE
    cmap_distinct = distinguishable_colors(size(SCE{ss}.SCE_unique_ROIs_sorted{sce_nb},2));
    
    %first 30 ROIs
    for ii=1:size(SCE{ss}.SCE_unique_ROIs_sorted{sce_nb},2)
        plot(traces(st_idx:end_idx,SCE{ss}.SCE_unique_ROIs_sorted{sce_nb}(ii))'+step, 'Color',0*[1 1 1], 'LineWidth', 1.5)
        %plot individual run events
        scatter(sce_event_onsets{cc}{ii}-st_idx,ones(1,size(sce_event_onsets{cc}{ii},1))+step ,14,'filled','MarkerFaceColor',...
            cmap_distinct(ii,:)) %dark green[0,128,0]/255
        
        %plot run state across all neuron traces
        plot((run_ones{cc}*0.5)+step, '-','Color',[34,139,34]./255,'LineWidth',0.5)
        
        % A reward pos
        plot((rew_A_pos{cc}*0.3)+step, '-','Color','b','LineWidth',0.5)
        plot((rew_B_pos{cc}*0.3)+step, '-','Color','r','LineWidth',0.5) 

        step = step - stepSize;
    end
    
    %label the y with relevant ROIs
    y_ticks = 0:-2:(step+2);
    y_labels = cell(1,size(y_ticks,2));

    for lab_idx=1:size(y_ticks,2)
        y_labels{lab_idx} = num2str(size(y_ticks,2)-lab_idx+1);
    end

    yticks(fliplr(y_ticks));
    yticklabels(y_labels)
    
    %mark 50 frames around SCE with dotted gray line
    plot([start_SCE_frame-st_idx-50,start_SCE_frame-st_idx-50],[step 5],'Color',[1 1 1]*0.5,'LineStyle','--');
    plot([start_SCE_frame-st_idx+50,start_SCE_frame-st_idx+50],[step 5],'Color',[1 1 1]*0.5,'LineStyle','--');
    
    %plot lap start and stop
    %get lap start indices
       
    if size(lap_start_idx2,2) == 1
        plot(repmat(lap_start_idx2,1,size(lap_start_idx2,2)*2)',repmat([step 5]',1,size(lap_start_idx2,2)),'k--')
    else
        plot(repmat(lap_start_idx2',1,size(lap_start_idx2,2))',repmat([step 5]',1,size(lap_start_idx2,2)),'k--')
    end

    %plot B reward zone start
    %plot(repmat(reward_B_3lap_loc,1,2)',repmat([step 5]',1,3),'r--')

    %plot A reward zone start
    %plot(repmat(reward_A_3lap_loc,1,2)',repmat([step 5]',1,3),'b--')

    %plot(repmat(run_on_after_B,1,2)',repmat([step 5]',1,3),'g--')

    %imagesc(events(st_idx:end_idx,SCE{ss}.SCE_unique_ROIs_sorted{sce_nb})')
    %hold on;
    %xlim([0 frame_range])
    %title('All SCE ROIs events sorted by onset time')
    %colormap('hot')
    %caxis([0 1])
    
    pause
    clf;
    
end

%% Plot raster , spiral, ROI FOV across  - 2 session comparison

figure('Position',[2230 30 780 960]);
hold on

%reward A/B vector
reward_A_vector = exp(i*2*pi*0.70);
reward_B_vector = exp(i*2*pi*0.30);
lap_start_vector = exp(i*2*pi*0);

for cc=14
    
    sce_nb = cc;
    %list of ROIs to plot on given session
    sce_rois = SCE{ss}.SCE_unique_ROIs_sorted{sce_nb};
    
    for ii=1:size(sce_rois,2) %with nans where no match
        
        disp(['Trial type: ', num2str(SCE{ss}.sce_assign_log(sce_nb))])
        
        %for each session
        for ss=ss
            
            ROI = sce_rois(ii);
            %skip of nan value
            if ~isnan(ROI)
                
                %spiral plot early in learning
                subplot(5,6,ii)
                polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
                hold on
                %title(cat_registered_cell{ii,ss});
                %plot A (2) trial events
                for ll=1:size(idxMin{ss}{1},2)
                    polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),14,'bo','MarkerFaceColor','b')
                    %place field center
                    %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
                end
                %plot reward and start locations
                polarplot([0+0i,20*reward_A_vector],'b--','LineWidth',2)
                polarplot([0+0i,20*lap_start_vector],'g--','LineWidth',2)
                polarplot([0+0i,20*reward_B_vector],'r--','LineWidth',2)

                %plot tuning specificity vector for all A trials (unfiltered)
                %polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(1)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
                
                %plot B (3) trial events
                for ll=1:size(idxMin{ss}{2},2)
                    polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),14,'ro','MarkerFaceColor','r')
                    %place field center
                    %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
                end
                
                %plot tuning specificity vector for all B trials
                %polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(2)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
                
                hold off
                
            end
            %pause
        end
        %pause
        disp(ii)
    end
    %pause()
    %clf;
end


%% Plot any arbitrary set of of neurons from list

%find SCE with mean TC below 0.3
sce_select = find(SCE{ss}.meanTC < 0.4)
sce_A_select = find(SCE{ss}.sce_assign_log == 3)

select_SCE =intersect(sce_select,sce_A_select)

[4;5;6;9;10;19;20;21;22;32;34;43;66;71;103;105;109;111]

nbROI_A_engage = sum(SCE{ss}.sce_activity.A,2)
nbROI_B_engage = sum(SCE{ss}.sce_activity.B,2)

merge_ROIs_engage = [nbROI_A_engage, nbROI_B_engage]

%reward A/B vector
reward_A_vector = exp(i*2*pi*0.70);
reward_B_vector = exp(i*2*pi*0.30);
lap_start_vector = exp(i*2*pi*0);

sce_nb = 27;
neuron_list =545;

%neuron_list = SCE{ss}.SCE_unique_ROIs_sorted{sce_nb};

figure('Position',[2230 30 780 960]);
hold on
  
subplot_counter = 1;
    for ii=1:size(neuron_list,2) %with nans where no match
        
        disp(['Trial type: ', num2str(SCE{ss}.sce_assign_log(sce_nb)),' ',...
                'SCE position start: ', num2str(SCE{1, 1}.SCE_pos_start(sce_nb))])
        
        %for each session
        for ss=ss
            
            ROI = neuron_list(ii);
            %skip of nan value
            if ~isnan(ROI)
                
                %spiral plot early in learning
                subplot(5,6,subplot_counter)
                polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
                hold on
                %title(cat_registered_cell{ii,ss});
                %plot A (2) trial events
                for ll=1:size(idxMin{ss}{1},2)
                    polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),14,'bo','MarkerFaceColor','b')
                    %place field center
                    %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
                end
                
                %plot tuning specificity vector for all A trials (unfiltered)
                %polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(1)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
                polarplot([0+0i,20*reward_A_vector],'b--','LineWidth',2)
                polarplot([0+0i,20*lap_start_vector],'g--','LineWidth',2)
                %plot B (3) trial events
                for ll=1:size(idxMin{ss}{2},2)
                    polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),14,'ro','MarkerFaceColor','r')
                    %place field center
                    %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
                end
                
                %plot tuning specificity vector for all B trials
                %polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(2)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
                polarplot([0+0i,20*reward_B_vector],'r--','LineWidth',2)
                polarplot([0+0i,20*lap_start_vector],'g--','LineWidth',2)

                hold off
                
            end
            %pause
        end
        %pause
        disp(ii)
subplot_counter = subplot_counter + 1;
    end


%% Plot all ROI events in two spiral
% if 0
% %which session
% ses=1;
% 
% figure('Position',[2230 30 780 960]);
% sce_rois = SCE{ses}.SCE_unique_ROIs_sorted{sce_nb};
% for dd=1:137
% 
%     sce_nb = dd;
%     %list of ROIs to plot on given session
%     for ii=1:size(sce_rois,2) %with nans where no match
%         ROI = sce_rois(ii);
%         disp(['Trial type: ', num2str(SCE{ses}.sce_assign_log(sce_nb))])
%         %s1 = subplot(2,1,1);
%         %spiral plot early in learning
%         polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
%         hold on
%         
%         %title(cat_registered_cell{ii,ss});
%         %plot A (2) trial events
%         for ll=1:size(idxMin{ss}{1},2)
%             polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
%             %place field center
%         end
%         
%         %plot tuning specificity vector for all A trials (unfiltered)
%         polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(1)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
%         
% end
% end 
% 
% 
%         s2= subplot(2,1,2);
%         polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
%         hold on
%         %plot B (3) trial events
%         for ll=1:size(idxMin{ss}{2},2)
%             polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
%             %place field center
%             %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
%         end
%         
%         %plot tuning specificity vector for all B trials
%         polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(2)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
%         
%         hold off
%     end
% 
%     pause
%     clf;
% end
% end





end
