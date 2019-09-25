function [outputArg1,outputArg2] = plot_spiral_raster_SCE(plot_raster_vars,session_vars,registered,cat_registered_cell,SCE, options)


%% Load variables

idxMin = plot_raster_vars.idxMin;
r_scaled = plot_raster_vars.r_scaled;
posVectorApprox = plot_raster_vars.posVectorApprox;
x = plot_raster_vars.x;

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;


%% Plot event maps one lap before and lap ahead of SCE got each SCE neuron
%% Organize to give a number that tells whether frame is correct A or B

ss=1
ses=1

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

%assign trial type
for ll=1:lap_count
    %if A trial
    if session_vars{ss}.Behavior.performance.trialOrder(ll) == 2
        frame_trial_assign(find(lap_frame ==ll)) = 2;
    elseif session_vars{ss}.Behavior.performance.trialOrder(ll) == 3
        frame_trial_assign(find(lap_frame ==ll)) = 3;
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
for cc=1:137
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

    for rc = 1:size(SCE{ss}.SCE_unique_ROIs_sorted{cc},2)
        extract_onset_idx_temp = find(session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1)  > st_idx &...
            session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(:,1) < end_idx);
        sce_event_onsets{cc}{rc} = session_vars{ss}.Events.Run.run_onset_offset{SCE{ss}.SCE_unique_ROIs_sorted{cc}(rc)}(extract_onset_idx_temp,1);
    end
end

% %2 - 3 - 3

lap_start_idxs(1)

%define A and B reward zone edges
pos_norm_3_lap_clipped = pos_norm(st_idx:end_idx);
%logical vector of run on state clipped
run_ones_clip = session_vars{ss}.Behavior.run_ones(st_idx:end_idx);

%all onsets
run_onset_loc = find(diff(run_ones_clip) == 1)+1;

%find all idxs within specified position
reward_B_3lap_loc = find(diff(pos_norm_3_lap_clipped >= 0.29 & pos_norm_3_lap_clipped <= 0.30) == 1)+1;
reward_A_3lap_loc = find(diff(pos_norm_3_lap_clipped >= 0.70 & pos_norm_3_lap_clipped <= 0.71) == 1)+1;

%for each reward B location find frame or first run epoch following
for ii=1:size(reward_B_3lap_loc,1)
    first_B_frame_epochs_idxs(ii) = find(run_onset_loc > reward_B_3lap_loc(ii),1);
end

run_on_after_B = run_onset_loc(first_B_frame_epochs_idxs)

%last frame mod for cont run
run_on_after_B(3) = 4114;

%% Figure out lap that event occurred on
figure;
for cc=75:137
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
    
    %get lap start indices
    lap_start_idxs =find(diff(frame_trial_assign(st_idx:end_idx)) ==1 | diff(frame_trial_assign(st_idx:end_idx)) == -1)+1;


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
    title('SCE ROIs traces sorted by onset time')
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
    %first 30 ROIs
    for ii=1:size(SCE{ss}.SCE_unique_ROIs_sorted{sce_nb},2)
        plot(traces(st_idx:end_idx,SCE{ss}.SCE_unique_ROIs_sorted{sce_nb}(ii))'+step, 'Color',0*[1 1 1], 'LineWidth', 1.5)
        %plot individual run events
        scatter(sce_event_onsets{cc}{ii}-st_idx,ones(1,size(sce_event_onsets{cc}{ii},1))+step ,14,'filled','MarkerFaceColor',...
        [0,128,0]/255)  
        step = step - stepSize;
    end
    %mark 50 frames around SCE with dotted gray line
    plot([start_SCE_frame-st_idx-50,start_SCE_frame-st_idx-50],[step 5],'Color',[1 1 1]*0.5,'LineStyle','--');
    plot([start_SCE_frame-st_idx+50,start_SCE_frame-st_idx+50],[step 5],'Color',[1 1 1]*0.5,'LineStyle','--');

    %plot lap start and stop
    %get lap start indices
    plot(repmat(lap_start_idxs,1,2)',repmat([step 5]',1,2),'k--')

    %plot B reward zone start
    %plot(repmat(reward_B_3lap_loc,1,2)',repmat([step 5]',1,3),'r--')

    %plot A reward zone start
    %plot(repmat(reward_A_3lap_loc,1,2)',repmat([step 5]',1,3),'b--')

    plot(repmat(run_on_after_B,1,2)',repmat([step 5]',1,3),'g--')

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

for cc=75
    
    sce_nb = cc;
    %list of ROIs to plot on given session
    sce_rois = SCE{ses}.SCE_unique_ROIs_sorted{sce_nb};
    
    for ii=1:size(sce_rois,2) %with nans where no match
        
        disp(['Trial type: ', num2str(SCE{ses}.sce_assign_log(sce_nb))])
        
        %for each session
        for ss=ses
            
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
        for ss=ses
            
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
