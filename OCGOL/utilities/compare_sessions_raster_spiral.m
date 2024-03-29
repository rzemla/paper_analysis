function [outputArg1,outputArg2] = compare_sessions_raster_spiral(session_vars,registered,ROI_outlines,ROI_zooms,options)

%% Define/load variables for each session

%for each session
for ii = 1:size(session_vars,2)
    % behavior and imaging related variables
    Behavior_split_lap{ii} = session_vars{ii}.Behavior_split_lap;
    Events_split_lap{ii} = session_vars{ii}.Events_split_lap;
    Behavior_split{ii} = session_vars{ii}.Behavior_split;
    Event_split{ii} = session_vars{ii}.Events_split;
    Imaging_split{ii} = session_vars{ii}.Imaging_split;
    Place_cell{ii} = session_vars{ii}.Place_cell;
    Behavior_full{ii} = session_vars{ii}.Behavior;
    
    %all within run domain
    position{ii}  = Behavior_split_lap{ii}.Run.position;
    time{ii} = Behavior_split_lap{ii}.Run.time;
    events_full{ii} = Events_split_lap{ii}.Run.run_onset_binary;
    run_intervals{ii} = Behavior_split_lap{ii}.run_ones;
    
    %global trial type order across restricted laps
    trialOrder{ii} = Behavior_full{ii}.performance.trialOrder;
end

%for each session
for ss=1:size(session_vars,2)
    %for each lap
    for ii=1:size(run_intervals{ss},2)
        events{ss}{ii} = events_full{ss}{ii}(logical(run_intervals{ss}{ii}),:);
    end
end

%lap_times = Behavior_split_lap.lap;
%how many trial types are there
%set to 2 - A and B for now
%trialTypes = 2;

%make this a condition for the type of trials that are compared
%all correct; all regardless of correct

%for each session
for ss = 1:size(session_vars,2)
    %find times from trials related to specific trial type
    %A trials
    trialTypeIdx{ss}{1} = find(trialOrder{ss} == 2 | trialOrder{ss} == 20);
    %B trials
    trialTypeIdx{ss}{2} = find(trialOrder{ss} == 3 | trialOrder{ss} == 30);
end

%% Define the spiral parameters according to the number of laps
%equivalent to number of laps for each session (turns = laps)
for ss = 1:size(session_vars,2)
    turns(ss) = size(trialOrder{ss},1); %The number of turns the spiral will have (how many laps)
    
    %x is the angle
    x{ss} = [-1*pi*turns(ss) : 0.01 : pi*turns(ss)];
    
    %r is that radius of the point
    r{ss} = [0:1/(length(x{ss})-1):1];
    
    %scale to lap length
    r_scaled{ss} = r{ss}.*turns(ss);
end

%all parameters in the run frame domain
%find the frames index of event and position
%for each session
for ss=1:size(session_vars,2)
    %for trial type (A or B)
    for ii=1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll= 1:size(trialTypeIdx{ss}{ii},1)
            %for each ROI
            for rr=1:size(events{ss}{trialTypeIdx{ss}{ii}(ll)},2)
                %event indices in run domain
                event_idx{ss}{ii}{ll}{rr} = find(events{ss}{trialTypeIdx{ss}{ii}(ll)}(:,rr) == 1);
                %position that corresponds to indices
                pos{ss}{ii}{ll}{rr} = position{ss}{trialTypeIdx{ss}{ii}(ll)}(event_idx{ss}{ii}{ll}{rr});
                %position vectors that will be used as input to spiral
                posVectors{ss}{ii}{ll}{rr} = trialTypeIdx{ss}{ii}(ll).*exp(1i.*((pos{ss}{ii}{ll}{rr}/200)*2*pi)).';
            end
        end
    end
end

%make the nearest approximation to a point along the spiral vector defined above
%predefine to avoid empty cells at the end
valMin = posVectors;
idxMin = posVectors;
posVectorApprox = posVectors;

%foreach session
for ss =1:size(session_vars,2)
    %for each trial type
    for kk = 1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll = 1:size(pos{ss}{kk},2)
            %for each ROI
            for rr = 1:size(events{ss}{1},2)
                %for each event
                for ee=1:size(pos{ss}{kk}{ll}{rr},1)
                    %fix empty cell clipping at endpoints
                    [valMin{ss}{kk}{ll}{rr}(ee),idxMin{ss}{kk}{ll}{rr}(ee)] = min(abs( (r_scaled{ss} - ( (trialTypeIdx{ss}{kk}(ll)-1) + (pos{ss}{kk}{ll}{rr}(ee)/200) ) ) ));
                    posVectorApprox{ss}{kk}{ll}{rr}(ee) = r_scaled{ss}(idxMin{ss}{kk}{ll}{rr}(ee))*exp(1i.*(pos{ss}{kk}{ll}{rr}(ee)/200)*2*pi);
                end
            end
        end
    end
end

%% Plot raster, event spiral and matching ROIs from FOV
%modify this to see events serially
if 0
figure('Position',[2600,300,1200,1000]);
for ii=options.idx_show%1:size(registered.multi.assigned_all,1)
    
    %ROI from session 1
    ROI = registered.multi.assigned_all(ii,1);
    
    subplot(2,3,1)
    imagesc(session_vars{1}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
    hold on;
    title(num2str(ROI));
    ylabel('Lap #'); 
    xlabel('Spatial bin');
    caxis([0 2])
    colormap(gca,'jet');
    hold off;
    
    
    
    %spiral plot early in learning
    subplot(2,3,2)
    polarplot(x{1},r_scaled{1},'k','Linewidth',1.5)
    hold on
    
    %plot A (2) trial events
    for ll=1:size(idxMin{1}{1},2)
        polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI}),r_scaled{1}(idxMin{1}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
        %place field center
        %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
    end
    
    %plot tuning specificity vector for all A trials
    polarplot([0+0i,15*session_vars{1, 1}.Place_cell{1, 4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)


    %plot B (3) trial events
    for ll=1:size(idxMin{1}{2},2)
        polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI}),r_scaled{1}(idxMin{1}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
        %place field center
        %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
    end
    
    %plot tuning specificity vector for all B trials
    polarplot([0+0i,15*session_vars{1, 1}.Place_cell{1, 5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)


    hold off
    
    subplot(2,3,3)
    imagesc(ROI_zooms{ii,1})
    hold on;
    colormap(gca, 'gray')
    xticks([])
    yticks([])
    b = bwboundaries(ROI_outlines{ii,1},'noholes');
    plot(b{1}(:,2),b{1}(:,1),'r')
    hold off
    
    %ROI from session 2
    ROI = registered.multi.assigned_all(ii,2);
    subplot(2,3,4)
    imagesc(session_vars{2}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
    hold on;
    title(num2str(ROI));
    ylabel('Lap #'); 
    xlabel('Spatial bin');
    caxis([0 2])
    colormap(gca, 'jet');
    hold off;
    
   
    
    subplot(2,3,5)
    polarplot(x{2},r_scaled{2},'k','Linewidth',1.5)
    hold on
    %plot A (2) trial events
    for ll=1:size(idxMin{2}{1},2)
        polarscatter(angle(posVectorApprox{2}{1}{ll}{ROI}),r_scaled{2}(idxMin{2}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
        %place field center
        %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
    end
    %plot tuning specificity vector for all A trials
    polarplot([0+0i,15*session_vars{2}.Place_cell{4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)

    
    %plot B (3) trial events
    for ll=1:size(idxMin{2}{2},2)
        polarscatter(angle(posVectorApprox{2}{2}{ll}{ROI}),r_scaled{2}(idxMin{2}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
        %place field center
        %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
    end
    
    %plot tuning specificity vector for all A trials
    polarplot([0+0i,15*session_vars{2}.Place_cell{5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)

    hold off
    
    subplot(2,3,6)
    imagesc(ROI_zooms{ii,2})
    %imagesc(ROI_zooms{registered.multi.assigned_all(ii,2),2})
    hold on;
    colormap(gca, 'gray')
    xticks([])
    yticks([])
    %b = bwboundaries(ROI_outlines{registered.multi.assigned_all(ii,2),2},'noholes');
    b = bwboundaries(ROI_outlines{ii,2},'noholes');
    plot(b{1}(:,2),b{1}(:,1),'r')
    hold off
    
    %pause;
    %clf;
    
end

%% Plot raster , spiral, ROI FOV across all 6 session

%find day1 and da

%order of plots
subplot_matrix = 1:18;
subplot_matrix = reshape(subplot_matrix,6,3)';

figure('Position',[1920 40 1920 960]);
for ii=1:size(registered.multi.assigned,1) %with nans where no match
    %size(registered.multi.assigned_all,1) - match on each ses %options.idx_show%
    
    %for each session
    for ss=1:6
        %ROI = registered.multi.assigned_all(ii,ss);
        ROI = registered.multi.assigned(ii,ss);
        %skip of nan value
        if ~isnan(ROI)
            subplot(3,6,subplot_matrix(1,ss))
            imagesc(session_vars{ss}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
            hold on;
            title(num2str(ROI));
            ylabel('Lap #');
            xlabel('Spatial bin');
            caxis([0 2])
            colormap(gca,'jet');
            hold off;
            
            %spiral plot early in learning
            subplot(3,6,subplot_matrix(2,ss))
            polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
            hold on
            
            %plot A (2) trial events
            for ll=1:size(idxMin{ss}{1},2)
                polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
                %place field center
                %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
            end
            
            %plot tuning specificity vector for all A trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
            
            %plot B (3) trial events
            for ll=1:size(idxMin{ss}{2},2)
                polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
                %place field center
                %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
            end
            
            %plot tuning specificity vector for all B trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
            
            hold off
            
            subplot(3,6,subplot_matrix(3,ss))
            imagesc(ROI_zooms{ss}{ROI})
            hold on;
            colormap(gca, 'gray')
            xticks([])
            yticks([])
            b = bwboundaries(ROI_outlines{ss}{ROI},'noholes');
            plot(b{1}(:,2),b{1}(:,1),'r')
            hold off
        end
    end
     mv_frame(ii) = getframe(gcf);
    pause(0.05);
    clf;
end
end

%% Plot raster , spiral, ROI FOV across  - 2 session comparison

%find day1 and da

%order of plots
subplot_matrix = 1:6;
subplot_matrix = reshape(subplot_matrix,2,3)';

figure('Position',[2230 30 780 960]);
for ii=1:size(registered.multi.assign_cell{1, 6},1) %with nans where no match
    %size(registered.multi.assigned_all,1) - match on each ses %options.idx_show%
    ses_nb = 1;
    %for each session
    for ss=[1,6]
        
        %ROI = registered.multi.assigned_all(ii,ss);
        %ROI = registered.multi.assigned(ii,ss);
        ROI = registered.multi.assign_cell{1,6}(ii,ses_nb);
        %skip of nan value
        if ~isnan(ROI)
            subplot(3,2,subplot_matrix(1,ses_nb))
            imagesc(session_vars{ss}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
            hold on;
            title(num2str(ROI));
            ylabel('Lap #');
            xlabel('Spatial bin');
            caxis([0 2])
            colormap(gca,'jet');
            hold off;
            
            %spiral plot early in learning
            subplot(3,2,subplot_matrix(2,ses_nb))
            polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
            hold on
            
            %plot A (2) trial events
            for ll=1:size(idxMin{ss}{1},2)
                polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
                %place field center
                %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
            end
            
            %plot tuning specificity vector for all A trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
            
            %plot B (3) trial events
            for ll=1:size(idxMin{ss}{2},2)
                polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
                %place field center
                %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
            end
            
            %plot tuning specificity vector for all B trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
            
            hold off
            
            subplot(3,2,subplot_matrix(3,ses_nb))
            imagesc(ROI_zooms{ss}{ROI})
            hold on;
            colormap(gca, 'gray')
            xticks([])
            yticks([])
            b = bwboundaries(ROI_outlines{ss}{ROI},'noholes');
            plot(b{1}(:,2),b{1}(:,1),'r')
            hold off
        end
        %update session number
         ses_nb = ses_nb +1;
    end
    pause(0.2)
    clf;
    disp(ii)
end


%% Write video of matching neurons across sessions

v = VideoWriter('I56_RTLS_all_matching_ROI.avi','Motion JPEG AVI');
open(v);
for ff=1:size(mv_frame,2)
    writeVideo(v,mv_frame(ff));
end
disp('done writing video');
close(v);

end

