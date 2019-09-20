function [outputArg1,outputArg2] = visualize_neuron_characteristics(plot_raster_vars,norm_events,registered,session_vars,task_selective_ROIs,cat_registered_cell,select_fields,options)

%% Define session/trial inputs
selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Spiral plot params
idxMin = plot_raster_vars.idxMin;
r_scaled = plot_raster_vars.r_scaled;
posVectorApprox = plot_raster_vars.posVectorApprox;
x = plot_raster_vars.x;

%% Load in position normalized events
event_norm_time = norm_events.event_norm_time;
event_norm_pos_run = norm_events.event_norm_pos_run;
event_lap_idx = norm_events.event_lap_idx;


%% Convert place field edges to angles and complex vector for polar plotting


for ss =sessionSelect
    for rr=1:size(session_vars{ss}.Place_cell{4}.placeField.edge,2)
        %unit vector
        session_vars{ss}.Place_cell{4}.placeField.edge_vec{rr} = exp(i*deg2rad((session_vars{ss}.Place_cell{4}.placeField.edge{rr}/100)*360));
        session_vars{ss}.Place_cell{5}.placeField.edge_vec{rr} = exp(i*deg2rad((session_vars{ss}.Place_cell{5}.placeField.edge{rr}/100)*360));
    end
end

%% Plot spiral and event map, place field edges and center

ss=4;

selectROI = task_selective_ROIs{ss}.A.idx;

figure('Position',[1961 60 938 902]);
for ii=1:size(selectROI,2) %with nans where no match
    %for each session

        ROI = selectROI(ii);
            
            %spiral plot early in learning
            subplot(2,2,1)
            polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
            hold on
            title(num2str(ROI))
            %title(cat_registered_cell{ii,ss});
            %plot A (2) trial events
            for ll=1:size(idxMin{ss}{1},2)
                polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
                %place field center
                %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
            end
            
            %plot tuning specificity vector for all A trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
            
            for pp=find(select_fields{ss}{4}{ROI} ==1)
            %plot place field edges as vectors
            polarplot([0+0i,15*session_vars{ss}.Place_cell{4}.placeField.edge_vec{ROI}(pp,1)],'b--','LineWidth',2)
            polarplot([0+0i,15*session_vars{ss}.Place_cell{4}.placeField.edge_vec{ROI}(pp,2)],'b--','LineWidth',2)
            end

            %plot B (3) trial events
            for ll=1:size(idxMin{ss}{2},2)
                polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
                %place field center
                %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
            end
            
            %plot tuning specificity vector for all B trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{1, 5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
            

            for pp=find(select_fields{ss}{5}{ROI} ==1)
            %plot place field edges as vectors
            polarplot([0+0i,15*session_vars{ss}.Place_cell{5}.placeField.edge_vec{ROI}(pp,1)],'r--','LineWidth',2)
            polarplot([0+0i,15*session_vars{ss}.Place_cell{5}.placeField.edge_vec{ROI}(pp,2)],'r--','LineWidth',2)
            end
            hold off

            %event map
            subplot(2,2,2)
            hold on
            plot(session_vars{ss}.Place_cell{4}.Spatial_Info.event_map{8}(:,ROI),'b')
            plot(session_vars{ss}.Place_cell{5}.Spatial_Info.event_map{8}(:,ROI),'r')
            hold off
            subplot(2,2,3)
            hold on
            title('Smoothed rate map')
            %plot(session_vars{1}.Place_cell{4}.Spatial_Info.rate_map{8}(:,ROI),'b')
            %plot(session_vars{1}.Place_cell{5}.Spatial_Info.rate_map{8}(:,ROI),'r')
            plot(session_vars{ss}.Place_cell{4}.Spatial_Info.rate_map_smooth{8}(:,ROI),'b') 
            plot(session_vars{ss}.Place_cell{5}.Spatial_Info.rate_map_smooth{8}(:,ROI),'r') 
            
            %A
            %plot place field edges
            for pp=find(select_fields{ss}{4}{ROI} ==1)
                plot([session_vars{ss}.Place_cell{4}.placeField.edge{ROI}(pp,1),session_vars{ss}.Place_cell{4}.placeField.edge{ROI}(pp,1)],[0 1],'LineStyle','-','Color',[0,0,1,0.5])
                plot([session_vars{ss}.Place_cell{4}.placeField.edge{ROI}(pp,2),session_vars{ss}.Place_cell{4}.placeField.edge{ROI}(pp,2)],[0 1],'LineStyle','-','Color',[0,0,1,0.5])
            end
            
        %B
            %plot place field edges
            for pp=find(select_fields{ss}{5}{ROI} ==1)
                plot([session_vars{ss}.Place_cell{5}.placeField.edge{ROI}(pp,1),session_vars{ss}.Place_cell{5}.placeField.edge{ROI}(pp,1)],[0 1],'LineStyle','--','Color',[1,0,0,0.5])
                plot([session_vars{ss}.Place_cell{5}.placeField.edge{ROI}(pp,2),session_vars{ss}.Place_cell{5}.placeField.edge{ROI}(pp,2)],[0 1],'LineStyle','--','Color',[1,0,0,0.5])
            end

            hold off
            subplot(2,2,4)
            title(['S.I.: ', num2str(session_vars{ss}.Place_cell{4}.Spatial_Info.Spatial_Info(8,ROI)),...
                        ' ',num2str(session_vars{ss}.Place_cell{5}.Spatial_Info.Spatial_Info(8,ROI)),'\newline',...
                        'p: ',num2str(session_vars{ss}.Place_cell{4}.Spatial_Info.ROI_pvalue(ROI)),' ',num2str(session_vars{ss}.Place_cell{5}.Spatial_Info.ROI_pvalue(ROI)),...
                        '\newline',...
                        'T.S.: ', num2str(session_vars{ss}.Place_cell{4}.Tuning_Specificity.tuning_specificity(ROI)),...
                                  ' ',num2str(session_vars{ss}.Place_cell{5}.Tuning_Specificity.tuning_specificity(ROI)),'\newline',...
                        'p: ', num2str(session_vars{ss}.Place_cell{4}.Tuning_Specificity.ROI_pvalue(ROI)),' ',...
                           num2str(session_vars{ss}.Place_cell{5}.Tuning_Specificity.ROI_pvalue(ROI))])
            
    pause()
    clf;
    disp(ii)
end


%% Calculate unweight circular variance (unline TS which is weighed by occupancy)



%% Plot event density vs.position

select_ROIs =task_selective_ROIs{1}.A.idx; 

for rr=1:size(select_ROIs,2)
    event_pos_temp = event_norm_pos_run{1}.A{select_ROIs(rr)}
    event_pos_temp_2 = event_norm_pos_run{1}.B{select_ROIs(rr)}
    
    event_pos_deg = circ_ang2rad(event_pos_temp*360)
    event_pos_deg_2 = circ_ang2rad(event_pos_temp_2*360)
    
    figure
    subplot(1,2,1)
    circ_plot(event_pos_deg,'hist',[],25,false)
    
    subplot(1,2,2)
    circ_plot(event_pos_deg_2,'hist',[],25,false)
pause
clf
end

end

