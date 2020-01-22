function [outputArg1,outputArg2] = plot_raster_spiral_STC_ROI_multi_ses(plot_raster_vars,session_vars,registered,cat_registered_cell,...
    ROI_zooms,ROI_outlines,options)


%% Load variables

idxMin = plot_raster_vars.idxMin;
r_scaled = plot_raster_vars.r_scaled;
posVectorApprox = plot_raster_vars.posVectorApprox;
x = plot_raster_vars.x;

selectTrial = options.selectTrial;
sessionSelect = options.sessionSelect;

%% Plot raster , spiral, ROI FOV across  - 2 session comparison

%order of plots
%subplot_matrix = 1:6;
%subplot_matrix = [1 2; 3 4;5 6];
%subplot_matrix = reshape(subplot_matrix,2,3)';
match_mat = registered.multi.assigned_filtered;

matrix_order = [1:7;8:14;15:21];

figure('Position',[560 327 1114 517]);
%for each matching sets of ROIs
for ii=1:size(match_mat,1) %with nans where no match

    %for each session
    for ss=sessionSelect
        
        ROI = match_mat(ii,ss);
        %skip of nan value
        if ~isnan(ROI)
            subplot(3,size(sessionSelect,2),matrix_order(1,ss))
            
            imagesc(ROI_zooms{ss}{match_mat(ii,ss)})
            %else %do nothing or make black in future
            %set background axis color to black
            set(gca,'color',0*[1 1 1]);
            %end
            hold on;
            axis square
            
            colormap('gray')
            xticks([])
            yticks([])
            
            %ROI outlines
            b= bwboundaries(ROI_outlines{ss}{match_mat(ii,ss)},'noholes');
            plot(b{1}(:,2),b{1}(:,1),'r')
            %else
            %end
            hold off
            
            
            %spiral plot early in learning
            subplot(3,size(sessionSelect,2),matrix_order(2,ss))
            polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
            hold on
            %plot remapping label and index number of the neuron
            title([cat_registered_cell{ii,ss}, '\newline', num2str(match_mat(ii,ss))]);
            %title([cat_registered_cell{ii,ss}, '\newline','SCE all: ',num2str(multi_ses_SCE_data.SCE_all_ROI_engage(ii,ss))]);
            %plot A (2) trial events
            for ll=1:size(idxMin{ss}{1},2)
                polarscatter(angle(posVectorApprox{ss}{1}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
                %place field center
                %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
            end
            
            %plot tuning specificity vector for all A trials (unfiltered)
            polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(1)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
            
            %plot B (3) trial events
            for ll=1:size(idxMin{ss}{2},2)
                polarscatter(angle(posVectorApprox{ss}{2}{ll}{ROI}),r_scaled{ss}(idxMin{ss}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
                %place field center
                %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
            end
            
            %plot tuning specificity vector for all B trials
            polarplot([0+0i,15*session_vars{ss}.Place_cell{selectTrial(2)}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
            
            hold off

        end

    end
    pause()
    clf;
    disp(ii)
end

%% Try plot of ROI zooms/outlines
figure


for ii = 1:7
    subplot(1,7,ii)
%if no match
%if ~isnan(match_mat(ROI(rr),ss))
    imagesc(ROI_zooms{ii}{match_mat(1,ii)})
%else %do nothing or make black in future
    %set background axis color to black
    set(gca,'color',0*[1 1 1]);
%end
hold on;

colormap('gray')
xticks([])
yticks([])

%if empty, do not draw boundaries
%if ~isnan(match_list_sorted(ROI(rr),ss))
    b= bwboundaries(ROI_outlines{ii}{match_mat(1,ii)},'noholes');
    plot(b{1}(:,2),b{1}(:,1),'r')
%else
%end
hold off

end

end
