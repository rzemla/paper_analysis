function [outputArg1,outputArg2] = plot_raster_spiral_multi_ses_label_check(plot_raster_vars,session_vars,registered,cat_registered_cell,options)


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

figure('Position',[2230 30 780 960]);
for ii=1:size(match_mat,1) %with nans where no match

    %for each session
    for ss=sessionSelect
        
        ROI = match_mat(ii,ss);
        %skip of nan value
        if ~isnan(ROI)
            
            %spiral plot early in learning
            subplot(1,size(sessionSelect,2),ss)
            polarplot(x{ss},r_scaled{ss},'k','Linewidth',1.5)
            hold on
            title(cat_registered_cell{ii,ss});
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




end
