function [outputArg1,outputArg2] = cent_diff_plotter(return_angle_diff)



%% Setup data for plotting - pooled data
%ts
%A
st_learn.pooled.ts.A = [1:7; return_angle_diff.ts.st_learn.pooled.mean.A(1:7); return_angle_diff.ts.st_learn.pooled.sem.A(1:7)];
st_recall.pooled.ts.A = [[1:3,6:9]; return_angle_diff.ts.st_recall.pooled.mean.A([1:3,6:9]); return_angle_diff.ts.st_recall.pooled.sem.A([1:3,6:9])];
%B
st_learn.pooled.ts.B = [1:7; return_angle_diff.ts.st_learn.pooled.mean.B(1:7); return_angle_diff.ts.st_learn.pooled.sem.B(1:7)];
st_recall.pooled.ts.B = [[1:3,6:9]; return_angle_diff.ts.st_recall.pooled.mean.B([1:3,6:9]); return_angle_diff.ts.st_recall.pooled.sem.B([1:3,6:9])];

%si
%A
st_learn.pooled.si.A = [1:7; return_angle_diff.si.st_learn.pooled.mean.A(1:7); return_angle_diff.si.st_learn.pooled.sem.A(1:7)];
st_recall.pooled.si.A = [[1:3,6:9]; return_angle_diff.si.st_recall.pooled.mean.A([1:3,6:9]); return_angle_diff.si.st_recall.pooled.sem.A([1:3,6:9])];
%B
st_learn.pooled.si.B = [1:7; return_angle_diff.si.st_learn.pooled.mean.B(1:7); return_angle_diff.si.st_learn.pooled.sem.B(1:7)];
st_recall.pooled.si.B = [[1:3,6:9]; return_angle_diff.si.st_recall.pooled.mean.B([1:3,6:9]); return_angle_diff.si.st_recall.pooled.sem.B([1:3,6:9])];

%% Setup data for plotting - by animal


%
%% Plotter

%A
input_data.A{1,1} = st_learn.pooled.ts.A;
input_data.A{1,2} = st_recall.pooled.ts.A;
%B
input_data.B{2,1} = st_learn.pooled.ts.B;
input_data.B{2,2} = st_recall.pooled.ts.B;
%A
input_data.A{3,1} = st_learn.pooled.si.A;
input_data.A{3,2} = st_recall.pooled.si.A;
%B
input_data.B{4,1} = st_learn.pooled.si.B;
input_data.B{4,2} = st_recall.pooled.si.B;

% input_data.AB{3} = TC_mean_sem.ABcorr_animal.si.st_learn.AB;
% input_data.AB{4} = TC_mean_sem.ABcorr_animal.si.st_recall.AB;

title_labels{1} = 'A Centroid diff pooled - TS';
title_labels{2} = 'B Centroid diff pooled - TS';

title_labels{3} = 'A Centroid diff pooled - SI';
title_labels{4} = 'B Centroid diff pooled - SI';

%rad to cm - ~196 = 2*pi 
cm2rad = (2*pi)./196;
%rad ticks that correspond to 0 25 50 cm
rad_ticks = [0 10 20 30 40].*cm2rad;

figure('Position', [2274 120 720 770])
for ii=1:4
    %learn vs raw PV
    subplot(2,2,ii)
    hold on
    axis square
    hold on
    xlabel('Relative day')
    ylim([0 1.4])
    yticks(rad_ticks)
    yticklabels({'0','10','20','30','40'});
    %if A (left side)
    if rem(ii,2) == 1
        xticks([1:9])
        xlim([0 10])
        xticklabels({'1','2','3','4','5','6','7','8','9'})
        %learn
    lA = plot_error_line(input_data.A{ii,1},'--',2,[65,105,225]./255);
    %if B (right side)
    rA = plot_error_line(input_data.A{ii,2},'-',2,[65,105,225]./255);
        
            %dashed 1 reference line
    %plot([0 8],[1 1],'--','Color',[0.5 0.5 0.5])
    else %recall
        xticks([1:9])
        xlim([0 10])
        xticklabels({'1','2','3','4','5','6','7','8','9'})
                %learn
    lB = plot_error_line(input_data.B{ii,1},'--',2,[220,20,60]./255);
    %recalll
    rB = plot_error_line(input_data.B{ii,2},'-',2,[220,20,60]./255);
            %dashed 1 reference line
    %plot([0 10],[1 1],'--','Color',[0.5 0.5 0.5])
    end
    
    %xtickangle(45)
    title(title_labels{ii})
    
    %plot correlation on left y axis

    %lB = plot_error_line(input_data.B{ii},'-',2,[220,20,60]/255);
    
    %plot correlation on right y axis
    
    set(gca,'FontSize',12)
    set(gca,'Linewidth',2)
    
    
    %legend([lA,lB],{'A','B'},'location','northeast')
    
end


end

