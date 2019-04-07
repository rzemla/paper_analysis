function [centroids] = centroid_diff_reward(session_vars,options)
%returns reward location difference for each session

%% Extract tuning spec structs and assign reward location
%tuning specificities struct for each session
for ii = 1:size(session_vars,2)
    ts{ii} = session_vars{ii}.Place_cell{1}.Tuning_Specificity;
    %cm
    rewardLoc(ii) = session_vars{ii}.Behavior.performance.reward_loc;
end
    

%% Get each tuning vector (complex and cartesian) for each ROI

%for each session
for ss=1:size(session_vars,2)
    
    %for each ROI in that session
    for ROI = 1:size(ts{ss}.tuning_vector,2)
        %get angle for that ROI from each session relative to 0
        %angle of each vector
        ang{ss}(ROI) = angle(ts{ss}.tuning_vector_specificity(ROI));
        %complex vector
        u{ss}(ROI) = ts{ss}.tuning_vector_specificity(ROI);
        %cartesian converted
        uCar{ss}{ROI} = [real(u{ss}(ROI)), imag(u{ss}(ROI))];
    end
    
end

%% Calculate difference from each centroid to 

%for each session
%largest angle difference should be pi
for ss=1:size(session_vars,2)
    [theta{ss}] = centroid_angle_diff(uCar{ss},rewardLoc(ss));
end

%input two complex vectors in neighboring cells
%get smallest angle between vectors
% theta(ROI) = atan2(uCar{1}(1)*uCar{2}(2) - uCar{1}(2)*uCar{2}(1), uCar{1}(1)*uCar{2}(1) + uCar{1}(2)*uCar{2}(2));
% theta(ROI) = abs(theta(ROI));

%using formula
% u = ts{1}.tuning_vector_specificity(ROI);
% v = ts{2}.tuning_vector_specificity(ROI);
% 
% %cartesian definition
% uCar = [real(u), imag(u)]; %[x1 , y1] v1
% vCar = [real(v), imag(v)]; %[x2 , y2] v2
%     
% %a = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
% %gives the angle in degrees between the vectors as measured in a counterclockwise 
% %direction from v1 to v2. If that angle would exceed 180 degrees, then the angle 
% %is measured in the clockwise direction but given a negative value. 
% %In other words, the output of 'atan2d' always ranges from -180 to +180 degrees.
% 
% %https://www.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
% theta = atan2(uCar(1)*vCar(2) - uCar(2)*vCar(1), uCar(1)*vCar(1) + uCar(2)*vCar(2));

%% ROIs that significantly TS tuned in each session

for ss=1:size(session_vars,2)
    theta_ts_sig{ss} = theta{ss}(ts{ss}.significant_ROI);
end


%% Save to structure

%difference between tuning vectors (absolute radians)
centroids.angleDiffReward = theta;
centroids.TS_sig.angleDiffReward = theta_ts_sig;
%centroids.TS_sig.significant_ROI = find;


%% Plot a histogram of the centroid difference between neurons with significantly tuned according to tuning specificity criterion

labels = {'Day 0 RF ', 'Day 1 GOL 1', 'Day 2 GOL 1', 'Day 3 GOL 1',...
    'Day 4 GOL 2','Day 5 GOL 2','Day 6 GOL 2'};
plot_order = [1:2:8, 2:2:8];
figure;
for ii=1:size(session_vars,2)
    subplot(4,2,plot_order(ii))
    hold on
    ax_h = axis;
    
    title(labels{ii});
    histogram(centroids.TS_sig.angleDiffReward{ii},20,'Normalization','probability');
    xticks([0.7854, 1.5708, 2.3562, 3.1416])
    xticklabels({'\pi/4','\pi/2','3\pi/4','\pi'})
    axis(ax_h, 'tight')
    ylim([0 0.2]);
    ylabel('Normalized count');
    hold off
end

%plot threshold line
%plot([options.centroidThres, options.centroidThres],[0 0.5],'r')

%% Plot angles and visualize difference for each ROI

if 0

    f = figure;
    
    for ROI = 1:size(Place_cell{1}.Tuned_ROI_mask,2)
        
        compass(ts{1}.tuning_vector_specificity(ROI),'b');
        hold on
        %compass(ts{2}.tuning_vector_specificity(ROI),'r');
        
        %plot unit difference vector as unit
        compass(exp(1i*theta{1}(107)),'g');
        
        pause;
        hold off
        
    end
end

end

