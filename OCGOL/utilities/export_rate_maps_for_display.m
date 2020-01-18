function [rate_maps_display] = export_rate_maps_for_display(session_vars)
%export rate maps for display - 
%for global near vs far examples
%for common vs. global 

%define Gaussian smoothing kernel
options.sigma_filter = 3;
gaussFilter = define_Gaussian_kernel(options);

%number of all ROIs
nb_ROI = size(session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm,2);

%expand to do circular convolution
non_norm_STC.A = session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm';
non_norm_STC.B = session_vars{1}.Place_cell{2}.Spatial_tuning_curve_no_norm';

%expand the edges
non_norm_STC_exp.A = [non_norm_STC.A(:,51:100),non_norm_STC.A,non_norm_STC.A(:,1:50)]';
non_norm_STC_exp.B = [non_norm_STC.B(:,51:100),non_norm_STC.B,non_norm_STC.B(:,1:50)]';

%for all ROIs
for rr=1:nb_ROI
    %smooth all A trial sessions
    A_smooth(:,rr) = conv(session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm(:,rr),gaussFilter, 'same');
    %circular convolution
    A_smooth_exp(:,rr) = conv(non_norm_STC_exp.A(:,rr),gaussFilter, 'same');
    
    %smooth all B trial sessions
    B_smooth(:,rr) = conv(session_vars{1}.Place_cell{2}.Spatial_tuning_curve_no_norm(:,rr),gaussFilter, 'same');
    %circular convolution
    B_smooth_exp(:,rr) = conv(non_norm_STC_exp.B(:,rr),gaussFilter, 'same');
    
end

%cicularly smoothed
A_smooth_circ = A_smooth_exp(51:150,:);
B_smooth_circ = B_smooth_exp(51:150,:);

%ROI - randi
%% Display random ROI
if 0
    figure
    hold on
    subplot(2,1,1)
    hold on
    plot(A_smooth_circ(:,1),'b')
    plot(B_smooth_circ(:,1),'r')
    
    subplot(2,1,2)
    hold on
    plot(A_smooth(:,1),'b')
    plot(B_smooth(:,1),'r')
end


%% Normalize to each other

% normalize from 0-1
smooth_maps_combined = [A_smooth' B_smooth'];

min_val_STC = min(smooth_maps_combined,[],2);
max_val_STC = max(smooth_maps_combined,[],2);

%normalize A/B matrix
A_smooth_norm = (A_smooth - min_val_STC')./(max_val_STC - min_val_STC)';
B_smooth_norm = (B_smooth - min_val_STC')./(max_val_STC - min_val_STC)';

% normalize from 0-1 (circularly
smooth_maps_combined_circ = [A_smooth_circ' B_smooth_circ'];

min_val_STC_circ = min(smooth_maps_combined_circ,[],2);
max_val_STC_circ = max(smooth_maps_combined_circ,[],2);

%normalize A/B matrix
A_smooth_norm_circ = (A_smooth_circ - min_val_STC_circ')./(max_val_STC_circ - min_val_STC_circ)';
B_smooth_norm_circ = (B_smooth_circ - min_val_STC_circ')./(max_val_STC_circ - min_val_STC_circ)';

%% Place into struct for export

%circularly smoothed and normalized to across A/B trials
rate_maps_display.A_smooth_norm_circ = A_smooth_norm_circ;
rate_maps_display.B_smooth_norm_circ = B_smooth_norm_circ;

%circularly smoothed and not normalized across A/B trials
rate_maps_display.A_smooth_circ = A_smooth_circ;
rate_maps_display.B_smooth_circ = B_smooth_circ;

%non normalized spatial tuning curves
rate_maps_display.non_norm_STC = non_norm_STC;

%% Display normalized values of circ smoothed ROIs
if 0
    figure
    hold on
    subplot(2,1,1)
    hold on
    plot(A_smooth_norm_circ(:,1),'b')
    plot(B_smooth_norm_circ(:,1),'r')
    
    subplot(2,1,2)
    hold on
    plot(A_smooth(:,1),'b')
    plot(B_smooth(:,1),'r')
end


end

