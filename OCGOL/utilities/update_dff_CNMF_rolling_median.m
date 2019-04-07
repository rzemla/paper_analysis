function [updated] = update_dff_CNMF_rolling_median(onset_offset, F_vars, Imaging, options)
%update dF/F with each iteration when events get detected

%% Select which inputs to use for dF/F recalculation (Rolling median)
if options.restrict == true
    F = F_vars.restricted.F;
    %F0_background = F_vars.restricted.F0_background;
elseif options.restrict == false
    F = F_vars.F;
    %F0_background = F_vars.F0_background;
end

%exponential smoothing constant (tau0?)
alpha = 1 - exp(-(0.033)/0.2); 

%update the percentile (for rolling median F0 baseline calculation)
options.prctile = 50;
%update the window (window over which to get the median) ~60s
options.window = 450; %450 ~15 seconds
%how many frames to shift as baseline is calculated
options.shift = 1;

%input F_df_exp (this is used when refinement is not done on F0 baseline
F_dff_exp_restricted = Imaging.trace_restricted;

%% Mask F where events were detected - used for both Jia/Danielson and rolling median portions of the code below

%mask trace with NaNs where event detected
%mask F to calculate F0_baseline
F_masked = F';

for rr = 1:size(F_masked,2)
    for ee=1:size(onset_offset{rr},1)
        %was not but changed to 0
        F_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = nan;
    end
end

%% Using Rolling median (manual) to recalculate the with masked events dF/F

%get matrix of nans for F_sm_masked
F_masked_nan = isnan(F_masked);

%make an indexing reference to allow for reconstruction of smoothed
%smoothed F mask
frame_idx = 1:size(F,2);

size(F)
size(F_masked_nan)

%make cell for each ROI F with event nans removed
for rr=1:size(F_masked_nan,2)
    F_nonan{rr} = F(rr,~F_masked_nan(:,rr));
    F_nonan_idx{rr} = frame_idx(~F_masked_nan(:,rr));
end

%feed with nans prcfilt can work with nan inputs and window is long enough
F0 = prctfilt(F_masked',options.prctile, options.window,options.shift,0);

%fill in nans - don't need to b/c long tau - movmin can work with nans in
%place b/c tau is long

%with added background for the component
F_df =(F - F0)./(F0);

%time by ROI inputs
F_df_exp = filter(alpha, [1 alpha-1], F_df');

%% Calculate sigma with events removed using Danielson/Jia dF/F for each ROI

%mask F_df_exp
F_df_exp_masked = F_df_exp;

for rr=1:size(F_df_exp_masked,2)
    for ee=1:size(onset_offset{rr},1)
        F_df_exp_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = nan;
    end
end

%calculate mean and sigma of the event masked dF/F
updated_mean_df = nanmean(F_df_exp_masked);
updated_std_df = nanstd(F_df_exp_masked);

%% Plot

ii=1;

figure
subplot(4,1,1)
hold on
title('F, F masked, F0')
p1 = plot(F(ii,:),'k');
p2 = plot(F_masked(:,ii),'r');
%updated baseline
p3 = plot(F0(ii,:),'b');
legend([p1 p2 p3],'F', 'F masked', 'new F0')

subplot(4,1,2)
hold on
title('Updated dF/F with new mean and std')
plot(F_df_exp(:,ii),'k')
%updated mean
refline(0,updated_mean_df(:,ii))
%updated mean+std
refline(0,updated_std_df(:,ii)+updated_mean_df(:,ii))
hold off


%% Plot changes before and after 
figure

ii=1;

subplot(3,1,1)
hold on
title('Original vs. Updated F0')
%raw F
p1 = plot(F_vars.restricted.F(ii,:),'k');
%original baseline
p2 = plot(F_vars.restricted.F0(ii,:),'r');
%new baseline
p3 = plot(F0(ii,:),'m');

legend([p1 p2 p3],'F', 'Original F0', 'Updated F0')
hold off

subplot(3,1,2)
hold on
title('Updated baseline')
%raw F
plot(F(ii,:),'k')
%redone baseline
plot(F0(ii,:),'r');
hold off

subplot(3,1,3)
hold on
%old dF/F exp
title('Incoming dF/F vs. updated dF/F')
p4 = plot(F_vars.restricted.F_df_exp(ii,:),'k');
%redone baseline
p5 = plot(F_df_exp(:,ii),'r');

legend([p4 p5],'Previous dF/F exp','Updated dF/F exp')
hold off

%% Recalculate rolling median dFF
% 
% %recalculate baseline F0
% F0_baseline_F_masked = prctfilt(F_masked',options.df_prctile,options.df_window,[],0); 
% 
% %used fixed background compomnent from background import 
% F_dff_baseline_masked = (F - F0_baseline_F_masked)./(F0_baseline_F_masked + F0_background);
% 
% %exponentially filter
% F_dff_exp_baseline_masked = filter(alpha, [1 alpha-1], F_dff_baseline_masked');
% 
% F_dff_exp = F_dff_exp_baseline_masked;
% 
% %% Calculate mean and sigma of the event masked dF/F for rolling median dFF
% 
% %mask F_df_exp
% F_dff_exp_masked = F_dff_exp;
% 
% for rr=1:size(F_dff_exp_masked,2)
%     for ee=1:size(onset_offset{rr},1)
%         F_dff_exp_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = NaN;
%     end
% end
% 
% updated_mean_dff = nanmean(F_dff_exp_masked);
% updated_std_dff = nanstd(F_dff_exp_masked);

%% Make structure
%output is transposed: time x ROI
updated.F_df_exp = F_df_exp; %Jia/rolling median (manual) updated
updated.F0 = F0; %updated raw F baseline
updated.mean_df = updated_mean_df;
updated.std_df = updated_std_df;



end

