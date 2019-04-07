function [updated] = update_dff_revised(onset_offset, F_vars, options)

%update dF/F with each iteration when events get detected

%select which inputs to use for dF/F recalculation (Rolling median)
if options.restrict == true
    F = F_vars.restricted.F;
    F0_background = F_vars.restricted.F0_background;
elseif options.restrict == false
    F = F_vars.F;
    F0_background = F_vars.F0_background;
end

%need a condition if restricted condition is met here

%exponential smoothing constant
alpha = 1 - exp(-(0.033)/0.2); 

%update the percentile (for rolling median F0 baseline calculation)
options.df_prctile = 50;

%update the window (window over which to get the median) ~60s
options.df_window = 450; %450 ~15 seconds

%% Mask F where events were detected - used for both Jia/Danielson and rolling median portions of the code below

%mask trace with NaNs where event detected
%mask F to calculate F0_baseline
F_masked = F';

for rr = 1:size(F_masked,2)
    for ee=1:size(onset_offset{rr},1)
        F_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = NaN;
    end
end

%% Using Danielson/Jia method recalculate the with masked events dF/F

%define smoothing and background constants
startTau = 200; %average min for the first frames before tau
tau2 =1820;   %1818 frames ==> ~60 seconds 
tau1 = 91; %frames) 90 frames ==> ~3 seconds

%F_sm_masked = smoothdata(F_masked,'movmean',tau1,'includenan');

%smooth Fsm on event masked F 
%all ROIs
tic;
for rr=1:size(F,1)
    
    F_sm_masked(:,rr) = smooth(F_masked(:,rr),tau1,'moving')';
end
toc;

%same approach but using a the filter function in matlab
% filter_offset = (tau1 - 1)/2;
% F_sm_f = filter((1/tau1)*ones(1,tau1),1,F(ROI,:));
%build in matlab function that does the same thing as above

%F0_start = movmin(F_sm_masked',startTau, 2);

%alternative way of dealing with start point
F0_start = movmin(F_sm_masked',[0, tau2], 2);

%one sample offset due to way movmin calculates the minimum
F0_end = movmin(F_sm_masked',[tau2, 0], 2);

%account for the shift due to movmin including the sample point
F0_end_offset = F0_end;
F0_end_offset(:,round(tau2/2):end) = F0_end(:,round(tau2/2)-1:end-1);

%combined baselines from the first tau2 frames with not enough time to
%calculate baseline with post-tau2 frames where there is enough time
F0 = [F0_start(:,1:tau2), F0_end_offset(:,tau2+1:end)];

%F0 = F0_end_offset;
%with added background for the component
F_df =(F - F0)./(F0+F0_background);

F_df_exp = filter(alpha, [1 alpha-1], F_df');

%% Calculate sigma with events removed using Danielson/Jia dF/F for each ROI

%mask F_df_exp
F_df_exp_masked = F_df_exp;

for rr=1:size(F_df_exp_masked,2)
    for ee=1:size(onset_offset{rr},1)
        F_df_exp_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = 0;
    end
end

%calculate mean and sigma of the event masked dF/F
updated_mean_df = nanmean(F_df_exp_masked);
updated_std_df = nanstd(F_df_exp_masked);


%% Recalculate rolling median dFF

%recalculate baseline F0
F0_baseline_F_masked = prctfilt(F_masked',options.df_prctile,options.df_window,[],0); 

%used fixed background compomnent from background import 
F_dff_baseline_masked = (F - F0_baseline_F_masked)./(F0_baseline_F_masked + F0_background);

%exponentially filter
F_dff_exp_baseline_masked = filter(alpha, [1 alpha-1], F_dff_baseline_masked');

F_dff_exp = F_dff_exp_baseline_masked;

%% Calculate mean and sigma of the event masked dF/F for rolling median dFF

%mask F_df_exp
F_dff_exp_masked = F_dff_exp;

for rr=1:size(F_dff_exp_masked,2)
    for ee=1:size(onset_offset{rr},1)
        F_dff_exp_masked(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr) = NaN;
    end
end

updated_mean_dff = nanmean(F_dff_exp_masked);
updated_std_dff = nanstd(F_dff_exp_masked);

%% Make structure

updated.F_df_exp = F_df_exp; %Jia updated
updated.mean_df = updated_mean_df;
updated.std_df = updated_std_df;

updated.F_dff_exp = F_dff_exp; %rolling median updated
updated.mean_dff = updated_mean_dff;
updated.std_dff = updated_std_dff;

%%

%%%%%%%SCRAP CODE %%%%%%%%%%%%%%%%%%%%%%%%%
% %time-dependent baseline calculation
% for ii=1:length(F_sm_masked)
%     if ii==1
%         F0_masked(ii) = min(F_sm_masked((ii+1):(tau2+1)));
%     elseif ii<(tau2+1)
%         F0_masked(ii) =  min([F_sm_masked(2:(ii-1)),F_sm_masked((ii+1):(tau2-ii-1))]);
%     else
%         F0_masked(ii) = min(F_sm_masked((ii-tau2):(ii-1)));
%     end
% end
% 
% %with added background for the component
% F_df_masked =(F(ROI,:) - F0_masked)./(F0_masked+F0_background(ROI,:));
% 
% %exponentially smooth it with tau = 0.2s
% 
% 
% F_df_exp_masked = filter(alpha, [1 alpha-1], F_df_masked);

end

