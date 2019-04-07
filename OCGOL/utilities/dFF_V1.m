function [F_vars] = dFF_V1(dff_vars,ROI)

%% Place struct vars into function workspace

options = dff_vars.options; 
A_keep = dff_vars.A_keep;
b = dff_vars.b;
FOV = dff_vars.FOV;
C_full = dff_vars.C_full;
f_full = dff_vars.f_full;
F_dark = dff_vars.F_dark;
T = dff_vars.T;
R_full = dff_vars.R_full;


%% Calculate the dF/F from CNMF with median percentile and 15s baseline window

%update the percentile
options.df_prctile = 50;

%update the window
options.df_window = 1800; %450 ~15 seconds %1800 - 60s

[F_dff,F0_d,F,Fd] = detrend_df_f_RZ_V1(A_keep,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);

%calculate the background for all of the componenents (used in Jia
%calculation below as well)
F0_background = prctfilt((A_keep'*[b,ones(prod(FOV),1)])*[f_full;-double(F_dark)*ones(1,T)],options.df_prctile,options.df_window,[],0);

%% Exponentially filter (smooth) the CNMF dF/F

%exponential filter
%alpha = 0.05; %window size of 39 frames ~1.2s

%smoothing coefficient for expoenential filter
%%https://en.wikipedia.org/wiki/Exponential_smoothing#Comparison_with_moving_average
%i.e. alpha =1 - e(-deltaT/tau); tau is time constant and deltaT is the sampling
%time interval; 
%sampling time much faster compared to time constant:
%alpha = deltaT/tau
%tau of ~0.2s based on Jia 2011/Danielson 2016b
alpha = 1 -exp(-(0.033)/0.2); 

F_dff_exp = filter(alpha, [1 alpha-1], F_dff');
F_dff_exp = F_dff_exp';

%OUTDATED - esssentially the same as given with filter function
%smoothing on all ROIS
%F_dff_exp2 = smoothts(F_dff(ROI,:),'e',6);


%% dF/F calculation based in Jia 2011/ Danielson 2016b
%F is the raw fluorescence from CNMF; Fd is the detrended fluorescence

%doesn't work well with the parameters established on CNMF F or Fd
%Danielson 2016/Jia 2011 dF/F calculation
%calculate in frame space
%calculate baseline F0
startTau = 200; %average min for the first frames before tau
tau2 =1820;   %1818 frames ==> ~60 seconds 
tau1 = 91; %frames) 90 frames ==> ~3 seconds

%smooth the raw F from CNMF for use in baseline calculation
%uses the points around the sample point to get average
%different handling of the endpoints compared to filter function
%deals with edges automatically and no need to manually offset the delay
%from filter function
%F_sm_ROI = smooth(F(ROI,:),tau1,'moving')';


%all ROIs
tic;
for rr=1:size(F,1)
    F_sm(rr,:) = smooth(F(rr,:)',tau1,'moving');
end
toc;

%same approach but using a the filter function in matlab
% filter_offset = (tau1 - 1)/2;
% F_sm_f = filter((1/tau1)*ones(1,tau1),1,F(ROI,:));
%build in matlab function that does the same thing as above

%F0_start = movmin(F_sm,startTau, 2);
%alternative way of dealing with start point
F0_start = movmin(F_sm,[0, tau2], 2);

%one sample offset due to way movmin calculates the minimum
F0_end = movmin(F_sm,[tau2, 0], 2);

%account for the shift from movmin (includes the time point sampled)
F0_end_offset = F0_end;
F0_end_offset(:,round(tau2/2):end) = F0_end(:,round(tau2/2)-1:end-1);

F0 = [F0_start(:,1:tau2), F0_end_offset(:,tau2+1:end)];

%F0 = F0_end_offset;
%with added background for the component
F_df =(F - F0)./(F0+F0_background);

%F_df =(F - F0)./(F0);

%exponentially smooth it with tau = 0.2s
alpha = 1 -exp(-(0.033)/0.2); 

F_df_exp = filter(alpha, [1 alpha-1], F_df');
F_df_exp = F_df_exp';

%% Plot a sample dF/F trace from each type of dF/F calculation

figure;

%plot Jia dF/F
subplot(3,1,1)
hold on
p1 = plot(F(ROI,:));
p2 = plot(F_sm(ROI,:));
p3 = plot(F0(ROI,:));
p3a = plot(F0_background(ROI,:) + F0(ROI,:),'k');
hold off
legend([p1,p2,p3], 'F', 'F sm', 'F0');

subplot(3,1,2)
hold on
p4 = plot(F(ROI,:));
p5 = plot(F0_background(ROI,:),'r') ;
p6 = plot(F(ROI,:) - Fd(ROI,:),'c');
p7 = plot(F0_d(ROI,:),'m');

hold off
legend([p4,p5,p6,p7], 'F', 'F0 background', 'F0 baseline', 'F0');

%plot exponenetially smoothed dF/F's
subplot(3,1,3)
hold on;
p8 = plot(F_df_exp(ROI,:));
p9 = plot(F_dff_exp(ROI,:));
hold off
legend([p8,p9], 'Jia/Danielson dFF', 'Rolling median dFF');

%both danielson and CNMF dF/F are nearly extractly the same when using Fd
%(just off by serveral orders of magnitude)
% %as input to Danielson/Jia dF/F calculation and 98% when using F as input

%corr score mostly different due way starting timepoint is treated
corr(F_df_exp(ROI,:)',F_dff_exp(ROI,:)')



%% Create output struct
F_vars.F = F;
F_vars.F_dff = F_dff;
F_vars.F0= F0;
F_vars.F0_d = F0_d;
F_vars.F0_background =F0_background;
F_vars.Fd = Fd;
F_vars.F_df_exp = F_df_exp;
F_vars.F_dff_exp = F_dff_exp;


end

