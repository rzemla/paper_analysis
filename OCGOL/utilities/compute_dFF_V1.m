function [F_vars] = compute_dFF_V1(CNMF_output,directory_name)

%% Place struct vars into function workspace

%CNMF parameters
options = CNMF_output.options; 

%final spatial and temporal components 
A_keep = CNMF_output.A_keep;
C_full = CNMF_output.C_full;

%field of view
FOV = CNMF_output.FOV;

%background related
b = CNMF_output.b;
f_full = CNMF_output.f_full;
F_dark = CNMF_output.F_dark;

%length of frames series
T = CNMF_output.T;
%residual
R_full = CNMF_output.R_full;


%% Calculate the dF/F from CNMF with median percentile and 15s baseline window

%update the percentile
options.df_prctile = 50;

%update the window
options.df_window = 450; %450 ~15 seconds %1800 - 60s

%windowing shift every
%options.shift = 1;

tic;
[F_dff,F0,F,Fd,F0_background, F0_baseline] = detrend_df_f_output_mod_V2(A_keep,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
toc;
%calculate the background for all of the componenents (used in Jia
%calculation below as well)
%F0_background = prctfilt((A_keep'*[b,ones(prod(FOV),1)])*[f_full;-double(F_dark)*ones(1,T)],options.df_prctile,options.df_window,[],0);

%this if the formula for the dF/F
F_dff_external = (F-F0_baseline)./(F0_baseline + F0_background);

%% Plot dF/F transformations


%% Exponentially filter (smooth) the CNMF dF/F

%smoothing coefficient for expoenential filter
%%https://en.wikipedia.org/wiki/Exponential_smoothing#Comparison_with_moving_average
%i.e. alpha =1 - e(-deltaT/tau); tau is time constant and deltaT is the sampling
%time interval; 
%sampling time much faster compared to time constant:
%alpha = deltaT/tau
%tau of ~0.2s based on Jia 2011/Danielson 2016b

%time constant for smoothing
tau0 = 0.2;

%alpha parameter for matlab filter function
alpha = 1 -exp(-(0.033)/tau0); 

%smoothed dF/F
F_dff_exp = filter(alpha, [1 alpha-1], F_dff');

%transposed smoothed dF/F (ROI x frames)
F_dff_exp = F_dff_exp';



%% Plot check
figure
for ii =1:20
    ROI =ii;
    %disp(ROI);
    
    subplot(4,1,1)
    hold on
    title(['ROI: ', num2str(ROI)])
    ylabel('Raw Fluorescence');
    p1 = plot(F(ROI,:),'k');
    p2 = plot(F0_baseline(ROI,:),'r');
    p3 = plot(F0_background(ROI,:),'y');
    p4 = plot(F0(ROI,:),'b');
    legend([p1 p2 p3 p4],{'Raw F', 'F0_baseline', 'F0 background', 'F0'})
    hold off
    
    subplot(4,1,2)
    hold on
    title('F - F0\_baseline and F0');
    ylabel('Raw Fluorescence');
    p5 = plot(Fd(ROI,:),'k');
    p6 = plot(F0(ROI,:),'r');
    legend([p5 p6], {'F-F0\_baseline','F0'})
    hold off
    
    %no de-noised with exp filter
    subplot(4,1,3)
    hold on
    ylabel('dF/F');
    title('dF/F without smoothing - background included');
    plot(F_dff(ROI,:),'k');
    hold off
    
    %denoised with exp filter
    subplot(4,1,4)
    hold on
    ylim([-0.2 2.5]);
    ylabel('dF/F')
    title('Exponentially smoothed/de-noised');
    plot(F_dff_exp(ROI,:),'k');
    hold off
    
    %pause (0.05);
    clf
end

%% Create output struct
%raw F
F_vars.F = F;
%dF/F (non smoothed/de-noised)
F_vars.F_dff = F_dff;
%F0 = F0_background + F0_baseline 
F_vars.F0 = F0;
%F0_baseline
F_vars.F0_baseline = F0_baseline;
%F0_background
F_vars.F0_background = F0_background;
%Fd = F - F0_baseline;
F_vars.Fd = Fd;
%dF/F (exp. smoothed/de-noised)
F_vars.F_dff_exp = F_dff_exp;

%% Save struct as workspace

%save the read_input directory of experiment directory
save(fullfile(directory_name,'read_inputs','F_vars.mat'),'F_vars');


end

