function [F_dff,F0,F,Fd,F0_background, F0_baseline] = detrend_df_f_output_mod_V2(A,b,C,f,YrA,options)

% V1 RZ - outputs raw fluorescence (F) and detrended F (Fd)

% detrend and extract DF/F values from processed fluorescence traces
% INPUTS
% A:        matrix of spatial components (2-d matrix, sparse)
% b:        matrix of spatial background (2-d matrix)
% C:        matrix of temporal components (2-d matrix)
% f:        matrix of temporal background (2-d matrix)
% YrA:      matrix of filtered residuals (2-d matrix, optional)
% options:  options structure used for specifying method for determining DF
%           default method is the median of the trace. By changing
%           options.df_prctile an arbitray percentile can be used (between 0 and 100).
%           a moving window can also be established by specifying options.df_window
%
% OUTPUTS
% F_dff:    detrended fluorescence in DF/F 
% F0:       baseline fluorescence for each trace

defoptions = CNMFSetParms;
if nargin < 5 || isempty(options)
    options = defoptions;
    disp('Options assigned by default. Check code!');
end

if ~isfield(options,'df_prctile') || isempty(options.df_prctile)
    options.df_prctile = defoptions.df_prctile;
end

if ~isfield(options,'df_window') || isempty(options.df_window)
    options.df_window = defoptions.df_window;
end

%reconstructed fluorescence from spatial, temporal, and residual components
F = diag(sum(A.^2))*(C + YrA);                                                     

% if not specific window size for smoothing or window size greater than #
% of frames (this should not run)
if isempty(options.df_window) || (options.df_window > size(C,2))
    Fd = prctile(F,options.df_prctile,2);
    F0 = repmat(prctile((A'*b)*f,options.df_prctile,2) + Fd,1,size(C,2));
    F_dff = (F - repmat(Fd,1,size(C,2)))./F0;
    
else
    % detrended fluorescence (Fd = F - baseline)
    %shift every frame (4th param); return F - baseline - 5th param
    Fd = prctfilt(F,options.df_prctile,options.df_window,1,1);
    
    %returns baseline of the background component (5th param = 0)
    F0_background = prctfilt((A'*b)*f,options.df_prctile,options.df_window,[],0);
    
    % background + baseline for each component
    %(F - Fd) = component baseline
    F0_baseline = F - Fd;
    
    F0 = F0_background + F0_baseline;  
    
    %dF/F = (F - F0_baseline)/(F0_baseline + F0_background)
    F_dff = Fd./F0;
    
    
end