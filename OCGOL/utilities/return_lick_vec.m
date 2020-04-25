function [lick_vec] = return_lick_vec(session_vars,ses)

%constuct binary lick vector from time and position
%this all lick data - need to align this data since it has not been
%resampled for performance analysis

%align these event to the imaging time
x = session_vars{ses}.Behavior.lick.time;
y = session_vars{ses}.Behavior.resampled.time;


%filter to not go beyond the edges of the complete lap data timepoints
%clip off timepoints that are outside of the complete lap time series
clip_idx = ~(x<y(1)| x>y(end));
lick_times_clipped =x(clip_idx);

%return the index for each lick timepoint that is closest to the imaging
%time
%find the minimum lick index
min_diff = @(f)min(abs(f-y));
[~,z] = arrayfun(min_diff,lick_times_clipped,'UniformOutput',true);

%generate binary lick vector of equivalent length to imaging frames
lick_vec = zeros(size(y,1),1);
%populate the vector
lick_vec(z) = 1;

end

