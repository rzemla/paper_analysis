function [Imaging, Behavior, F_vars]= restrict_data_training(Behavior, XML, options)

% if alignment not restricted to full laps only - RZ;  restrict trace to full lap
if options.restrict == 0 
    
    %assign the whole C_df matrix as the trace variable for the struct
    %output; frame x ROI
    %Imaging.trace = C_df;
    
    %get framelength of imaging session from XML - RZ
    %length should be 1 frames longer than C_df trace as the first frame is
    %removed from the CNMF analyzed and MC data
    frameLength = length(XML.PVScan.Sequence.Frame);
    
    %preallocate array - saves time for each read from XML by preallocating
    %memory - RZ
    timeStampsXML = zeros(frameLength,1);
    
    %extract the relative timestamps for each frame acquired in session
    %iterate through each XML index/Frame structure and extract relative
    %time of each frame - convert to double and save as vector
    %time is in seconds
    
    for ii=1:frameLength
        timeStampsXML(ii,1) = str2double(XML.PVScan.Sequence.Frame{1,ii}.Attributes.relativeTime);
    end
    
    %end frame based on XML parsing
    %remove the first timepoint associated with the first frames that is
    %removed for CNMF analysis
    imaging_time = timeStampsXML(2:end);
    
    %save the C_df aligned time (with first frame removed) into imaging
    %struct
    Imaging.time = imaging_time;
    
    %time step based on imaging timepoints - pass on this dt to
    Imaging.dt = median(diff(Imaging.time));
    
    %save the associated options into the Behavior matrix
    Behavior.options = options;
    
    %import corresponding lap indices
    lapNb=Behavior.lapNb;
    
end

% restict traces for full laps
if options.restrict==1 
    
    %% Import data
    
    %start and stop time for each complete lap
    lap_start_stop = Behavior.lap;
    
    %absolute behavior time
    behavior_time = Behavior.time;
    
    %position (cm)
    position=Behavior.position;
    
    %import corresponding lap indices on behavior time domain
    lapNb=Behavior.lapNb;
    
    %same as above, but the position is normalzied from 0 - 1
    norm_position = Behavior.normalizedposition;
    
    %cumulative position since the beginning of the lap
    cum_position = Behavior.cumulativeposition;
    
    %load the dF/F signal into the imaging struct
    %Imaging struct gets set up here
    %Imaging.trace = C_df;
    
    %select which complete lap is the start lap for analysis
    switch options.startlap
        %first complete lap
        case 'first'
            startlap = 1;  
        %which lap is start lap
        otherwise
            startlap = options.startlap;
    end
    
    %select which lap is the end lap for the analysis
    switch options.endlap
        %last complete lap
        case 'last'
            endlap = size(Behavior.lap,2);
        
        %use the last lap that is selected
        otherwise
            endlap=options.endlap;
    end
    
    %import texture locations if available
    %modify this
%     if options.textures==true
%         texture=Behavior.texture_signal;
%     end
%     
    
    %% Extract Timestamps from XML file
    %total # of frames (including first that is excluded in MC/CNMF)
    frameLength = length(XML.PVScan.Sequence.Frame);
    
    %preallocate array for imaging ttimestamps
    timeStampsXML = zeros(frameLength,1);
    
    %extract timestamps (s) including first
    for ii=1:frameLength
        timeStampsXML(ii,1) = str2double(XML.PVScan.Sequence.Frame{1,ii}.Attributes.relativeTime);
    end
    
    %trim first timepoint and save in imaging_time variable
    imaging_time = timeStampsXML(2:end);
    
    %% Restrict imaging time and traces
    
    % Set start and end time based on first and last lap
    %take the start time correspond to the first complete lap (from
    %behavior)
    startT = lap_start_stop{startlap}(:,1);
    
    %take the end time corresponding to the last complete lap 
    endT = lap_start_stop{endlap}(:,2);
    
    %indices for lap-restricted imaging times 
    t_in = find(imaging_time >= startT & imaging_time <= endT);
    
    % restrict calcium trace using t_in indices
    %C_df_R = C_df(t_in,:);
    
    %restrict all CNMF/F variables for later dF/F recalculation
    %raw
%     F_vars.restricted.F = F_vars.F(:,t_in);
%     
%     %F0's
%     F_vars.restricted.F0_background = F_vars.F0_background(:,t_in);
%     F_vars.restricted.F0_baseline = F_vars.F0_baseline(:,t_in);
%     F_vars.restricted.F0 = F_vars.F0(:,t_in);
%     F_vars.restricted.Fd = F_vars.Fd(:,t_in);
%     
%     %dF/F
%     F_vars.restricted.F_dff = F_vars.F_dff(:,t_in);
%     F_vars.restricted.F_dff_exp = F_vars.F_dff_exp(:,t_in);
    
    % restrict imaging time
    imaging_time_R = imaging_time(t_in);

    %% Restrict behavior time and voltages
    
    % indices for restricted behavior time
    %b_in - indices of the corresponding timepoints
    b_in = find(behavior_time >= startT & behavior_time <= endT);
    
    %restrict beahvior time
    behavior_time_R = behavior_time(b_in);
    
    %restrict behavior position (cm)
    position_R = position(b_in);
    
    %restrict normalized position
    norm_position_R = norm_position(b_in);
    
    %restrict cumulative position
    cum_position_R = cum_position(b_in);
    
    %restrict lap idx in behavior time domain
    lapNb_R = lapNb(b_in);
    
    %restrict texture
    %WORK ON THIS
%     if options.textures == true
%         texture_R = texture(b_in,:);
%         Behavior.restricted.texture_signal=texture_R;
%     end
    
    % Restrict laps (same unless narrowed range is chosen
    for ii = startlap:endlap
        lap_start_stop_restricted{ii} = lap_start_stop{ii};
    end
    
    %check to see if there are any empty cells that do not have a start or
    %end and remove them - should not occur 
    lap_start_stop_restricted = lap_start_stop_restricted(~cellfun(@isempty, lap_start_stop_restricted));
    
    %% Make structure
    
    %time corresponding to complete laps
    Behavior.restricted.time=behavior_time_R;
    
    %offset position correspond to complete laps
    Behavior.restricted.position=position_R;
    
    %lap index for complete laps
    Behavior.restricted.lapNb = lapNb_R;
    
    %absolute cumulative position correspondiing to complete laps
    Behavior.restricted.cumulativeposition=cum_position_R;
    
    %normalized (0-1) position corresponding to complete laps
    Behavior.restricted.normalizedposition=norm_position_R;
    
    %start and stop times for each lap for complete laps
    Behavior.restricted.lap =lap_start_stop_restricted;
    
    %time corespondong to frames imaged (with first frame removed)
    Imaging.time = imaging_time;
    
    %time step based on imaging timepoints - pass on this dt to
    Imaging.dt = median(diff(Imaging.time));
    
    %frame imaging times that correspond to complete laps
    Imaging.time_restricted = imaging_time_R;
    
    %imaging signal that that corresponds to complete laps
    %Imaging.trace_restricted = C_df_R;
    
    %export the options for this script to Behavior struct
    Behavior.options=options;
    
    %% Plot
    
    figure;
    if options.dispfig==true
        subplot(2,1,1)
        %plot the entire calcium signal with respective imaging times along with
        %entire normalized position of the animal with respective
        %behavioral time for selected cell - c2plot defined in main script
        plot(imaging_time, C_df(:,options.ROI)); 
        hold on; 
        title('Non-restricted data');
        plot(Behavior.time, Behavior.normalizedposition)
        %set max of x axis to max behavior time
        xlim([0 max(Behavior.time)]);
        hold off
        
        subplot(2,1,2)
        %plot restricted calcium trace with restricted normalized position for
        %selected cell in c2plot
        hold on
        title('Restricted data');
        plot(Imaging.time_restricted, Imaging.trace_restricted(:,options.ROI)); hold on; plot(Behavior.restricted.time, Behavior.restricted.normalizedposition)
        %set max of x axis to max behavior time
        xlim([0 max(Behavior.time)]);
        hold off
        
    end
    
end

end
