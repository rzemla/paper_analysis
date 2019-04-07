function [Events] = event_properties(updated_dff,Events,options)

%calculate the event properties of the dF/F signal
%think about ways of dealing with multipeaked events
onset_offset = Events.onset_offset;

%assign updated calcium trace  based on evetn detection ('Rolling median')
if strcmpi(options.dff_type, 'Jia/Danielson')
    dff = updated_dff.F_df_exp;
    
elseif  strcmpi(options.dff_type, 'Rolling median')
    dff= updated_dff.F_dff_exp;
    
elseif  strcmpi(options.dff_type, 'Rolling median (manual)')
    dff= updated_dff.F_df_exp;
end

%assign time of calcium trace
dff_time = Events.dff_time;


%% Measure calcium event properties (including multipeaked events)

%preallocate event traces
event_traces = cell(1,size(dff,2));

%imaging period (for AUC and event duration calculation)
%all values should be the same (mean = median = min = max)
%make sure dff does not arrive disjointed here through separate trials
imagingPeriod = median(diff(dff_time));

%run this for all, run, no_run epochs
for tt=1:3
    %change onset_offset variable to correspond to all, run, no_run events
    %in sequential order
    if tt == 1
        onset_offset = Events.onset_offset;
    elseif tt == 2
        onset_offset = Events.Run.run_onset_offset;
    elseif tt == 3
        onset_offset = Events.NoRun.norun_onset_offset;
    end
    
    %for each ROI
    for rr=1:size(dff,2)
        
        %check that event is no empty; if empty, set to empty
        if ~isempty(onset_offset{rr})
            %for each event
            for ee=1:size(onset_offset{rr},1)
                %extract each calcium trace that corresponds to each event
                event_traces{tt}{rr}{ee} = dff(onset_offset{rr}(ee,1):onset_offset{rr}(ee,2),rr);
                %event duration (imaging time (sec))
                event_dur{tt}{rr}(ee) = dff_time(onset_offset{rr}(ee,2)) - dff_time(onset_offset{rr}(ee,1));
                %mean dF/F of event
                event_mean{tt}{rr}(ee) = mean(event_traces{tt}{rr}{ee});
                %event AUC (dF/F*s)
                event_AUC{tt}{rr}(ee) = imagingPeriod*trapz(event_traces{tt}{rr}{ee});
                %event peak
                event_peak{tt}{rr}(ee) = max(event_traces{tt}{rr}{ee});
                %event amplitude (peak - onset dF/F)
                event_amp{tt}{rr}(ee) = event_peak{tt}{rr}(ee) - dff(onset_offset{rr}(ee,1),rr);
                %event width?
            end
            
        else
            %extract each calcium trace that corresponds to each event
            event_traces{tt}{rr} = [];
            %event duration (imaging time (sec))
            event_dur{tt}{rr} = [];
            %mean dF/F of event
            event_mean{tt}{rr} = [];
            %event AUC (dF/F*s)
            event_AUC{tt}{rr}= [];
            %event peak
            event_peak{tt}{rr} = [];
            %event amplitude (peak - onset dF/F)
            event_amp{tt}{rr} = [];
        end
        %number of events per ROI
        nb_events{tt}(rr) = size(onset_offset{rr},1);
        
    end
    
    %total event count - should be the same as # of statistically significant
    %events
    nb_events_total{tt} = sum(nb_events{tt});
    
end

%% Section for non-analyzed events and multi-peaked event detection/separation?


%% Save to structure
%properties of all events (no distinction between run and no-run events
Events.properties.trace = event_traces{1};
Events.properties.duration = event_dur{1};
Events.properties.peak = event_peak{1};
Events.properties.amplitude = event_amp{1};
Events.properties.mean = event_mean{1};
Events.properties.AUC = event_AUC{1};
Events.properties.nb_events = nb_events{1};
Events.properties.nb_events_total = nb_events_total{1};

%run events
Events.Run.properties.trace = event_traces{2};
Events.Run.properties.duration = event_dur{2};
Events.Run.properties.peak = event_peak{2};
Events.Run.properties.amplitude = event_amp{2};
Events.Run.properties.mean = event_mean{2};
Events.Run.properties.AUC = event_AUC{2};
Events.Run.properties.nb_events = nb_events{2};
Events.Run.properties.nb_events_total = nb_events_total{2};

%norun events
Events.NoRun.properties.trace = event_traces{3};
Events.NoRun.properties.duration = event_dur{3};
Events.NoRun.properties.peak = event_peak{3};
Events.NoRun.properties.amplitude = event_amp{3};
Events.NoRun.properties.mean = event_mean{3};
Events.NoRun.properties.AUC = event_AUC{3};
Events.NoRun.properties.nb_events = nb_events{3};
Events.NoRun.properties.nb_events_total = nb_events_total{3};

%Events.properties.width=event_width;
end

