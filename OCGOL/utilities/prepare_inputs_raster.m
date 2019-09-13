function [plot_raster_vars] = prepare_inputs_raster(session_vars,CNMF_vars,removeROI)


%% Define/load variables for each session

%for each session
for ii = 1:size(session_vars,2)
    % behavior and imaging related variables
    Behavior_split_lap{ii} = session_vars{ii}.Behavior_split_lap;
    Events_split_lap{ii} = session_vars{ii}.Events_split_lap;
    Behavior_split{ii} = session_vars{ii}.Behavior_split;
    Event_split{ii} = session_vars{ii}.Events_split;
    Imaging_split{ii} = session_vars{ii}.Imaging_split;
    Place_cell{ii} = session_vars{ii}.Place_cell;
    Behavior_full{ii} = session_vars{ii}.Behavior;
    
    %all within run domain
    position{ii}  = Behavior_split_lap{ii}.Run.position;
    time{ii} = Behavior_split_lap{ii}.Run.time;
    events_full{ii} = Events_split_lap{ii}.Run.run_onset_binary;
    run_intervals{ii} = Behavior_split_lap{ii}.run_ones;
    
    %global trial type order across restricted laps
    trialOrder{ii} = Behavior_full{ii}.performance.trialOrder;
end

%for each session
for ss=1:size(session_vars,2)
    %for each lap
    for ii=1:size(run_intervals{ss},2)
        events{ss}{ii} = events_full{ss}{ii}(logical(run_intervals{ss}{ii}),:);
    end
end

%lap_times = Behavior_split_lap.lap;
%how many trial types are there
%set to 2 - A and B for now
%trialTypes = 2;

%make this a condition for the type of trials that are compared
%all correct; all regardless of correct

%for each session (all A or B regardless of correct)
for ss = 1:size(session_vars,2)
    %find times from trials related to specific trial type
    %A trials
    trialTypeIdx{ss}{1} = find(trialOrder{ss} == 2 | trialOrder{ss} == 20);
    %B trials
    trialTypeIdx{ss}{2} = find(trialOrder{ss} == 3 | trialOrder{ss} == 30);
end

%get the chosen spatial components of each ROI (A matrix from CNMF output)
A_keep_sel = CNMF_vars{1, 1}.A_keep(:,removeROI{1}.compSelect);  

%center of mass of each component
A_keep_sel_com = com(A_keep_sel,512,512);

%outline
coor_keep_sel = CNMF_vars{1}.Coor_kp(removeROI{1}.compSelect);


%% Define the spiral parameters according to the number of laps
%equivalent to number of laps for each session (turns = laps)
for ss = 1:size(session_vars,2)
    turns(ss) = size(trialOrder{ss},1); %The number of turns the spiral will have (how many laps)
    
    %x is the angle
    x{ss} = [-1*pi*turns(ss) : 0.01 : pi*turns(ss)];
    
    %r is that radius of the point
    r{ss} = [0:1/(length(x{ss})-1):1];
    
    %scale to lap length
    r_scaled{ss} = r{ss}.*turns(ss);
end

%all parameters in the run frame domain
%find the frames index of event and position
%for each session
for ss=1:size(session_vars,2)
    %for trial type (A or B)
    for ii=1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll= 1:size(trialTypeIdx{ss}{ii},1)
            %for each ROI
            for rr=1:size(events{ss}{trialTypeIdx{ss}{ii}(ll)},2)
                %event indices in run domain
                event_idx{ss}{ii}{ll}{rr} = find(events{ss}{trialTypeIdx{ss}{ii}(ll)}(:,rr) == 1);
                %position that corresponds to indices
                pos{ss}{ii}{ll}{rr} = position{ss}{trialTypeIdx{ss}{ii}(ll)}(event_idx{ss}{ii}{ll}{rr});
                %position vectors that will be used as input to spiral
                posVectors{ss}{ii}{ll}{rr} = trialTypeIdx{ss}{ii}(ll).*exp(1i.*((pos{ss}{ii}{ll}{rr}/200)*2*pi)).';
            end
        end
    end
end

%make the nearest approximation to a point along the spiral vector defined above
%predefine to avoid empty cells at the end
valMin = posVectors;
idxMin = posVectors;
posVectorApprox = posVectors;

%foreach session
for ss =1:size(session_vars,2)
    %for each trial type
    for kk = 1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll = 1:size(pos{ss}{kk},2)
            %for each ROI
            for rr = 1:size(events{ss}{1},2)
                %for each event
                for ee=1:size(pos{ss}{kk}{ll}{rr},1)
                    %fix empty cell clipping at endpoints
                    [valMin{ss}{kk}{ll}{rr}(ee),idxMin{ss}{kk}{ll}{rr}(ee)] = min(abs( (r_scaled{ss} - ( (trialTypeIdx{ss}{kk}(ll)-1) + (pos{ss}{kk}{ll}{rr}(ee)/200) ) ) ));
                    posVectorApprox{ss}{kk}{ll}{rr}(ee) = r_scaled{ss}(idxMin{ss}{kk}{ll}{rr}(ee))*exp(1i.*(pos{ss}{kk}{ll}{rr}(ee)/200)*2*pi);
                end
            end
        end
    end
end

%% Split calcium traces into intervals to avoid joining by plot into continuous plot

%breakpoints
diff_A_trials_idx = find(diff(Imaging_split{1}{4}.time_restricted)> 1);
diff_B_trials_idx = find(diff(Imaging_split{1}{5}.time_restricted)> 1);

%construct start and end idx matrices - correct
start_idx.A = [1; diff_A_trials_idx+1];
end_idx.A = [diff_A_trials_idx; size(Imaging_split{1}{4}.time_restricted,1)];
%combined start and end idx
start_end_idx.A = [start_idx.A, end_idx.A]; 
%for B trials
start_idx.B = [1; diff_B_trials_idx+1];
end_idx.B = [diff_B_trials_idx; size(Imaging_split{1}{5}.time_restricted,1)];
%combined start and end idx
start_end_idx.B = [start_idx.B, end_idx.B]; 

%% Get lap indices for each lap in all B or B trials

%get unique lap indices
lapA_idxs = unique(Behavior_split{1}{4}.resampled.lapNb);
lapB_idxs = unique(Behavior_split{1}{5}.resampled.lapNb);

%get lap start and end indices for all A or B trials
%all A
for ll=1:size(lapA_idxs,1)
    lap_idxs.A(ll,1) = find(Behavior_split{1}{4}.resampled.lapNb == lapA_idxs(ll),1,'first');
    lap_idxs.A(ll,2) = find(Behavior_split{1}{4}.resampled.lapNb == lapA_idxs(ll),1,'last');
end

%all B
for ll=1:size(lapB_idxs,1)
    lap_idxs.B(ll,1) = find(Behavior_split{1}{5}.resampled.lapNb == lapB_idxs(ll),1,'first');
    lap_idxs.B(ll,2) = find(Behavior_split{1}{5}.resampled.lapNb == lapB_idxs(ll),1,'last');
end

%% Event onsets in run interval
%for each ROI
for rr=1:size(events{ss}{1},2)
    %time of significant run events in A
    event_norm_time.A{rr} =Imaging_split{1}{4}.time_restricted(find(Event_split{1}{4}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in A
    event_norm_pos_run.A{rr} = Behavior_split{1}{4}.resampled.position_norm(find(Event_split{1}{4}.Run.run_onset_binary(:,rr) == 1));
    %time of significant run events in B
    event_norm_time.B{rr} =Imaging_split{1}{5}.time_restricted(find(Event_split{1}{5}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in B
    event_norm_pos_run.B{rr} = Behavior_split{1}{5}.resampled.position_norm(find(Event_split{1}{5}.Run.run_onset_binary(:,rr) == 1));
end

%% Export variables as part of struct

plot_raster_vars.start_end_idx = start_end_idx;
plot_raster_vars.lap_idxs = lap_idxs;

plot_raster_vars.x = x;
plot_raster_vars.r_scaled = r_scaled;

plot_raster_vars.idxMin = idxMin;
plot_raster_vars.posVectorApprox = posVectorApprox;

plot_raster_vars.event_norm_time = event_norm_time;
plot_raster_vars.event_norm_pos_run = event_norm_pos_run;


plot_raster_vars.A_keep_sel_com = A_keep_sel_com;
plot_raster_vars.coor_keep_sel = coor_keep_sel;

end

