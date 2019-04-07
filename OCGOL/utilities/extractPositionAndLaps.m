function [Behavior] = extractPositionAndLaps(CSV,options)

%V3 - adds behavioral tag extraction in each lap

%% Assign parameters

%distance min before counting a new lap (s)
mindist=options.mindist;

%frequency at which behavioral voltage was acquired
acqfr=options.acqHz; % behavior acquisition frequency (Hz)

%CSV matrix with behavioral voltages
inputCSV=CSV;

%% Separate time vector from signal channels and digitize voltages
%1 position step = 0.155986  mm = 0.000

%Isolate signals from time (remove the first time column)
csvStimRaw = inputCSV(:,2:end);

%Digitize the voltage signals
%set ~+5V high voltage to 1 and ~0V low voltage to 0
csvStimRaw(find((csvStimRaw > 1))) = 1;
csvStimRaw(find((csvStimRaw < 1))) = 0;

%take out digitized A and B signal
csvABdigitized = csvStimRaw(:,3:4);

%Place time in separate variable and convert to seconds
%time in miliseconds
timeOutput = inputCSV(:,1);
%convert from ms to sec
timeOutputSec = timeOutput./1000;

%% Decode position from quadrature A/B rotary encoder channels

%past logical state of rotary channels - initialization to algorithm
pastStateA = csvABdigitized(1,1);
pastStateB = csvABdigitized(1,2);

%current state - starting at the second position
currentStateA = 0;
currentStateB = 0;

%position counter for algorithm
positionCount = 0;
%vector giving the position of the animal at each timepoint
positionVector = zeros(size(csvABdigitized,1),1);

%position decoder

for idx = 2:size(csvABdigitized,1)
    currentStateA = csvABdigitized(idx,1);
    currentStateB = csvABdigitized(idx,2);
    if (currentStateA == pastStateA && currentStateB == pastStateB)
        positionVector(idx) =  positionCount;
    elseif (currentStateA ~= pastStateA || currentStateB ~= pastStateB)
        if (currentStateA == 0 && currentStateB == 0)
            if pastStateA == 0
                positionCount = positionCount + 1;
            elseif pastStateA == 1
                positionCount = positionCount - 1;
            end
            positionVector(idx) =  positionCount;
        elseif (currentStateA == 1 && currentStateB == 0)
            if pastStateA == 0
                positionCount = positionCount + 1;
            elseif pastStateA == 1
                positionCount = positionCount - 1;
            end
            positionVector(idx) =  positionCount;
        elseif (currentStateA == 1 && currentStateB == 1)
            if pastStateA == 0
                positionCount = positionCount - 1;
            elseif pastStateA == 1
                positionCount = positionCount + 1;
            end
            
            positionVector(idx) =  positionCount;
            
        elseif (currentStateA == 0 && currentStateB == 1)
            if pastStateA == 0
                positionCount = positionCount - 1;
                
            elseif pastStateA == 1
                positionCount = positionCount + 1;
            end
            
            positionVector(idx) =  positionCount;
        end
        
    end
    pastStateA = currentStateA;
    pastStateB = currentStateB;
end

%convert position indices to centimeters
cum_position=positionVector*0.0155986;

%% Plot cumulative and converted position against time

figure;
subplot(4,1,1)
hold on
ylabel('Encoder Channel A')
plot(timeOutputSec,csvABdigitized(:,1),'r')
hold off

subplot(4,1,2)
hold on
ylabel('Encoder Channel B')
plot(timeOutputSec,csvABdigitized(:,2),'b')
hold off

subplot(4,1,3)
hold on
ylabel('Position ')
plot(timeOutputSec,positionVector,'k')
hold off

subplot(4,1,4)
hold on
ylabel('Position (cm) ')
xlabel('Time (s)')
plot(timeOutputSec,cum_position,'k')
hold off


%% Separate cumulative position into laps

%separate logical signals corresponding to lap crossing
lap = csvStimRaw(:,8);

%returns all indices where the RFID tag signal is high
lap_RFID = find(lap >= 1);

%get cm positions of where the RFID tag is high
RFID_POS = cum_position(lap_RFID);

%select only tag positions with greater than 10 cm separatation
RFID_keep = [1; diff(RFID_POS) >= mindist];
%convert to logical vector
RFID_keep = logical(RFID_keep);

%get absolute cm position of each RFID crossing 
RFID_POS_keep = RFID_POS(RFID_keep);

%get absolute lap indices where RFID tag onset begings
lap_RFID_keep = lap_RFID(RFID_keep);

%% Convert lap onset to binary representation (vector)

%create a zero vector with size equal to length of recording session
lap_binary = zeros(size(lap,1),1);

%set onset indices of each lap to 1
lap_binary(lap_RFID_keep) = 1;

%copy indices of lap onset
lap_row = lap_RFID_keep;

%copy the equivalent time in seconds
time = timeOutputSec;


%% Insert plotter here


%% Segment the laps

%may need to increase 1 point separation here so there is no overlap
%between laps - 0.1 ms difference that should have significant effect
%however
%start of lap - first timepoint when tag goes high
%end of lap - 1 point right before tag goes high 

%if at least 1 lap completed
if length(lap_row)>1 
    %display number of complete laps
    disp(['Total laps run: ' num2str(length(lap_row)-1)])
    
    %for each complete lap (total RFID crossings -1)
    for i=1:(length(lap_row)-1)
        
        %get absolute position in cm at each time between two lap
        %crossings
        pos_cumul{i}=cum_position(lap_row(i):(lap_row(i+1)-1),1);
        
        %take sbsolute starting position (cm) of each lap and subtract
        %get relative position from the start of the lap
        pos_reset{i}=pos_cumul{i}-cum_position(lap_row(i));
        
        %normalize postion from 0 - 1
        norm_pos{i} = (pos_reset{i} - min(pos_reset{i})) / ( max(pos_reset{i}) - min(pos_reset{i}) );
        
    end

    
% length(cell2mat(pos_cumul'))
% length(cum_position(lap_row(1):lap_row(end)-1))

%start of lap is at crossing of RFID tag
%absolute indices from recording session
lap_start = lap_row(1:(length(lap_row)-1));
lap_start = lap_start';

%lap stop endpoints (1 index before tag onset)
lap_stop = lap_row-1;
lap_stop = lap_stop';
    
%%  Add first and last incomplete lap information

    %Absolute position for the first incomplete lap (until 1st tag)
    lap1 = cum_position(1:lap_row(1)-1,1);
    
    %Position for the last incomplete lap - relative since start of last
    %lap (added a 1 shift here)
    lastlap = cum_position((lap_row(end)):end)-cum_position(lap_row(end));
    
    %normalize last incomplete lap relative to last complete lap
    lastlap_norm=(lastlap - min(pos_reset{end})) / ( max(pos_reset{end}) - min(pos_reset{end}) );
    
    %normalize first incomplete lap relative to first complete lap - RZ
    lap1_norm=(lap1 - min(pos_reset{1})) / ( max(pos_reset{1}) - min(pos_reset{1}) );
    
    %normalize relative to the end of the first complete lap - RZ
    lap1_norm=lap1_norm+1-max(lap1_norm);
    
    %convert lap position in relation to the last position from the next
    %complete lap
    %get all end positions
    for ii=1:size(pos_reset,2)
        endLapPositions(ii) = pos_reset{ii}(end);
    end
    
    %median end postion
    endPosMed = median(endLapPositions);
    
    %arrange end position in relation to of median of end lap positions
    %(cm)
    lap1 = lap1 + endPosMed - max(lap1);
    
    %add the first and last incomplete laps - non-normalized, but scaled
    %ralative to start of each lap in cm
    pos_reset=[lap1 pos_reset lastlap];
    
    %concatenate normalized positions together
    norm_pos=[lap1_norm norm_pos lastlap_norm];
    
    %check if concatenation of the index equals in size to recording
    %session
%      x = cell2mat(pos_reset');   
%      x = cell2mat(norm_pos'); 

%% Get table with start and stop indices for each complete lap

    %only endpoints of complete laps
    lap_stop=lap_stop(2:end);
    
    %convert indices to time in seconds
    lap_start_time = time(lap_start)';
    lap_stop_time = time(lap_stop)';
    
    %for each complete lap
    for i=1:length(lap_start)
        
        %start and stop times for each lap in a separate cell
        lap_start_stop_time{i}=[lap_start_time(i) lap_stop_time(i)];
    end
    
    %assign a lap index to each timepoint
    for ii = 1:size(pos_reset,2)
        %create empty ones vector corresponding to # recorded indices in that lap
        %that has a number corresponding to the lap sequence
        %0 index is first incomplete laps
       lapNb_idx{ii} = (ii-1)*ones(size(pos_reset{ii},1),1);
    end   
    
    %reshape relative cm position into cells that can be converted to
    %vector
    B=reshape(pos_reset,[],1);
    
    %do the same for normalized position
    N=reshape(norm_pos,[],1);
    
    %reshape the lap indices - RZ
    L=reshape(lapNb_idx,[],1);
    
    %convert to continuous vectors
    %original # of samples acquired
    position = cell2mat(B);
    
    %do the same for the normalized position
    position_norm = cell2mat(N);
    
    %convert lap associated with each index to continuous column vector
    lapNb = cell2mat(L);
    
end

%% Quality Control %%

%% Check that each lap is >190 cm and < 200 cm

%Get lap position (cm) at each start and stop time point
%complete laps
for ii=1:size(lap_start_stop_time,2)
    %time index of lap start
    time_idx_lap{ii}(1) = find(time == lap_start_stop_time{ii}(1));
    %time index of lap stop
    time_idx_lap{ii}(2) = find(time == lap_start_stop_time{ii}(2));
    %start and stop position of each lap
    position_lap(ii,1) = position(time_idx_lap{ii}(1));
    position_lap(ii,2) = position(time_idx_lap{ii}(2));
end

%total lap distance
position_lap_diff = position_lap(:,2) - position_lap(:,1);

%find laps that have distances out of expected boundaries
out_of_bounds_lap_idx = find(position_lap_diff < 190 & position_lap_diff > 200);

%generate warning if any laps out-of-bounds
if isempty(out_of_bounds_lap_idx)
    disp('All laps are within bounds');
    allLapsInBounds = 1;
else
    disp('Some laps are out of bounds !!!');
    allLapsInBounds = 0;
    %save to struct lap indices that are out of bounds
    Behavior.out_of_bounds_lap_idx = out_of_bounds_lap_idx;
end


%save raw lick signal across recording to lick variable
lick = csvStimRaw(:,7);
 
%% Make structure

%flag for whether all laps are in bound (no lap tag missed)
Behavior.allLapsInBounds = allLapsInBounds;

%start and stop position of each lap
Behavior.position_lap = position_lap;

%time index of each lap (start and stop)
Behavior.time_idx_lap = time_idx_lap;

%absolute time vector
Behavior.time = time;

%cumulative offset position after each lap crossing
Behavior.position = position;

%raw lick signal
Behavior.lick = lick;

%lab number corresponding to each index
Behavior.lapNb = lapNb;

%binary high signals when lap crossing occurs on behavioral time domain
Behavior.lap_binary = lap_binary;

%input options that are output - subcategorize (into struct) each options for
%for each function
Behavior.options = options;

%cumulative position since the start of the recording
Behavior.cumulativeposition = cum_position;

%at least one lap has to be run by the mouse
%lap_row corresponds to all the lap cross indices 
if length(lap_row)>1
    %save the  vector of normalized positions
    Behavior.normalizedposition = position_norm;
    
    %save the corresponding absolute time start and time stop for each lap
    Behavior.lap = lap_start_stop_time;
    
    %duplicate from above- can remove - RZ - CHECK
    Behavior.normalizedposition = position_norm;
end

%% Display figure

%if at least 1 complete lap run
if length(lap_row)>1
    
    %display figure option
    if options.dispfig==1
        
        %plot all laps
        figure; 
        subplot(2,1,1)
        hold on
        title('Position');
        ylabel('Position [cm]');
        plot(Behavior.time,Behavior.position,'k');
        hold off
        
        subplot(2,1,2)
        hold on
        title('Normalized position');
        %normalized position to each lap
        plot(Behavior.time,Behavior.normalizedposition,'k');
        
        %plot RFID crossing indices against time
        plot(Behavior.time,lap_binary, 'r');
        
        %plot the lap index - lap # scaled down by a factor of 20 - each
        %step 0.05
        plot(Behavior.time,Behavior.lapNb./20,'r');
        
    end
end

end



