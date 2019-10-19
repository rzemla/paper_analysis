function [trialType] = defineTrialSignal_punish(trialRanges, trialTypeCh, Behavior)


%% Load in position and time variables from Behavior struct

position = Behavior.position;
position_norm = Behavior.normalizedposition;
time = Behavior.time;
lap = Behavior.lap;

%% Find all of the trial signals

%find each microtexture 
for ii= 1:size(trialRanges,1)
    texOnIdx{ii} = find(trialTypeCh >= trialRanges(ii,1) & trialTypeCh <= trialRanges(ii,2)); %OK
    texOnsetOnlyIdx{ii} = find(diff(texOnIdx{ii}) >  1) + 1; %OK
    filterTexIndices{ii} = [1; texOnsetOnlyIdx{ii}];
    rawTexIndices{ii} =  texOnIdx{ii}(filterTexIndices{ii});
end

%scan the previous 2 indices and forward 2 indices to see that voltage
%does not exceed or fall below the range for that voltage - remove indices
%that have any values that exceed that range

removeTexIdx = cell(1,size(rawTexIndices,2));
%for each microtexture
for ii=1:size(rawTexIndices,2)
    for jj=1:size(rawTexIndices{ii},1)
        %get neighboring voltages for the each index associated with the
        %texture
        neighborVoltages{ii}{jj} = trialTypeCh((rawTexIndices{ii}(jj)-2):(rawTexIndices{ii}(jj)+2));
        %look at previous 2 voltage samples - should be less than 0.3V
        if sum((neighborVoltages{ii}{jj}(1:2) < trialRanges(ii,1))) ~=2
            removeTexIdx{ii} =  [removeTexIdx{ii}, jj];
        end
        %look at next 1 voltage samples should not be out of range
        if ((neighborVoltages{ii}{jj}(end-1) < trialRanges(ii,1)) | (neighborVoltages{ii}{jj}(end-1) > trialRanges(ii,2)))
            %add index for removal
            removeTexIdx{ii} =  [removeTexIdx{ii}, jj];
        end

    end
end

%copy the texture indices
rawTexIdxClean = rawTexIndices;

%remove false positive textures
for ii=1:size(rawTexIndices,2)
    %if found false positve textures for those indices
    if isempty(removeTexIdx{ii}) == 0
        rawTexIdxClean{ii}(removeTexIdx{ii}) = [];
    end
end


%% Find time and position of each texture
 
 %get norm position and offset position
 %for each texture
 for ii = 1:size(rawTexIndices,2)
    position_norm_tex{ii} = position_norm(rawTexIdxClean{ii});
    position_tex{ii} = position(rawTexIdxClean{ii});
    %median of position
    position_med(ii) = median(position(rawTexIdxClean{ii}));
    time_tex{ii} = time(rawTexIdxClean{ii});
    %voltage values at each index
    voltage_tex{ii} = trialTypeCh(rawTexIdxClean{ii});
 end
 
 %% Enforce order of texture 
 %high voltage A trials far - later) (index 1)
 %mid voltage B trials first near; (index 2)
 %low voltage punish trials (last) (index 3)
 
%get mean of voltages of the signal to determine type
 for ii =1:size(voltage_tex,2)
     mean_voltage(ii) = mean(voltage_tex{ii});
 end
 
% A trial is high voltage (2) and B trial is mid voltage (3); punish = 1
%get existing order
[~,idx_P] = min(mean_voltage);
[~,idx_A] = max(mean_voltage);
[~,idx_B] = find(median(mean_voltage) == mean_voltage);

%B, A, P
current_trial_order = [idx_B, idx_A, idx_P];

%rearrange cell order of trials such that B is first; A is second; P is last
position_norm_tex([1,2,3]) =  position_norm_tex([current_trial_order]);
position_tex([1,2,3]) =  position_tex([current_trial_order]);
time_tex([1,2,3]) =  time_tex([current_trial_order]);
position_med([1,2,3]) = position_med([current_trial_order]);
voltage_tex([1,2,3]) =  voltage_tex([current_trial_order]);


%trial name based on voltage output A or B
trialName{1} = 'Mid voltage (B trial - near)';
trialName{2} = 'High voltage (A trial - far)';
trialName{3} = 'Low voltage (punish)';

%% Only include trial signals that are within complete laps
%start time of 1st complete lap
start_full_lap_time = lap{1}(1);
%end time of last complete lap
end_full_lap_time = lap{end}(2);

%for first set of trials (B)
first_set_exclude = find(time_tex{1} < start_full_lap_time | time_tex{1} > end_full_lap_time);
%for second set of trials (A)
second_set_exclude = find(time_tex{2} < start_full_lap_time | time_tex{2} > end_full_lap_time);
%for third set of trials (P)
third_set_exclude = find(time_tex{3} < start_full_lap_time | time_tex{3} > end_full_lap_time);

%remove trial signals that are outside of complete laps (B)
if ~isempty(first_set_exclude)
    time_tex{1}(first_set_exclude) = [];
    position_tex{1}(first_set_exclude) = [];
    position_norm_tex{1}(first_set_exclude) = [];
end

%(A)
if ~isempty(second_set_exclude)
    time_tex{2}(second_set_exclude) = [];
    position_tex{2}(second_set_exclude) = [];
    position_norm_tex{2}(second_set_exclude) = [];
end

%(P)
if ~isempty(third_set_exclude)
    time_tex{3}(third_set_exclude) = [];
    position_tex{3}(third_set_exclude) = [];
    position_norm_tex{3}(third_set_exclude) = [];
end

%% Assign trial type info to struct

trialType.time = time_tex;
trialType.position = position_tex;
trialType.position_norm = position_norm_tex;
trialType.voltage = voltage_tex;
trialType.trialName = trialName;


end

