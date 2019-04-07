function [Behavior] = OCGOL_performance(Behavior,CSV)

%Get performance of the animal on the OCGOL trials

%reward zone size (cm)
rewardZoneSize = 10;

%range of y axis of histogram
yRange = 1;

%bin size
binSize = 2;

%bin location of rewards (cm)
binReward_B = 59;
binReward_A = 139;

%calculate performance based on reward location recorded on each lap vs.
%median/mean of those locations?

%% lick index, position, and time 
%raw CSV lick signal
lickSignal = CSV(:,8);

%votlage range for lick on
lickRange = [3.5, 4.2];

%find all indices within range
%these correspond to the lick indices in the CSV cell
idx = find(lickSignal > lickRange(1) & lickSignal < lickRange(2));

%digitize the lick signal
lickSignalDigital = zeros(1,size(lickSignal,1));
lickSignalDigital(idx) = 1;

%indices of lick onset
%takes only the indices that go from low to high
diffLickSignal = diff(lickSignalDigital);
%get indices and compensate for 1 point offset from diff output
%indices correspond to the raw Lick signal (non-restricted)
keepIdxOnset = find(diffLickSignal>0) + 1;

%keepIdx = [1; find(diff(idx)>1)+1];

%keep only onset of each lick (count);
%idxFilter = idx(keepIdx);
idxFilter = keepIdxOnset';
%lick position time matrix concatenated laps
lickMat = idxFilter;

%find position associated with each lick onset
lickMat(:,2) = Behavior.position(idxFilter);

%find absolute time associated with each lick (3rd column)
lickMat(:,3) = Behavior.time(idxFilter);

%round the bins and sort
%binFilterSort = sort(round(lickMat(:,2)));

%% reward location - similar to lick location extraction%%%%%%%
%%%find location where animal received the rewards associated with
%%%reward zone

%voltage range for reward location signal from channel 2 (DAC #1)
rewardLocRange_Ch2 = [3.04 3.09];
%isolate the raw voltage trace from the DAC channel #1 (all texture
%signals)
textureSignal = CSV(:,3);

%find all indices within range that correspond to the reward location in the CSV cell
%reward location indices
idxRL = find(textureSignal > rewardLocRange_Ch2(1) & textureSignal < rewardLocRange_Ch2(2));

%binarize the lick signal
textureRsignalDigital = zeros(1,size(textureSignal,1));
textureRsignalDigital(idxRL) = 1;

%takes only the indices that go from low to high
diffRewardLocSignal = diff(textureRsignalDigital);
%get indices and compensate for 1 point offset from diff output
%indices correspond to the reward location on index axis (non-restricted)
keepIdxOnsetRewLoc = find(diffRewardLocSignal>0) + 1;

rewardLocMat(:,1) = keepIdxOnsetRewLoc';
%find bin associated with each reward location
rewardLocMat(:,2) = Behavior.position(keepIdxOnsetRewLoc');

%time association with each reward - when RFID location was crossed
rewardLocMat(:,3) = Behavior.time(keepIdxOnsetRewLoc');

%clean up reward location signals
removeRewIdx = [];
for rr=1:size(rewardLocMat,1)
    %get neighboring voltages for the each index associated with the
    %texture
    neighborVoltagesRew{rr} = textureSignal((rewardLocMat(rr,1)-2):(rewardLocMat(rr,1)+2));
    
    %look at previous 2 voltage samples - should be less than 0.3V
    if sum((neighborVoltagesRew{rr}(1:2) < rewardLocRange_Ch2(1))) ~=2
        removeRewIdx =  [removeRewIdx, rr];
    end
    %look at next 1 voltage samples should not be out of range
    if ((neighborVoltagesRew{rr}(end-1) < rewardLocRange_Ch2(1)) || (neighborVoltagesRew{rr}(end-1) > rewardLocRange_Ch2(2)))
        %add index for removal
        removeRewIdx =  [removeRewIdx,rr];
    end
    
end

%take unique indices
removeRewIdx = unique(removeRewIdx);

%copy the texture indices
rewardLocMatClean = rewardLocMat;

%remove false positive remward locations
if isempty(removeRewIdx) == 0
    rewardLocMatClean(removeRewIdx,:) = [];
end

%assign the trial type to the 4th column - 2 if A trial (far location);
%3 if B trial (near location)
%B trials
%100 is halfway through the track
rewardLocMatClean(find(rewardLocMatClean(:,2) <100),4) = 3;
rewardLocMatClean(find(rewardLocMatClean(:,2) >100),4) = 2;

%structure with reward location statistics from each lap/trial
RewardLoc.mat = rewardLocMatClean;
RewardLoc.matLabels = {'Idx reward','Location reward', 'Time reward','Trial type'};
RewardLoc.medALoc = median(rewardLocMatClean(find(rewardLocMatClean(:,4) == 2),2));
RewardLoc.medBLoc = median(rewardLocMatClean(find(rewardLocMatClean(:,4) == 3),2));

%assign reward bin locations based on DAC output vs estimated output
binReward_B = RewardLoc.medALoc;
binReward_A = RewardLoc.medBLoc;

%rewardRange in cm (10 cm)
binRewardRange1 = [binReward_B, (binReward_B+rewardZoneSize)];
binRewardRange2 = [binReward_A, (binReward_A+rewardZoneSize)];

%anticipatory range in cm
binAntRange1 = [(binReward_B -1 -rewardZoneSize),(binReward_B -1)];
binAntRange2 = [(binReward_A -1 -rewardZoneSize),(binReward_A -1)];

%trial order based on reporting of reward site from RFID tag
trialOrderTag = RewardLoc.mat(:,4)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% trial type based on DAC #2 output
%%%%%reward location - similar to lick location extraction%%%%%%%

%voltage range for reward location signal from channel 2 (DAC #1)
trialTypeA_Range_Ch1 = [4.24 4.29];
trialTypeB_Range_Ch1 = [2.43 2.49];

%concatenate above into one matrix for looping purposes
trialTypeRange = [trialTypeA_Range_Ch1;trialTypeB_Range_Ch1];

%isolate the raw voltage trace from the DAC channel #2(output 2 trial
%signal)
trialSignal = CSV(:,2);

%find all indices within range that correspond to the trial type in the CSV cell
idxTA = find(trialSignal > trialTypeA_Range_Ch1(1) & trialSignal < trialTypeA_Range_Ch1(2));
idxTB = find(trialSignal > trialTypeB_Range_Ch1(1) & trialSignal < trialTypeB_Range_Ch1(2));

%binarize the responses signal
trialSignalDigital = zeros(1,size(trialSignal,1));
%set A trials to 2
trialSignalDigital(idxTA) = 2;
%set B trials to 3
trialSignalDigital(idxTB) = 3;

%takes only the indices that go from low to high
diffTrialSignal = diff(trialSignalDigital);
%get indices and compensate for 1 point offset from diff output
%indices correspond to the reward location on index axis (non-restricted)
keepIdxOnsetTA = find(diffTrialSignal>1 & diffTrialSignal<3) + 1;
keepIdxOnsetTB = find(diffTrialSignal>2 & diffTrialSignal<4) + 1;

%put the onset signal indices into cells for each trial type
trialIdxOnsets{1} = keepIdxOnsetTA;
trialIdxOnsets{2} = keepIdxOnsetTB;

%add index filter to A and B trials here
%clean up reward location signals
removeTrialIdx = cell(1,size(trialIdxOnsets,2));

%remove trial type signals that are outside of full lap range
%get the time of id'd trial signals and remove those signals outside of
%complete lap range
for tt=1:size(trialIdxOnsets,2)
    %get the times
    trialIdxOnsetTime{tt} = Behavior.time(trialIdxOnsets{tt});
    trialRemove{tt} = find(trialIdxOnsetTime{tt} < Behavior.lap{1}(1) | trialIdxOnsetTime{tt} > Behavior.lap{end}(2));
end

%remove before full lap start or after full lap end signals if detected
for tt=1:size(trialRemove,2)
    %remove that signal from that trial subtype
    if ~isempty(trialRemove{tt})
        trialIdxOnsets{tt}(trialRemove{tt}) = [];
    end
end


%for both trial types
for tt=1:size(trialIdxOnsets,2)
    %for each index
    for rr=1:size(trialIdxOnsets{tt},2)
        %get neighboring voltages for the each index associated with the
        %texture
        neighborVoltagesTrial{tt}{rr} = trialSignal((trialIdxOnsets{tt}(rr)-2):(trialIdxOnsets{tt}(rr)+2));
        
        %are at least two voltage points within the range of the signal
        %if not, move to the other trials
        ptsAbove = numel(find(neighborVoltagesTrial{tt}{rr} >= trialTypeRange(tt,1) & neighborVoltagesTrial{tt}{rr} <= trialTypeRange(tt,2)));
        if ptsAbove < 2
            removeTrialIdx{tt} = [removeTrialIdx{tt}, rr];
        end
       
    end
end

%take unique indices
for tt=1:size(removeTrialIdx,2)
    removeTrialIdx{tt} = unique(removeTrialIdx{tt});
end

%copy the texture indices
trialIdxOnsetsClean = trialIdxOnsets;

%predefine removed onsets
removedOnsets = cell(1,size(trialIdxOnsetsClean,2));

%remove falsly identified trial signals
for tt=1:size(removeTrialIdx,2)
    if ~isempty(removeTrialIdx{tt})
        removedOnsets{tt} = trialIdxOnsetsClean{tt}(removeTrialIdx{tt});
        trialIdxOnsetsClean{tt}(removeTrialIdx{tt}) = [];
    end
end

%add the falsely identified trial signals to the other set of trials
%flip the cell onsets
removedOnsetsFlip = fliplr(removedOnsets);

%add the onsets to the opposite trials
for tt=1:size(trialIdxOnsetsClean,2)
    trialIdxOnsetsClean{tt} = [trialIdxOnsetsClean{tt}, removedOnsetsFlip{tt}];
end

%extract the trial type based on the voltage signal
%merge the two isolated indices together from code above
%merge and sort the indices into 1 matrix
%trialTypeMat = sort([keepIdxOnsetTA, keepIdxOnsetTB]);
trialTypeMat = sort([trialIdxOnsetsClean{1}, trialIdxOnsetsClean{2}]);

%located the indices of the respective A and B trials and set trial
%type in row below
[~,idxTypeMatA,~] = intersect(trialTypeMat(1,:),trialIdxOnsetsClean{1});
trialTypeMat(2,idxTypeMatA) = 2;
[~,idxTypeMatB,~] = intersect(trialTypeMat(1,:),trialIdxOnsetsClean{2});
trialTypeMat(2,idxTypeMatB) = 3;

%transpose the trialType matrix
%by column: indices, trial type, time, position
trialTypeMat = trialTypeMat';

%add associated raw time
trialTypeMat(:,3) = Behavior.time(trialTypeMat(:,1));
%add associated raw position
trialTypeMat(:,4) = Behavior.position(trialTypeMat(:,1));


%% reward delivery location

%raw CSV lick signal
rewardSignal = CSV(:,7);

%voltage range for lick on
rewardRange = [4.5, 4.7];

%find all indices where reward signal is high
idxReward = find(rewardSignal > rewardRange(1) & rewardSignal < rewardRange(2));

%digitize the lick signal
rewardSignalDigital = zeros(1,size(rewardSignal,1));
rewardSignalDigital(idxReward) = 1;

%indices of lick onset
%takes only the indices that go from low to high
diffRewardSignal = diff(rewardSignalDigital);
%get indices and compensate for 1 point offset from diff output
%indices correspond to the raw Lick signal (non-restricted)
keepIdxOnsetReward = find(diffRewardSignal>0) + 1;

%keepIdx = [1; find(diff(idx)>1)+1];

%keep only onset of each lick (count);
%idxFilter = idx(keepIdx);
idxFilterReward = keepIdxOnsetReward';
%reward position and time matrix concatenated laps

%find time associated with each lick (3rd column)
rewardMat(:,1) = Behavior.time(idxFilterReward);

%find bin associated with each lick onset
rewardMat(:,2) = Behavior.position(idxFilterReward);

%reward indices from CSV
rewardMat(:,3) = idxFilterReward;

%filter the reward delivered signals
%clean up reward location signals
removeRewOnIdx = [];
for rr=1:size(rewardMat,1)
    %get neighboring voltages for the each index associated with the
    %texture
    neighborVoltagesRewOn{rr} = rewardSignal((rewardMat(rr,3)-2):(rewardMat(rr,3)+2));
    
    %look at previous 2 voltage samples - should be less than 0.3V
    if sum((neighborVoltagesRewOn{rr}(1:2) < rewardRange(1))) ~=2
        removeRewOnIdx =  [removeRewOnIdx, rr];
    end
    %look at next 1 voltage samples should not be out of range
    if ((neighborVoltagesRewOn{rr}(end-1) < rewardRange(1)) || (neighborVoltagesRewOn{rr}(end-1) > rewardRange(2)))
        %add index for removal
        removeRewOnIdx =  [removeRewOnIdx,rr];
    end
    
end

%take unique indices
removeRewOnIdx = unique(removeRewOnIdx);

%copy the texture indices
rewardMatClean = rewardMat;

%remove false positive remward locations
if isempty(removeRewOnIdx) == 0
    rewardMatClean(removeRewOnIdx,:) = [];
end

%for two reward locations on the track
%find the minimum and maximum locations of delivered rewards
maxRewLoc = max(rewardMatClean(:,2));
minRewLoc = min(rewardMatClean(:,2));
%midway separation point of two rewards
splitLocPoint = (minRewLoc + maxRewLoc)/2;

%% figure out the trial order based on the reward signal

%laps already determined by behavior_lap_RZ script

%split each part of lick matrix into cell based on lap times
%(restricted /full laps)
lapNb = size(Behavior.lap,2);

%reshape Behavior.laps cell into matrix for search purposes
lapTimes = cell2mat(reshape(Behavior.lap,lapNb,[]));

%assign corresponding laps to DAC selected trials based on timing of
%signal from DAC in relation to lap start
%for each lap
for ii=1:size(lapTimes,1)
    %take difference from all lap start times and get index to which
    %the signal is closest
    [~,Im{ii}] = min(abs((lapTimes(:,1) - trialTypeMat(ii,3))));
    %assign lap correspondence to 5th columns
    trialTypeMat(ii,5) = Im{ii}(1);
end

%select only complete laps
%corresponding indices
selectLaps = find(trialTypeMat(:,5) >= 1 & trialTypeMat(:,5)<= lapNb);

%copy matrix
trialTypeMatFinal = trialTypeMat;

%select laps based on filter above
trialTypeMatFinal = trialTypeMatFinal(selectLaps,:);

%% separate the lick matrix into separate cells based on lapTimes matrix

%for each lap
for ll=1:size(lapTimes,1)
    %find licks within time range of each lap
    lickLapIdx = find(lickMat(:,3) >= lapTimes(ll,1) & lickMat(:,3) < lapTimes(ll,2));
    %find rewards within time range of each lap
    rewardLapIdx = find(rewardMatClean(:,1) >= lapTimes(ll,1) & rewardMatClean(:,1) < lapTimes(ll,2));
    %assign matrix to lap cell (licks)
    licks{ll} = lickMat(lickLapIdx,:);
    %assign matrix to lap cell (rewards)
    rewards{ll} = rewardMatClean(rewardLapIdx,:);
end

%bin edges
edges = (0:2:200);

%split licks into bins
for ll = 1:lapNb
    [N_licks(ll,:),~,~] = histcounts(licks{ll}(:,2),edges,'Normalization','count');
end

%trial order - starting from first lap
trialOrder = trialTypeMatFinal(:,2);

%assign to string names for plotting
for ii=1:size(trialOrder,1)
    if (trialOrder(ii) == 2)
        trialName{ii} = 'A';
    elseif (trialOrder(ii) == 3)
        trialName{ii} = 'B';
    end
end

%%%% - variables needed going forward
%rewardLocMatClean - RFID locations where a reward zone was activated
%detemine start of reward zone
%rewardMatClean - signal every time mouse collected a reward
%use to determine if any missed trials
%rewards - rewards that the animal collected broken down by lap
%trialTypeMatFinal - trial types according to trial type signal
%use to determine trial type regardless of performance of anumals
%licks - lick locations broken down by laps
%use to determine if animal performed correctly on given trial
%%%%

%% performance calculations

%determine which laps were correct/wrong
for ii = 1:lapNb
    %if A trial - check if licks in B ant or reward zone
    if trialOrder(ii) == 2
        if isempty(find(licks{ii}(:,2) >= binAntRange2(1) & licks{ii}(:,2) <= binRewardRange2(2)))
            trialCorrect(ii) = 1;
            trialCorrName{ii} = 'Y';
            
        else
            trialCorrect(ii) = 0;
            trialCorrName{ii} = 'N';
            trialOrder(ii) = 20;
        end
        %if B trial
    elseif trialOrder(ii) == 3
        if isempty(find(licks{ii}(:,2) >= binAntRange1(1) & licks{ii}(:,2) <= binRewardRange1(2)))
            trialCorrect(ii) = 1;
            trialCorrName{ii} = 'Y';
        else
            trialCorrect(ii) = 0;
            trialCorrName{ii} = 'N';
            trialOrder(ii) = 30;
        end
        
    end
end

%% Technical checks - work on this

%1) make sure reward signal comes on in each zone at each lap - reset
%signal on first cue occured;
%2) make sure that each lap is less than 200/205 cm - RFID tag cleared 
%3) make a matrix with time and position of each texture on each lap 


%% plot
%as single plot with trials increase along y axis

%raster plot indicating trial type and whether correct or not
figure('Position',[150 150 900 600]);
imagesc(N_licks);
hold on;
%trial type
text(105*ones(lapNb,1),1:lapNb,trialName);
%correct/wrong
tc = text(110*ones(lapNb,1),1:lapNb,trialCorrName);
%add red color to wrong trials

title('Lick map')
ylabel('Lap #');
xlabel('Spatial bin');
colormap('jet');
hold off

%% Save to structure
Behavior.performance.trialOrder = trialOrder;
Behavior.performance.trialCorrect = trialCorrect';

end

