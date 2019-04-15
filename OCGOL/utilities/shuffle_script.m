
function [on_map_shuffle, shuffle_tuning_specificity]=shuffle_script(binary,Place_cell,Behavior,Events,options)
%% Import data

%ones where the animals is considered running across all restricted frames
run_ones=Behavior.run_ones;
%which bin each running frame corresponds to
bin=Place_cell.Bin;
%number of shuffles
Nshuffle=options.Nshuffle;
%which bin sizes to use
Nbin=options.Nbin;
%smoothing sigma for gaussian smoothing of rate maps
sigma=options.sigma_filter;

%% Create null distribution of onset

%binary is event onset across all running frames

%frames, frames, column vector of ROI indices
%p = randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n inclusive
%for each ROI permute the frames indices of
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1:size(binary,2))','un',0));

%ROIs becomes columns
binary_shuffle_idx=binary_shuffle_idx';

%for each ROI, shuffle the onsets to different frames (temporal shuffle?)
for i=1:size(binary,2)
    binary_shuffle(:,i)=binary(binary_shuffle_idx(:,i),i);
end

%% shuffle map
%run_binary_shuffle=binary_shuffle(run_ones==1,:);
run_binary_shuffle=binary_shuffle;

%for each bin range
for i=1:length(Nbin)
    %for each ROI
    for n=1:size(run_binary_shuffle,2)
        %bin = which position bin each frames belongs to
        %for each, ROI find the shuffles event onsets
        on_shuffle{i}{n}=bin{i}(run_binary_shuffle(:,n)==1);
        
        %for each bin, find how many events occured there
        for binN=1:Nbin(i)
            %map of number of events at each location for each ROI
            on_map_shuffle{i}(binN,n)=numel(find(on_shuffle{i}{n}==binN));
        end
    end
end

%% Shuffle tuning
%save null onset in structure to feed to shuffle function
Events_shuffle.options=Events.options;
%Events_shuffle.Run.run_onset_binary=binary_shuffle;
%temporal shuffle
Events_shuffle.Run.run_onset_binary=run_binary_shuffle;

%function tuning_specificity with null distribution
%feeds the shuffle events

[shuffle]=tuning_specificity_shuffle_RZ_V1(Place_cell,Behavior,Events_shuffle,options);

%tuning specificity values for each ROI (returns vector)
shuffle_tuning_specificity=shuffle.Tuning_Specificity.tuning_specificity;

end