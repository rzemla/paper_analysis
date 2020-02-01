function [excl_day_combined,excl_day_combined_day_nan] = generate_filter_LUTs(short_term_learn, short_term_recall, long_term_recall)

%% Get number of animals and sessions per experiment (use performance matrix)

%number of animals and sessions for each experiment
%stick to this order --> st_learn, st_recall, lt_recall

[nb_animal(1), nb_ses{1}] = return_nb_animal_ses(short_term_learn);
[nb_animal(2), nb_ses{2}] = return_nb_animal_ses(short_term_recall);
[nb_animal(3), nb_ses{3}] = return_nb_animal_ses(long_term_recall);

%% Generate session discard LUTs for each animal as part of struct
exclusion_vecs = cell(max(nb_animal),3);

%for each animal
%for each set of sessions

%for each experiment
%stick to this order --> st_learn, st_recall, lt_recall
for ee=1:3
    for aa=1:nb_animal(ee)
        %include all by default (1)
        exclusion_vecs{aa,ee} = ones(1,nb_ses{ee}(aa));
        %exclude the following from short term learning
        %if short term learning
        if ee==1
            if aa==1
                exclusion_vecs{aa,ee}([3,4]) = 0;
            elseif aa==3
                exclusion_vecs{aa,ee}([2]) = 0;
            elseif aa==5
                exclusion_vecs{aa,ee}([6,7]) = 0;
            end
            %long term recall
        elseif ee==3
            if aa==3
                exclusion_vecs{aa,ee}([2,3,4]) = 0;
            end
        end
    end
end
    
%% Generate session to day correspondence for each animal    
ses2day_match = cell(max(nb_animal),3);

%absolute assignment (without excluding low quality days)
%short term recall assignment
%animal, then experiment type
ses2day_match{1,1} = [1 2 3 4 6 7 8];
ses2day_match{2,1} = [1 2 3 4 6 7 8];
ses2day_match{3,1} = [1 2 3 4 5 6 7 8];
ses2day_match{4,1} = [1:9];
ses2day_match{5,1} = [1:6,8];
ses2day_match{6,1} = [1:6];

%short term recall
for aa=1:nb_animal(2)
    ses2day_match{aa,2} = [1 2 3 6 7 8 9];
end

%long term recall
for aa=1:nb_animal(3)
    ses2day_match{aa,3} = [1 6 16 20 25 30];
end

%% Combine select vectors with matching days into 1 cell
%top vector is which session to be excluded
%bottom vector is which day does the session correspond to
%preallocate
excl_day_combined = cell(max(nb_animal),3);

for ii = 1:size(excl_day_combined,1)
    for jj = 1:size(excl_day_combined,2)
        excl_day_combined{ii,jj} = [exclusion_vecs{ii,jj}; ses2day_match{ii,jj}];
    end
end


%% Set the excluded days to nan
%make copy
excl_day_combined_day_nan = excl_day_combined;

for ii = 1:size(excl_day_combined,1)
    for jj = 1:size(excl_day_combined,2)
        %if not an empty cell
        if ~isempty(excl_day_combined_day_nan{ii,jj})
            %find zeros in exclusion vector
            removeIdx_temp = excl_day_combined_day_nan{ii,jj}(1,:) == 0;
            excl_day_combined_day_nan{ii,jj}(2,removeIdx_temp) = nan;
        end
    end
end

end

