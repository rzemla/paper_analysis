function [nb_animal,nb_ses] = return_nb_animal_ses(input_struct)
%% Get number of animals and absolute sessions for each class of animals
%stick to this order --> st_learn, st_recall, lt_recall

%number of animals in experiment
nb_animal = size(input_struct.perf,2);

%for every animal, get number of experiments
for xx=1:nb_animal
    nb_ses(xx) = size(input_struct.perf{xx}.ses_perf,2);
end

end

