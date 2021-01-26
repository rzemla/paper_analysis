function [place_structs] = extract_place_struct(path,range_animal,range_days,select_all)
%extracts and returns the existing place cell struct for give set of
%animals
%range_animal and range_days are vectors of all indices

%return data from all animals and sessions
if select_all == 1
    for ii = 1:numel(path)
        for jj = 1:numel(path{ii})
            disp([ii,jj])
            ca_analysis_out = fullfile(path{ii}{jj},'output','*.mat');
            ca_analysis_file = dir(ca_analysis_out);
            
            mat_file = fullfile(ca_analysis_file.folder,ca_analysis_file.name);
            
            place_structs(ii,jj) = load(mat_file, 'Place_cell');
            
        end
    end
else
    %return data from select animal on select sessions
    %use default values provided in function
    for ii = range_animal
        for jj = range_days
            disp([ii,jj])
            ca_analysis_out = fullfile(path{ii}{jj},'output','*.mat');
            ca_analysis_file = dir(ca_analysis_out);
            
            mat_file = fullfile(ca_analysis_file.folder,ca_analysis_file.name);
            
            place_structs(ii,jj) = load(mat_file, 'Place_cell');
            
        end
    end    
    
end




end
