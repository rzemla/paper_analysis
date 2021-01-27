%% Load in all the animal directories

[path_dir_learn, crossdir_learn,...
            path_dir_st_recall, crossdir_st_recall,...
            path_dir_lt_recall, crossdir_lt_recall,...
            path_dir_task_sel] = return_animal_dirs();
        
        
%% Place all functions below into another function
%set path_dir directory here
%directories for task selective neurons in Fig. 2 and 3
path_dir = path_dir_task_sel;
%directories for short term_learning
%path_dir = path_dir_learn;
%directories for short term recall
%path_dir = path_dir_st_recall;
%directories for long term recall
%path_dir = path dir_lt_recall;

%set which trials should be used here
%correct trials 1 - A correct; 2 - B correct
%all trials 4 - all A; 5 - all B;
trial_type = [1,2];
%trial_type = [4,5];

%% Make copy of Place_cell_struct and save in archive folder

%create archive folder for each imaging session for each animal
%navigate to animal session folder
for ii=1:numel(path_dir)
    for jj=1:numel(path_dir{ii})
        cd(path_dir{ii}{jj})
        %check if archive folder already exisits
        if ~logical(exist(fullfile(path_dir{ii}{jj},'archive'),'dir'))
            %create folder if it does not exist
            mkdir(path_dir{ii}{jj},'archive')
        end
    end
end

%% Extract Place_cell structs for each session
%flag to select partcular struct or all structs
select_all = 1; % uses all animals and sessions

%select_all = 0; %uses subset defined by vectors below 
range_animal = 1;
range_days = 1;

[place_structs] = extract_place_struct(path_dir,range_animal,range_days,select_all);

%% Save in archive folder for each animal
for ii=1:numel(path_dir)
    for jj=1:numel(path_dir{ii})
        cd(fullfile(path_dir{ii}{jj},'archive'))
        %display animal and session index
        disp([ii, jj]);
        Place_cell = place_structs(ii, jj).Place_cell;
        %check if archived mat file already exists
        if isempty(dir('*.mat'))
            %append current date and time to the struct
            save_str = ['Place_cell_',datestr(now,'dd_mmm_yyyy_HH_MM_SS'),'.mat'];
            save(save_str,'Place_cell','-v7.3')
        else
            disp('Archived Place_cell struct already exists.');
        end
    end
end

%% Run the correction code from Jason 
%run Jason's correction code
%original SI tuned indices
%orig_ts_idx = find(Place_cell{1, 1}.Tuning_Specificity.significant_ROI ==1);

for ii =1:numel(path_dir)
    for jj=1:numel(path_dir{ii})
        disp([ii,jj])
        %make this parameter adjustable
        for tt=trial_type %for only correct A or B trials
            %input place struct; 1 - A corr, 2 - B corr, 4 - all A, 5 - all B
            place_struct_input = place_structs(ii, jj).Place_cell{1, tt};
            
            %original indices
            orig_si_idx = find(place_struct_input.Spatial_Info.significant_ROI ==1);
            %corrected indices of tuned ROIs
            out_temp = whichSignificant2(place_struct_input);
            if tt==1
                [animal(ii).sum_rowA(jj,:), animal(ii).si_dataA(jj)] = export_si_diff(orig_si_idx,out_temp,place_struct_input);
            elseif tt==2
                [animal(ii).sum_rowB(jj,:), animal(ii).si_dataB(jj)] = export_si_diff(orig_si_idx,out_temp,place_struct_input);
            end
        end
        
    end
end

%% Generate new Place_struct with updated variables - CONTINUE HERE
%replace p-values and significant ROIs
%place the sig ROIs, p-values 
%(8th row - 100 bins) - check how it compare to max score obtained from
%Jason - not the same b/c these are adjusted compared to max

for ii =1%:numel(path_dir)
    %check that this works for the other animal/session combos (should
    %work)
    for jj=1%:numel(path_dir{ii})
        %for each trial type
        for tt=trial_type
            %table of spatial info for each spatial bin size
            si_info_table_ori_temp = place_structs(ii,jj).Place_cell{1, tt}.Spatial_Info.Spatial_Info ;
            %original p-values
            si_pval_ori_temp = place_structs(1).Place_cell{1, 1}.Spatial_Info.ROI_pvalue;
            
        end
    end
end





place_structs(1).Place_cell{1, 1}.Spatial_Info.significant_ROI

place_structs(1).Place_cell{1, 1}.Tuned_ROI

place_structs(1).Place_cell{1,1}.Tuned_ROI_mak

%make a copy within the struct itself in archive

%recalculate shared ROIs



%% Overwrite exisiting Place_cell struct


