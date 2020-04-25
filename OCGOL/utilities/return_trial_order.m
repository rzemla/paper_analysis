function [trial_type_vec] = return_trial_order(t_order,lap_nb)


%lap indices of A trials
A_lap_idx = find(t_order == 2 | t_order == 20);
%lap indices of B trials
B_lap_idx = find(t_order == 3 | t_order == 30);
%find punish laps
punish_lap_idx = find(t_order == -10);
%technical laps
tech_laps_idx = find(t_order == -1);

%combined indices into one cell with lap indices
%order: A,B,P,T (A trials, B trials, punish trials, technical trials)
lap_idx_split = {A_lap_idx,B_lap_idx,punish_lap_idx,tech_laps_idx};

%for each type of laps, assign numerical values
%A =2; B =3; punish = -10; technical = -1
%for each class of laps
%preallocate trial type vector
trial_type_vec = zeros(size(lap_nb,1),1);

for tt=1:4
    %check that cell is not empty (not laps of this type detected)
    if ~isempty(lap_idx_split{tt})
        %each selected lap index
        for ll=lap_idx_split{tt}'
            temp_idx = find( lap_nb == ll);
            %A
            if tt==1
                trial_type_vec(temp_idx) = 2;
                %B
            elseif tt==2
                trial_type_vec(temp_idx) = 3;
                %P
            elseif tt==3
                %T
                trial_type_vec(temp_idx) = -10;
            elseif tt==4
                trial_type_vec(temp_idx) = -1;
            end
        end
    end
end



end

