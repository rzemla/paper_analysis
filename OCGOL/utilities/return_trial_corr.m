function [trial_corr_vec] = return_trial_corr(t_corr,lap_nb)

%find correct laps
corr_lap_idx = find(t_corr == 1);
%find incorrect laps
incorr_lap_idx = find(t_corr == 0);
%find all else trials
all_else_lap_idx = find(t_corr == -1 | t_corr == -2);


%lap idx bundles by correct state
%corr, incorr, all else
corr_idx_split = {corr_lap_idx,incorr_lap_idx,all_else_lap_idx};

%preallocate correct vector
trial_corr_vec = zeros(size(lap_nb,1),1);

%correct, incorrect, all else
for tt=1:3
    %check that cell is not empty (not laps of this type detected)
    if ~isempty(corr_idx_split{tt})
        %each selected lap index
        for ll=corr_idx_split{tt}'
            temp_idx = find(lap_nb == ll);
            %corr
            if tt==1
                trial_corr_vec(temp_idx) = 1;
                %incorr = 0
            elseif tt==2
                trial_corr_vec(temp_idx) = 0;
                %all else = nan
            elseif tt==3
                %all else (punish or tech) = -1 (tech) or -2 (punish)
                trial_corr_vec(temp_idx) = nan;

            end
        end
    end
end


end

