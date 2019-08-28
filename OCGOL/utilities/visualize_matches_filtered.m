function [outputArg1,outputArg2] = visualize_matches_filtered(rows,cols,registered,ROI_zooms,ROI_outlines,nb_ses)

%number of session to look at
%nb_ses = 6;

%% Match list (filtered)

%matching list across all sessions (post manual filtering)
matching_list = registered.multi.assigned_filtered;

%remove no match entries
remove_no_match_idx = find(sum(isnan(matching_list),2) == 0);
%copy match list
match_list_clean = matching_list;
%clear all nan entiries
match_list_clean(remove_no_match_idx,:) = [];

%sort matching list according to to number of matches
[match_sort_nb, match_sort_idx] = sort(sum(~isnan(match_list_clean),2), 'descend');

%get resorted match list
match_list_sorted = match_list_clean(match_sort_idx,:);

%%  Visualize matching ROIs across sessions

%matrix corresponding to display for each set of 20
display_matrix = reshape([1:(rows*cols)],cols,rows)';

%extend display matrix to show 100 ROIs at a time
%batchs of 20 (for 100 ROIs)
batch_nb = 5; %100 ROIs per cluster of subplots
display_matrix_ex = reshape([1:(batch_nb*rows*cols)],batch_nb*cols,rows)';
%reshape to set of 20 per session column (20x7, 20x7 ....) for 100 total
%ROIs
% for bb=1:batch_nb
%     display_matrix_celled{bb} = display_matrix_ex(1+((bb-1)*20):20+((bb-1)*20),:);
% end
% %transform to continuous matrix
% display_matrix_ex = cell2mat(display_matrix_celled);

%total number of matched ROIs to display (size of matching list)
displayROInb = size(match_list_sorted ,1);

%ROI input list (in chunks of 20 start to finish
start_end_ROI = [(1:20:20*floor(displayROInb/20))',(20:20:20*floor(displayROInb/20))'];
start_end_ROI = [start_end_ROI; [start_end_ROI(end)+1, displayROInb]];

%break into groups of 5
nb_plots = ceil(size(start_end_ROI,1)/batch_nb);

%ROI lists into batches (cell) - batch of ROI ranges to plot together
start_batch = 1:batch_nb:size(start_end_ROI);
end_batch = [start_batch(2:end)-1 size(start_end_ROI,1)];
comb_batch = [start_batch; end_batch];

%number of batch plots
for pp=1:nb_plots
    
    figure('renderer','painters','Position', [2200 100 1900 900]) %scale the width by 300 for each ses
    %assign the roi range for display (group of ROI list
    %iterate through each range
    start_batch_idx = 0;
    for iter=comb_batch(1,pp):comb_batch(2,pp)
        %batch counter
        
        ROI = start_end_ROI(iter,1):start_end_ROI(iter,2);
        
        for rr=1:size(display_matrix_ex,1)
            %for each session
            for ss=1:nb_ses %size(display_matrix,2)
                %convert to position on extended plot
                %subplot(rows,batch_nb*cols,display_matrix_ex(rr,ss+(start_batch_idx*6)))
                disp(display_matrix_ex(rr,ss+(start_batch_idx*6)))
                %determine spacing between columns and rows
                %vertically space every batch
                %             if ss+(start_batch_idx*6) == 6 || ss+(start_batch_idx*6) == 12
                %                 subaxis(rows, batch_nb*cols, display_matrix_ex(rr,ss+(start_batch_idx*6)), 'SpacingHorizontal', 0.001,...
                %                     'SpacingVertical',0.001,...
                %                     'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0,'MarginBottom',0.1)
                %             else
                subaxis(rows, batch_nb*cols, display_matrix_ex(rr,ss+(start_batch_idx*6)), 'SpacingHorizontal', 0.0015,...
                    'SpacingVertical',0.001,...
                    'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.1);
                % end
                %prevent out of bounds error
                if ~(rr> size(ROI,2))
                    %if no match
                    if ~isnan(match_list_sorted(ROI(rr),ss))
                        imagesc(ROI_zooms{ss}{match_list_sorted(ROI(rr),ss)})
                    else %do nothing or make black in future
                        %set background axis color to black
                        set(gca,'color',0*[1 1 1]);
                    end
                    hold on;
                    
                    colormap('gray')
                    xticks([])
                    yticks([])
                    %if empty, do not draw boundaries
                    if ~isnan(match_list_sorted(ROI(rr),ss))
                        b= bwboundaries(ROI_outlines{ss}{match_list_sorted(ROI(rr),ss)},'noholes');
                        plot(b{1}(:,2),b{1}(:,1),'r')
                    else
                    end
                    hold off
                end
            end
        end
        
        %advance batch index
        start_batch_idx = start_batch_idx + 1;
    end
    
end

end

