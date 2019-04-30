function [outputArg1,outputArg2] = visualize_matches(rows,cols,ROI_zooms,ROI_outlines)

%visualize matching ROIs across sessions

%matrix corresponding to display
display_matrix = reshape([1:(rows*cols)],cols,rows)';

%total number of matched ROIs to display
displayROInb = size(ROI_zooms,1);

%ROI input list (in chunks of 20 start to finish
start_end_ROI = [(1:20:20*floor(displayROInb/20))',(20:20:20*floor(displayROInb/20))'];
start_end_ROI = [start_end_ROI; [start_end_ROI(end)+1, displayROInb]];

%iterate through each range
for iter=1:size(start_end_ROI,1)
    figure('renderer','painters','Position', [2200 100 200 900]) %scale the width by 300 for each ses
    %assign the roi range for display
ROI = start_end_ROI(iter,1):start_end_ROI(iter,2);

%for each sessions
for rr=1:size(display_matrix,1)
    for ss=1:size(display_matrix,2)
        subplot(rows,cols,display_matrix(rr,ss))
        subaxis(rows, cols, display_matrix(rr,ss), 'SpacingHorizontal', 0.01,...
            'SpacingVertical',0.01,...
            'MarginLeft',0.05,'MarginRight',0.05,'MarginTop',0,'MarginBottom',0.1);
        %prevent out of bounds error
        if ~(rr> size(ROI,2))
            imagesc(ROI_zooms{ROI(rr) ,ss})
            hold on;
            colormap('gray')
            xticks([])
            yticks([])
            b= bwboundaries(ROI_outlines{ROI(rr),ss},'noholes');
            plot(b{1}(:,2),b{1}(:,1),'r')
            hold off
        end
    end
end

end



end

