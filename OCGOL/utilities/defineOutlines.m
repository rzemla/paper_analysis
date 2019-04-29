function [ROI_zooms,ROI_outlines] = defineOutlines(assign_all_matching,nbSes,session_vars)

%figure out what variables you need for this again - assemble into function


%turn this into GUI that allows exclusion of mismatched componenets

%plot zoomed in spatial component across 3 sessions
%plot outline of the id'd component
%plot the  associated trace 
%filter out neurons that are not matched

%concatenated components
compConcat = [];
%initialize flag variable for correction at edges of image
flagged = 0;

combineMatch = assign_all_matching;

%for each selected ROI
for jj=1:size(combineMatch,1)
    %assign the matched ROIs across sessions
    ROI = combineMatch(jj,:);
    
    %for plotting each (debug)
    %figure;
    for ss =1:nbSes
        %for plotting each (debug)
        %subplot(1,3,ss)
        %display template
        %imagesc(session_vars{ss}.template);
    
        %hold on
        %centers
        %scatter(session_vars{ss}.centers(ROI(ss),2),session_vars{ss}.centers(ROI(ss),1),'y*');
        %component outline
        %plot(session_vars{ss}.Coor_kp{ROI(ss)}(1,:),session_vars{ss}.Coor_kp{ROI(ss)}(2,:),'r', 'LineWidth',1);
        
        %zoom
        %xlim([session_vars{ss}.centers(ROI(ss),2)-60, session_vars{ss}.centers(ROI(ss),2)+60])
        %ylim([session_vars{ss}.centers(ROI(ss),1)-60, session_vars{ss}.centers(ROI(ss),1)+60])
        
        %get the rounded x and y range of the ROI of interest across
        %sessions
        yrange = [ceil(session_vars{ss}.centers(ROI(ss),2))-15, ceil(session_vars{ss}.centers(ROI(ss),2))+15];
        xrange = [ceil(session_vars{ss}.centers(ROI(ss),1))-15, ceil(session_vars{ss}.centers(ROI(ss),1))+15];
        
        %place each ROI zoom into a separate cell
        %if the area of interest if within range of the template
        if ~((sum(xrange < 1) > 0) || (sum(xrange > 512) > 0)) 
            if ~((sum(yrange < 1) > 0) || (sum(yrange > 512) > 0))
                disp(jj)
                ROI_zooms{jj,ss} = session_vars{ss}.template(xrange(1):xrange(2), yrange(1):yrange(2));
                %ROI outlines
                %create full-scale outline in a  sea of nans;
                ROI_outlines{jj,ss} = zeros(512,512);
                %create the binary outline
                %ROI_outlines{jj,ss}(session_vars{ss}.Coor_kp{ROI(ss)}(1,:)',session_vars{ss}.Coor_kp{ROI(ss)}(2,:)') = 1;
                ROI_outlines{jj,ss}(sub2ind([512,512],session_vars{ss}.Coor_kp{ROI(ss)}(2,:),session_vars{ss}.Coor_kp{ROI(ss)}(1,:))) =1;
                %clip the outline to match the field of of the ROI
                ROI_outlines{jj,ss} = ROI_outlines{jj,ss}(xrange(1):xrange(2), yrange(1):yrange(2));
                
            end
        end
        
        %adjust the xrange and yranges if beyond 1 or 512 
        %create blank 121x121 matrix and move the selected area into top
        %corner
        %determine which ranges are out of range
        xlow = xrange < 1;
        xhigh =  xrange > 512;
        ylow = yrange < 1;
        yhigh = yrange > 512;
        
        
        %set up conditional
        if xlow(1) ==1
            %xlowdiff = 
            xrange(1) = 1;
            flagged = 1;
        end
        if xlow(2) == 1
            xrange(2) = 1;
            flagged = 1;
        end
        if xhigh(1) ==1
            xrange(1) = 512;
            flagged = 1;
        end
        if xhigh(2) == 1
            xrange(2) = 512;
            flagged = 1;
        end
        if ylow(1) ==1
            yrange(1) = 1;
            flagged = 1;
        end
        if ylow(2) == 1
            yrange(2) = 1;
            flagged = 1;
        end
        if yhigh(1) ==1
            yrange(1) = 512;
            flagged = 1;
        end
        if yhigh(2) == 1
            yrange(2) = 512;
            flagged = 1;
        end
        
        %insert adjusted field
        if flagged == 1
            ROI_zooms{jj,ss} = zeros(15*2+1,15*2+1);
            ROI_zooms{jj,ss}(1:xrange(2)-xrange(1)+1, 1:yrange(2)-yrange(1)+1) = session_vars{ss}.template(xrange(1):xrange(2), yrange(1):yrange(2));
            %ROI outlines
            %create full-scale outline in a  sea of nans;
            %ROI_outlines{jj,ss} = nan(15*2+1,15*2+1);
            %create the binary outline
            ROI_outlines{jj,ss} = zeros(512,512);
            %create the binary outline
            %ROI_outlines{jj,ss}(session_vars{ss}.Coor_kp{ROI(ss)}(1,:)',session_vars{ss}.Coor_kp{ROI(ss)}(2,:)') = 1;
            ROI_outlines{jj,ss}(sub2ind([512,512],session_vars{ss}.Coor_kp{ROI(ss)}(2,:),session_vars{ss}.Coor_kp{ROI(ss)}(1,:))) =1;
            %store temp outline matrix
            temp = ROI_outlines{jj,ss};
            %set size with zeros
            ROI_outlines{jj,ss} = zeros(15*2+1,15*2+1);
            
            ROI_outlines{jj,ss}(1:xrange(2)-xrange(1)+1, 1:yrange(2)-yrange(1)+1) = temp(xrange(1):xrange(2), yrange(1):yrange(2));
        end
                   
                   
        %colormap
        %colormap( 'gray');
        %hold off
        
        %reset trim flag
        flagged = 0;
        
    end
end

end

