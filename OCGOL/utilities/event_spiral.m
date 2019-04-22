function [spiralEvents] = event_spiral(animal_data, ROIrange, options)

%% Import variables

Behavior_split_lap = animal_data.Behavior_split_lap;
Events_split_lap = animal_data.Events_split_lap;
Behavior_split =  animal_data.Behavior_split;
Event_split = animal_data.Events_split;
Imaging_split = animal_data.Imaging_split;
Place_cell =animal_data.Place_cell;
Behavior_full =animal_data.Behavior;


%%

%V3 - add mean dF/F occupancy by lap to the plots

%all within run domain
position  = Behavior_split_lap.Run.position;
time = Behavior_split_lap.Run.time;
events_full = Events_split_lap.Run.run_onset_binary;
run_intervals = Behavior_split_lap.run_ones;

for ii=1:size(run_intervals,2)
    events{ii} = events_full{ii}(logical(run_intervals{ii}),:);
end

lap_times = Behavior_split_lap.lap;
 
%how many trial types are there
%set to 2 - A and B for now
trialTypes = 2;

%global trial type order across restricted laps
%first 20 trials for 20 - RZ CHECK LATER
trialOrder = Behavior_full.performance.trialOrder;

%find times from trials related to specific trial type
%A trials
trialTypeIdx{1} = find(trialOrder == 2);
%B trials
trialTypeIdx{2} = find(trialOrder == 3);

%create common placeFieldNb matrix
% for ii=1:size(Place_cell,2)
%     placeFieldNb(:,ii) = Place_cell{ii}.Field.placeFieldNb;
% end

%% Define the spiral parameters according to the number of laps
%equivalenet to number of laps shown
turns = size(trialOrder,1); %The number of turns the spiral will have (how many laps)

%x is the angle
x=[-1*pi*turns : 0.01 : pi*turns];
%r is that radius of the point
r=[0:1/(length(x)-1):1];

%scale to lap length
r_scaled = r.*turns;

%all parameters in the run frame domain
%find the frames index of event and position
%for trial type
for ii=1:size(trialTypeIdx,2)
    %for each lap belonging to that trial
    for ll= 1:size(trialTypeIdx{ii},1)
        %for each ROI
        for rr=1:size(events{trialTypeIdx{ii}(ll)},2)
            %event indices in run domain
            event_idx{ii}{ll}{rr} = find(events{trialTypeIdx{ii}(ll)}(:,rr) == 1);
            %position that corresponds to indices
            pos{ii}{ll}{rr} = position{trialTypeIdx{ii}(ll)}(event_idx{ii}{ll}{rr});
            %position vectors that will be used as input to spiral
            posVectors{ii}{ll}{rr} = trialTypeIdx{ii}(ll).*exp(1i.*((pos{ii}{ll}{rr}/200)*2*pi)).';
        end
    end
end

%make the nearest approximation to a point along the spiral vector defined above
%predefine to avoid empty cells at the end
valMin = posVectors;
idxMin = posVectors;
posVectorApprox = posVectors;
%for each ROI
%for each trial type
for kk = 1:size(trialTypeIdx,2)
    %for each lap belonging to that trial
    for ll = 1:size(pos{kk},2)
        %for each ROI
        for rr = 1:size(events{1},2)
            %for each event
            for ee=1:size(pos{kk}{ll}{rr},1)
                %fix empty cell clipping at endpoints
                [valMin{kk}{ll}{rr}(ee),idxMin{kk}{ll}{rr}(ee)] = min(abs( (r_scaled - ( (trialTypeIdx{kk}(ll)-1) + (pos{kk}{ll}{rr}(ee)/200) ) ) ));
                posVectorApprox{kk}{ll}{rr}(ee) = r_scaled(idxMin{kk}{ll}{rr}(ee))*exp(1i.*(pos{kk}{ll}{rr}(ee)/200)*2*pi);
            end
        end
    end
end

%% Collapse the data here - todo

%% Place field centers (test)

%centers of place fields in bins
% centerA = cell2mat(Place_cell{1}.placeField.center(ROIrange));
% centerB = cell2mat(Place_cell{2}.placeField.center(ROIrange));
% 
% %convert to radians
% centerA_angle = deg2rad((centerA/100)*360);
% centerB_angle = deg2rad((centerB/100)*360);


%% Plot spirals and dF/F maps

%Tuned = find(Place_cell{2}.Spatial_Info.significant_ROI==1);

for ii=1:size(ROIrange,2)
    
ROI = ROIrange(ii);
%in polar coordinate system without arrows
f = figure('rend','painters','pos',[0 0 1800 900]);

%turn off figure visibility
set(f, 'Visible', 'on');

%show index of frame
disp(ii);
disp('Displaying ROI:');
disp(ROI);

subplot(2,3,1)
polarplot(x,r_scaled,'k','Linewidth',1.5)
hold on

%plot A (2) trial events
for ll=1:size(idxMin{1},2)
polarscatter(angle(posVectorApprox{1}{ll}{ROI}),r_scaled(idxMin{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
%place field center
%polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
end

%plot B (3) trial events
for ll=1:size(idxMin{2},2)
polarscatter(angle(posVectorApprox{2}{ll}{ROI}),r_scaled(idxMin{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
%place field center
%polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
end
hold off

%dF/F maps by trial laps
%plot A
subplot(2,3,2)
imagesc(Place_cell{1}.dF_lap_map_ROI{ROI});
title('A trials','Color','b','FontSize', 14)
xlabel('Spatial bin')
ylabel('Lap')

%plot B trials
subplot(2,3,3)
imagesc(Place_cell{2}.dF_lap_map_ROI{ROI});
title('B trials','Color','r','FontSize', 14)
xlabel('Spatial bin')
ylabel('Lap')

subplots = get(f,'Children'); % Get each subplot in the figure

colormap( 'jet');
%reverse order (by most recent added)
for i=1:2%length(subplots) % for each subplot
    caxis(subplots(i),[0,2.5]); % set the clim
    
end

%annotation boxes
%ROI index
%what to display in the annotation box
labelStr = ['ROI: ', num2str(ROI)];
%add ROI labels
a = annotation(f,'textbox',[0.2 0.3 0 0 ],'String',labelStr,'FitBoxToText','on');
a.FontSize = 16;
a.Color = 'red';

%# of place fields
%labelStr = ['Place field #: ', num2str(placeFieldNb(ROI,:))];
%add ROI labels
a = annotation(f,'textbox',[0.2 0.25 0 0],'String',labelStr,'FitBoxToText','on');
a.FontSize = 16;
a.Color = 'blue';

%Spatial info
%labelStr = ['Spatial info: ', num2str([Place_cell{1}.Spatial_Info.significant_ROI(ROI), Place_cell{2}.Spatial_Info.significant_ROI(ROI)])];
%add ROI labels
a = annotation(f,'textbox',[0.2 0.2 0 0],'String',labelStr,'FitBoxToText','on');
a.FontSize = 16;
a.Color = 'black';



subplot(2,3,5)
hold on
ylim([-1.5 2.5])
stem(find(Event_split{1}.Run.run_onset_binary(:,ROI) ==1),ones(size(find(Event_split{1}.Run.run_onset_binary(:,ROI) ==1),1),1), 'b');
plot(Imaging_split{1}.trace_restricted(:,ROI),'k','Linewidth',0.5);
plot((Behavior_split{1}.resampled.position/200)-1.2,'Color', [0.4 0.4 0.4]);

hold off

subplot(2,3,6)
hold on
ylim([-1.5 2.5])
stem(find(Event_split{2}.Run.run_onset_binary(:,ROI) ==1),ones(size(find(Event_split{2}.Run.run_onset_binary(:,ROI) ==1),1),1), 'r');
plot(Imaging_split{2}.trace_restricted(:,ROI),'k','Linewidth',0.5);
plot((Behavior_split{2}.resampled.position/200)-1.2,'Color', [0.4 0.4 0.4]);

hold off



%%
    %get frames for movie
    mvFrame(ii) = getframe(f);
    
    if options.manualAdvance == 1
        pause;
    else
        %pause(options.plotSpeed);
    end
    
    if options.hold == 0
        hold off
    end
    
end

%% Make movie

if options.makeVideo == 1
    %movie playback of the frames generated above
    %mF = figure;
    %movie(mF,mvFrame,1,0.25);
    
    %save video
    %video writer object
    %v = VideoWriter(options.videoName,'Uncompressed AVI');
    %works with powerpoint with quality loss
    v = VideoWriter(options.videoName,'MPEG-4');
    
    %set frame rate
    v.FrameRate = 1.5;
    
    %open the video object
    open(v);
    
    disp('Making movie...')
    %write frames to video object
    writeVideo(v, mvFrame);
    
    %close video object;
    close(v);
    
    disp('Done');
    
end

%% Save to structure

spiralEvents.position = position;
%spiralEvents.lap = lap;
spiralEvents.event_idx = event_idx;
spiralEvents.idxMin = idxMin;
spiralEvents.positionVectorAppx = posVectorApprox;
spiralEvents.r_scaled = r_scaled;



end



