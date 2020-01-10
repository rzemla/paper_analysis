function [outputArg1,outputArg2] = raster_spiral_partial_maps(session_vars,CNMF_vars,removeROI,templates,partial_idx_by_animal_zone,options)

%% Define/load variables for each session

%for each session
for ii = 1:size(session_vars,2)
    % behavior and imaging related variables
    Behavior_split_lap{ii} = session_vars{ii}.Behavior_split_lap;
    Events_split_lap{ii} = session_vars{ii}.Events_split_lap;
    Behavior_split{ii} = session_vars{ii}.Behavior_split;
    Event_split{ii} = session_vars{ii}.Events_split;
    Imaging_split{ii} = session_vars{ii}.Imaging_split;
    Place_cell{ii} = session_vars{ii}.Place_cell;
    Behavior_full{ii} = session_vars{ii}.Behavior;
    
    %all within run domain
    position{ii}  = Behavior_split_lap{ii}.Run.position;
    time{ii} = Behavior_split_lap{ii}.Run.time;
    events_full{ii} = Events_split_lap{ii}.Run.run_onset_binary;
    run_intervals{ii} = Behavior_split_lap{ii}.run_ones;
    
    %global trial type order across restricted laps
    trialOrder{ii} = Behavior_full{ii}.performance.trialOrder;
end

%for each session
for ss=1:size(session_vars,2)
    %for each lap
    for ii=1:size(run_intervals{ss},2)
        events{ss}{ii} = events_full{ss}{ii}(logical(run_intervals{ss}{ii}),:);
    end
end

%lap_times = Behavior_split_lap.lap;
%how many trial types are there
%set to 2 - A and B for now
%trialTypes = 2;

%make this a condition for the type of trials that are compared
%all correct; all regardless of correct

%for each session (all A or B regardless of correct)
for ss = 1:size(session_vars,2)
    %find times from trials related to specific trial type
    %A trials
    trialTypeIdx{ss}{1} = find(trialOrder{ss} == 2 | trialOrder{ss} == 20);
    %B trials
    trialTypeIdx{ss}{2} = find(trialOrder{ss} == 3 | trialOrder{ss} == 30);
end

%get the chosen spatial components of each ROI (A matrix from CNMF output)
A_keep_sel = CNMF_vars{1, 1}.A_keep(:,removeROI{1}.compSelect);  

%center of mass of each component
A_keep_sel_com = com(A_keep_sel,512,512);

%outline
coor_keep_sel = CNMF_vars{1}.Coor_kp(removeROI{1}.compSelect);


%% Define the spiral parameters according to the number of laps
%equivalent to number of laps for each session (turns = laps)
for ss = 1:size(session_vars,2)
    turns(ss) = size(trialOrder{ss},1); %The number of turns the spiral will have (how many laps)
    
    %x is the angle
    x{ss} = [-1*pi*turns(ss) : 0.01 : pi*turns(ss)];
    
    %r is that radius of the point
    r{ss} = [0:1/(length(x{ss})-1):1];
    
    %scale to lap length
    r_scaled{ss} = r{ss}.*turns(ss);
end

%all parameters in the run frame domain
%find the frames index of event and position
%for each session
for ss=1:size(session_vars,2)
    %for trial type (A or B)
    for ii=1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll= 1:size(trialTypeIdx{ss}{ii},1)
            %for each ROI
            for rr=1:size(events{ss}{trialTypeIdx{ss}{ii}(ll)},2)
                %event indices in run domain
                event_idx{ss}{ii}{ll}{rr} = find(events{ss}{trialTypeIdx{ss}{ii}(ll)}(:,rr) == 1);
                %position that corresponds to indices
                pos{ss}{ii}{ll}{rr} = position{ss}{trialTypeIdx{ss}{ii}(ll)}(event_idx{ss}{ii}{ll}{rr});
                %position vectors that will be used as input to spiral
                posVectors{ss}{ii}{ll}{rr} = trialTypeIdx{ss}{ii}(ll).*exp(1i.*((pos{ss}{ii}{ll}{rr}/200)*2*pi)).';
            end
        end
    end
end

%make the nearest approximation to a point along the spiral vector defined above
%predefine to avoid empty cells at the end
valMin = posVectors;
idxMin = posVectors;
posVectorApprox = posVectors;

%foreach session
for ss =1:size(session_vars,2)
    %for each trial type
    for kk = 1:size(trialTypeIdx{ss},2)
        %for each lap belonging to that trial
        for ll = 1:size(pos{ss}{kk},2)
            %for each ROI
            for rr = 1:size(events{ss}{1},2)
                %for each event
                for ee=1:size(pos{ss}{kk}{ll}{rr},1)
                    %fix empty cell clipping at endpoints
                    [valMin{ss}{kk}{ll}{rr}(ee),idxMin{ss}{kk}{ll}{rr}(ee)] = min(abs( (r_scaled{ss} - ( (trialTypeIdx{ss}{kk}(ll)-1) + (pos{ss}{kk}{ll}{rr}(ee)/200) ) ) ));
                    posVectorApprox{ss}{kk}{ll}{rr}(ee) = r_scaled{ss}(idxMin{ss}{kk}{ll}{rr}(ee))*exp(1i.*(pos{ss}{kk}{ll}{rr}(ee)/200)*2*pi);
                end
            end
        end
    end
end

%% Split calcium traces into intervals to avoid joining by plot into continuous plot

%breakpoints
diff_A_trials_idx = find(diff(Imaging_split{1}{4}.time_restricted)> 1);
diff_B_trials_idx = find(diff(Imaging_split{1}{5}.time_restricted)> 1);

%construct start and end idx matrices - correct
start_idx.A = [1; diff_A_trials_idx+1];
end_idx.A = [diff_A_trials_idx; size(Imaging_split{1}{4}.time_restricted,1)];
%combined start and end idx
start_end_idx.A = [start_idx.A, end_idx.A]; 
%for B trials
start_idx.B = [1; diff_B_trials_idx+1];
end_idx.B = [diff_B_trials_idx; size(Imaging_split{1}{5}.time_restricted,1)];
%combined start and end idx
start_end_idx.B = [start_idx.B, end_idx.B]; 

%% Get lap indices for each lap in all B or B trials

%get unique lap indices
lapA_idxs = unique(Behavior_split{1}{4}.resampled.lapNb);
lapB_idxs = unique(Behavior_split{1}{5}.resampled.lapNb);

%get lap start and end indices for all A or B trials
%all A
for ll=1:size(lapA_idxs,1)
    lap_idxs.A(ll,1) = find(Behavior_split{1}{4}.resampled.lapNb == lapA_idxs(ll),1,'first');
    lap_idxs.A(ll,2) = find(Behavior_split{1}{4}.resampled.lapNb == lapA_idxs(ll),1,'last');
end

%all B
for ll=1:size(lapB_idxs,1)
    lap_idxs.B(ll,1) = find(Behavior_split{1}{5}.resampled.lapNb == lapB_idxs(ll),1,'first');
    lap_idxs.B(ll,2) = find(Behavior_split{1}{5}.resampled.lapNb == lapB_idxs(ll),1,'last');
end

%% Event onsets in run interval
%for each ROI
for rr=1:size(events{ss}{1},2)
    %time of significant run events in A
    event_norm_time.A{rr} =Imaging_split{1}{4}.time_restricted(find(Event_split{1}{4}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in A
    event_norm_pos_run.A{rr} = Behavior_split{1}{4}.resampled.position_norm(find(Event_split{1}{4}.Run.run_onset_binary(:,rr) == 1));
    %time of significant run events in B
    event_norm_time.B{rr} =Imaging_split{1}{5}.time_restricted(find(Event_split{1}{5}.Run.run_onset_binary(:,rr) == 1))/60;
    %normalizesd position of significant run events in B
    event_norm_pos_run.B{rr} = Behavior_split{1}{5}.resampled.position_norm(find(Event_split{1}{5}.Run.run_onset_binary(:,rr) == 1));
end

%% Plot partial remapping traces 

%get STCs for both A and B trials
figure
Aselect = [402 479 349 346 91 134 387 562];

subplot(6,2,ii)

%% Gaussian smooth non_normalized STCs (visualization purposes)
options.sigma_filter = 3;
gaussFilter = define_Gaussian_kernel(options);

nb_ROI = size(session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm,2);

%expand to do circular convolution
non_norm_STC.A = session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm';
non_norm_STC.B = session_vars{1}.Place_cell{2}.Spatial_tuning_curve_no_norm';

%expand the edges
non_norm_STC_exp.A = [non_norm_STC.A(:,51:100),non_norm_STC.A,non_norm_STC.A(:,1:50)]';
non_norm_STC_exp.B = [non_norm_STC.B(:,51:100),non_norm_STC.B,non_norm_STC.B(:,1:50)]';

%for all ROIs
for rr=1:nb_ROI
    %smooth all A trial sessions
    A_smooth(:,rr) = conv(session_vars{1}.Place_cell{1}.Spatial_tuning_curve_no_norm(:,rr),gaussFilter, 'same');
    %circular convolution
    A_smooth_exp(:,rr) = conv(non_norm_STC_exp.A(:,rr),gaussFilter, 'same');
    
    %smooth all B trial sessions
    B_smooth(:,rr) = conv(session_vars{1}.Place_cell{2}.Spatial_tuning_curve_no_norm(:,rr),gaussFilter, 'same');
    %circular convolution
    B_smooth_exp(:,rr) = conv(non_norm_STC_exp.B(:,rr),gaussFilter, 'same');
    
end

%cicularly smoothed
A_smooth_circ = A_smooth_exp(51:150,:);
B_smooth_circ = B_smooth_exp(51:150,:);

%single view match
figure
hold on
%zone II A (far)
subplot(2,1,1)
hold on
plot(A_smooth_circ(:,346),'b')
plot(B_smooth(:,346),'r')

subplot(2,1,2)
hold on
plot(A_smooth(:,205),'b')
plot(B_smooth(:,205),'r')

% normalize from 0-1
smooth_maps_combined = [A_smooth' B_smooth'];

min_val_STC = min(smooth_maps_combined,[],2);
max_val_STC = max(smooth_maps_combined,[],2);

%normalize A/B matrix
A_smooth_norm = (A_smooth - min_val_STC')./(max_val_STC - min_val_STC)';
B_smooth_norm = (B_smooth - min_val_STC')./(max_val_STC - min_val_STC)';


% normalize from 0-1 (circularly
smooth_maps_combined_circ = [A_smooth_circ' B_smooth_circ'];

min_val_STC_circ = min(smooth_maps_combined_circ,[],2);
max_val_STC_circ = max(smooth_maps_combined_circ,[],2);

%normalize A/B matrix
A_smooth_norm_circ = (A_smooth_circ - min_val_STC_circ')./(max_val_STC_circ - min_val_STC_circ)';
B_smooth_norm_circ = (B_smooth_circ - min_val_STC_circ')./(max_val_STC_circ - min_val_STC_circ)';

figure
hold on
%zone II A (far)
subplot(2,1,1)
hold on
plot(A_smooth_norm(:,346),'b')
plot(B_smooth_norm(:,346),'r')


%zone II B (near)
figure
hold on
plot(session_vars{1, 1}.Place_cell{1}.Spatial_tuning_curve(:,205),'b')
plot(session_vars{1, 1}.Place_cell{2}.Spatial_tuning_curve(:,205),'r')

Asel = [91 134 387 562];
Bsel = [166 332 421];

for rr=Bsel
figure
hold on
title(num2str(rr))
%zone III A (far)
plot(session_vars{1, 1}.Place_cell{1}.Spatial_tuning_curve(:,rr),'b')
plot(session_vars{1, 1}.Place_cell{2}.Spatial_tuning_curve(:,rr),'r')
end

%zone III B (near)
figure
hold on
plot(session_vars{1, 1}.Place_cell{1}.Spatial_tuning_curve(:,205),'b')
plot(session_vars{1, 1}.Place_cell{2}.Spatial_tuning_curve(:,205),'r')
%zone III

figure
hold on
plot(session_vars{1, 1}.Place_cell{1}.Spatial_Info.rate_map_sm_plot{8}(:,205),'b' )
plot(session_vars{1, 1}.Place_cell{2}.Spatial_Info.rate_map_sm_plot{8}(:,205),'r' )

%% Plot as one big plot
marker_size = 12;

%get colors
paper_colors = return_paper_colormap;

for kk=1
    if options.partial == 1
                
        %number of each type of ROI - A,B,common
        % zone I (first 2 ROIs)
        %ROI_all = [309,326,602];
        % zone II (first 2 ROIs)
        %ROI_all = [346, 147, 147];
      
        % zone III (first 2 ROIs)
        %ROI_all =[349,300,346,205,91,166];
        %animal #8 (zone I, zone II, zone III)
        ROI_all =[349,300,346,147,387,166];
        
        %6x7
        linear_idx = 1:42;
        matrix_order = reshape(linear_idx,6,7)';
        
        ROI_colors = [139,0,139]./255;
        ROI_outline_order = matrix_order(:,1);
        
        figure('Position',[1920 40 1920 960],'Renderer','painters')
        %Show spatial outline of the component
        for rr=1:6
            subplot(7,6,ROI_outline_order(rr))
            imagesc(templates{1}.template);
            hold on
            title(num2str(ROI_all(rr)))
            axes(gca);
            axis square
            xticks(gca,[])
            yticks(gca,[])
            
            grayMap = brighten(gray,0.5);
            colormap(gca,grayMap)
            %plot componenet outline
            plot(coor_keep_sel{ROI_all(rr)}(1,:),coor_keep_sel{ROI_all(rr)}(2,:),'Color',ROI_colors, 'LineWidth',1);
            %zoom into component based on center of mass calculation
            xlim([A_keep_sel_com(ROI_all(rr),2)-30, A_keep_sel_com(ROI_all(rr),2)+30])
            ylim([A_keep_sel_com(ROI_all(rr),1)-30, A_keep_sel_com(ROI_all(rr),1)+30])
        end
        
        %show color-coded (by trial) calcium traces of the component
        ROI_trace_order = matrix_order(1:6,2:4);
        for rr=1:6
            subplot(7,6,ROI_trace_order(rr,:))
            hold on
            xlim([0 12.5])
            ylim([-0.2 5])
            yticks([0 1 2])
            ylabel('dF/F');
            xticks(0:3:12.5);
            set(gca,'FontSize',14)
            set(gca,'LineWidth',1)
            %convert to minuntes and offset to start at 0 min (0s)
            %for each split (cell above) % for A trial traces
            for ii=1:size(start_end_idx.A,1)
                plot(Imaging_split{1}{4}.time_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2))/60,...
                    Imaging_split{1}{4}.trace_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2),ROI_all(rr)),'Color',paper_colors(1,:),...
                    'LineWidth',1)
            end
            %for each split (cell above) % for B trial traces
            for ii=1:size(start_end_idx.B,1)
                plot(Imaging_split{1}{5}.time_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2))/60,...
                    Imaging_split{1}{5}.trace_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2),ROI_all(rr)),'Color',paper_colors(2,:),...
                    'LineWidth',1)
            end
        end
        
        %polar plots
        ROI_polar_order = matrix_order(1:6,5);
        for rr=1:6
            %tuning vectors of the component
            pax1 = subplot(7,6,ROI_polar_order(rr),polaraxes);
            hold on
            pax1.FontSize = 14;
            angles = 0;
            pax1.ThetaTick = angles;
            thetaticklabels(pax1,{'lap start'});
            rticks(pax1, []);
            rlim(pax1,[0 25]);
            
            polarplot(x{1},r_scaled{1},'k','Linewidth',1.0)
            
            %plot A (2) trial events
            for ll=1:size(idxMin{1}{1},2)
                polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{1}{ll}{ROI_all(rr)}),marker_size,'o','MarkerFaceColor',paper_colors(1,:),'MarkerEdgeColor',paper_colors(1,:))
                %place field center
                %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
            end
            %plot B (3) trial events
            for ll=1:size(idxMin{1}{2},2)
                polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{2}{ll}{ROI_all(rr)}),marker_size,'o','MarkerFaceColor',paper_colors(2,:),'MarkerEdgeColor',paper_colors(2,:))
                %place field center
                %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
            end
        end
        
        ROI_STC_order = matrix_order(1:6,6);
        %plot STCs
        for rr=1:6
            subplot(7,6,ROI_STC_order(rr))
            hold on
            plot(A_smooth_norm_circ(:,ROI_all(rr)),'Color',paper_colors(1,:),'LineWidth',2)
            plot(B_smooth_norm_circ(:,ROI_all(rr)),'Color',paper_colors(2,:),'LineWidth',2)
            xticks([0 100])
            xticklabels({'0','1'})
            yticks([0 1])
            %set(gca,'FontSize',10)
            if rr==6
            xlabel('Normalized position')
            ylabel('Normalized activity')
            end
        end
        
        %plot normalized position
        subplot(7,6,matrix_order(7,2:4))
        hold on
        yticks([0 0.5 1])
        ylabel('Normalized position')
        xlabel('Time [min]');
        xticks(0:3:12.5);
        xlim([0.1 12.5])
        ylim([0 1])
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %A laps
        for ii=1:size(lap_idxs.A,1)
            plot(Imaging_split{1}{4}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
                Behavior_split{1}{4}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
                'Color',paper_colors(1,:),'LineWidth',1.5)
        end
        %B laps
        for ii=1:size(lap_idxs.B,1)
            plot(Imaging_split{1}{5}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
                Behavior_split{1}{5}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
                'Color',paper_colors(2,:),'LineWidth',1.5)
        end
    end
end

%% Plot raster, event spiral and matching ROIs from FOV
%21 partial remap
%AB_list =[2	5	6	8	11	14	15	16	18	19	21	23	24	30	33	35	36	38	39	41	42	43	46	47	49	50	56	57	58	59	60	61	62	63	64	69	70	71	77	78	80	82	84	86	87	88	91	92	93	95	104	107	110	112	114	115	116	117	118	119	122	125	127	128	130	132	133	134	135	136	137	139	142	143	148	154	155	162	164	165	166	168	169	175	176	177	179	180	181	186	187	188	189	190	191	194	195	197	198	199	201	207	208	213	214	216	217	218	223	226	227	228	231	233	234	235	237	241	242	244	245	246	247	254	256	260	263	264	268	269	270	271	272	274	279	280	282	283	284	285	286	288	291	292	296	298	299	304	306	312	318	320	321	326	329	333	334	343	348	350	357	361	362	363	365	366	367	368	369	370	372	373	374	375	379	380	381	382	384	386	387	389	391	392	393	394	398	399	400	401	402	405	407	408	410	411	412	413	415	416	417	419	423	427	429	431	432	433	436	437	438	440	441	443	447	448	449	450	451	453	457	458	460	461	462	463	466	469	470	471	472	474	475	476	477	478	480	481	482	483	486	489	490	491	492	493	494	495	498	500	506	510	512	513	514	515	519	520	521	522	523	524	525	526	527	529	531	534	535	538	540	541	542	543	544	545	549	554	555	556	557	562	563	564	565	566	568	570	571	574	575	576	578	580	581	583	589	590	591	592	593	594	596	597	598	600	603	604	605	606	607	610	612	614	616	620	621	624	627	629	633	635	639	640	641	642	643	649	651	652	653	655	657	660	665	668	669	672	676	682	687	690	692	701];
if 0
    %ROI = 38;
    %split cells
    %365,371, 452, 371
    %for ROI=AB_list
    for ROI=45:size(session_vars{1}.Place_cell{1}.Spatial_tuning_curve,2)
        figure('Position',[1920 40 1920 960])
        %Show spatial outline of the component
        subplot(3,5,1)
        imagesc(templates{1}.template);
        hold on
        title(num2str(ROI))
        axes(gca);
        axis square
        xticks(gca,[])
        yticks(gca,[])
        grayMap = brighten(gray,0.5);
        colormap(gca,grayMap)
        %plot componenet outline
        plot(coor_keep_sel{ROI}(1,:),coor_keep_sel{ROI}(2,:),'m', 'LineWidth',1);
        %zoom into component based on center of mass calculation
        xlim([A_keep_sel_com(ROI,2)-30, A_keep_sel_com(ROI,2)+30])
        ylim([A_keep_sel_com(ROI,1)-30, A_keep_sel_com(ROI,1)+30])
        
        %show color-coded (by trial) calcium traces of the component
        subplot(3,5,[2 3 4])
        hold on
        ylim([-0.2 2.5])
        yticks([0 1 2])
        ylabel('dF/F');
        xticks(0:3:12);
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %convert to minuntes and offset to start at 0 min (0s)
        %for each split (cell above) % for A trial traces
        for ii=1:size(start_end_idx.A,1)
            plot(Imaging_split{1}{4}.time_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2))/60,...
                Imaging_split{1}{4}.trace_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2),ROI),'Color',[0 0 1 1],...
                'LineWidth',1)
        end
        %for each split (cell above) % for B trial traces
        for ii=1:size(start_end_idx.B,1)
            plot(Imaging_split{1}{5}.time_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2))/60,...
                Imaging_split{1}{5}.trace_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2),ROI),'Color',[1 0 0 1],...
                'LineWidth',1)
        end
        
        %tuning vectors of the component
        pax1 = subplot(3,5,5,polaraxes);
        hold on
        pax1.FontSize = 14;
        angles = 0;
        pax1.ThetaTick = angles;
        thetaticklabels(pax1,{'lap start'});
        rticks(pax1, []);
        rlim(pax1,[0 30]);
        
        polarplot(x{1},r_scaled{1},'k','Linewidth',1.5)
        
        %plot A (2) trial events
        for ll=1:size(idxMin{1}{1},2)
            polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI}),r_scaled{1}(idxMin{1}{1}{ll}{ROI}),marker_size,'bo','MarkerFaceColor','b')
            %place field center
            %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
        end
        %plot B (3) trial events
        for ll=1:size(idxMin{1}{2},2)
            polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI}),r_scaled{1}(idxMin{1}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
            %place field center
            %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
        end
        
        %spiral map of the componenet
        subplot(3,5,10)
        %correct A
        compass(Place_cell{1}{1}.Tuning_Specificity.tuning_vector{ROI},'b');
        hold on
        %correct B
        compass(Place_cell{1}{2}.Tuning_Specificity.tuning_vector{ROI},'r');
        %deal with labels
        labels = findall(gca,'type','text');
        set(labels,'visible','off');
        %set(findall(gcf, 'String', '30', '-or','String','60') ,'String', '  ');
        
        %only correct
        %only works if using polar axis; plot in cartesian ref frame when using
        %compass
        %polarplot([0+0i,Place_cell{1}{1}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
        
        %plot normalized position
        subplot(3,5,[7 8 9])
        hold on
        yticks([0 0.5 1])
        ylabel('Normalized position')
        xlabel('Time [min]');
        xticks(0:3:12);
        ylim([0 1])
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %A laps
        for ii=1:size(lap_idxs.A,1)
            plot(Imaging_split{1}{4}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
                Behavior_split{1}{4}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
                'Color',[0 0 1 0.6],'LineWidth',1.5)
        end
        %B laps
        for ii=1:size(lap_idxs.B,1)
            plot(Imaging_split{1}{5}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
                Behavior_split{1}{5}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
                'Color',[1 0 0 0.6],'LineWidth',1.5)
        end
        %overlay significant calcium run events
        %A
        scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
        %B
        scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
        pause;
        clf;
    end
end

%% Plot figure 2b traces 
if options.plotFigure2 ==1
    
    %number of each type of ROI - A,B,common
    %ROI_all = [234,246,197] 
    %ROI_all = [398	404	455];
    %ROI_all = [111, 246, 197];
    ROI_colors = {'b','r','m'};
    ROI_outline_order = [1,6,11];
    
    figure('Position',[1920 40 1920 960],'Renderer','painters')
    %Show spatial outline of the component
    for rr=1:3
        subplot(4,5,ROI_outline_order(rr))
        imagesc(templates{1}.template);
        hold on
        title(num2str(ROI_all(rr)))
        axes(gca);
        axis square
        xticks(gca,[])
        yticks(gca,[])
        
        grayMap = brighten(gray,0.5);
        colormap(gca,grayMap)
        %plot componenet outline
        plot(coor_keep_sel{ROI_all(rr)}(1,:),coor_keep_sel{ROI_all(rr)}(2,:),ROI_colors{rr}, 'LineWidth',1);
        %zoom into component based on center of mass calculation
        xlim([A_keep_sel_com(ROI_all(rr),2)-30, A_keep_sel_com(ROI_all(rr),2)+30])
        ylim([A_keep_sel_com(ROI_all(rr),1)-30, A_keep_sel_com(ROI_all(rr),1)+30])
    end

%show color-coded (by trial) calcium traces of the component
ROI_trace_order = [2,3,4; 7 8 9; 12 13 14];
for rr=1:3
    subplot(4,5,ROI_trace_order(rr,:))
    hold on
    xlim([0.1 18.2])
    ylim([-0.2 4])
    yticks([0 1 2])
    ylabel('dF/F');
    xticks(0:3:15.5);
    set(gca,'FontSize',14)
    set(gca,'LineWidth',1)
    %convert to minuntes and offset to start at 0 min (0s)
    %for each split (cell above) % for A trial traces
    for ii=1:size(start_end_idx.A,1)
        plot(Imaging_split{1}{4}.time_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2))/60,...
            Imaging_split{1}{4}.trace_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2),ROI_all(rr)),'Color',[0 0 1 1],...
            'LineWidth',1)
    end
    %for each split (cell above) % for B trial traces
    for ii=1:size(start_end_idx.B,1)
        plot(Imaging_split{1}{5}.time_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2))/60,...
            Imaging_split{1}{5}.trace_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2),ROI_all(rr)),'Color',[1 0 0 1],...
            'LineWidth',1)
    end
end

ROI_polar_order = [5,10,15];
for rr=1:3
    %tuning vectors of the component
    pax1 = subplot(4,5,ROI_polar_order(rr),polaraxes);
    hold on
    pax1.FontSize = 14;
    angles = 0;
    pax1.ThetaTick = angles;
    thetaticklabels(pax1,{'lap start'});
    rticks(pax1, []);
    rlim(pax1,[0 25]);
    
    polarplot(x{1},r_scaled{1},'k','Linewidth',1.5)
    
    %plot A (2) trial events
    for ll=1:size(idxMin{1}{1},2)
        polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{1}{ll}{ROI_all(rr)}),'bo','MarkerFaceColor','b')
        %place field center
        %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
    end
    %plot B (3) trial events
    for ll=1:size(idxMin{1}{2},2)
        polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{2}{ll}{ROI_all(rr)}),'ro','MarkerFaceColor','r')
        %place field center
        %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
    end
end

%spiral map of the componenet
%subplot(4,5,10)

%correct A
% compass(Place_cell{1}{1}.Tuning_Specificity.tuning_vector{ROI},'b');
% hold on
% %correct B
% compass(Place_cell{1}{2}.Tuning_Specificity.tuning_vector{ROI},'r');
% %deal with labels
% labels = findall(gca,'type','text');
% set(labels,'visible','off');
%set(findall(gcf, 'String', '30', '-or','String','60') ,'String', '  ');

%only correct
%only works if using polar axis; plot in cartesian ref frame when using
%compass
%polarplot([0+0i,Place_cell{1}{1}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)

%plot normalized position
subplot(4,5,[17 18 19])
hold on
yticks([0 0.5 1])
ylabel('Normalized position')
xlabel('Time [min]');
xticks(0:3:9.5);
xlim([0.1 10.2])
ylim([0 1])
set(gca,'FontSize',14)
set(gca,'LineWidth',1)
%A laps
for ii=1:size(lap_idxs.A,1)
    plot(Imaging_split{1}{4}.time_restricted(lap_idxs.A(ii,1):lap_idxs.A(ii,2))/60,...
        Behavior_split{1}{4}.resampled.position_norm(lap_idxs.A(ii,1):lap_idxs.A(ii,2)),...
        'Color',[0 0 1 0.6],'LineWidth',1.5)
end
%B laps
for ii=1:size(lap_idxs.B,1)
    plot(Imaging_split{1}{5}.time_restricted(lap_idxs.B(ii,1):lap_idxs.B(ii,2))/60,...
        Behavior_split{1}{5}.resampled.position_norm(lap_idxs.B(ii,1):lap_idxs.B(ii,2)),...
        'Color',[1 0 0 0.6],'LineWidth',1.5)
end
%overlay significant calcium run events
%A
%scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
%B 
%scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
%pause;
%clf;

    
end

%% NOT USED IMPORTED CODE
% figure('Position',[2600,300,1200,1000]);
% for ii=options.idx_show%1:size(registered.multi.assigned_all,1)
%     
%     %ROI from session 1
%     ROI = registered.multi.assigned_all(ii,1);
%     
%     subplot(2,3,1)
%     imagesc(session_vars{1}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
%     hold on;
%     title(num2str(ROI));
%     ylabel('Lap #'); 
%     xlabel('Spatial bin');
%     caxis([0 2])
%     colormap(gca,'jet');
%     hold off;
%     
%     
%     
%     %spiral plot early in learning
%     subplot(2,3,2)
%     polarplot(x{1},r_scaled{1},'k','Linewidth',1.5)
%     hold on
%     
%     %plot A (2) trial events
%     for ll=1:size(idxMin{1}{1},2)
%         polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI}),r_scaled{1}(idxMin{1}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
%         %place field center
%         %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
%     end
%     
%     %plot tuning specificity vector for all A trials
%     polarplot([0+0i,15*session_vars{1, 1}.Place_cell{1, 4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
% 
% 
%     %plot B (3) trial events
%     for ll=1:size(idxMin{1}{2},2)
%         polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI}),r_scaled{1}(idxMin{1}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
%         %place field center
%         %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
%     end
%     
%     %plot tuning specificity vector for all B trials
%     polarplot([0+0i,15*session_vars{1, 1}.Place_cell{1, 5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
% 
% 
%     hold off
%     
%     subplot(2,3,3)
%     imagesc(ROI_zooms{ii,1})
%     hold on;
%     colormap(gca, 'gray')
%     xticks([])
%     yticks([])
%     b = bwboundaries(ROI_outlines{ii,1},'noholes');
%     plot(b{1}(:,2),b{1}(:,1),'r')
%     hold off
%     
%     %ROI from session 2
%     ROI = registered.multi.assigned_all(ii,2);
%     subplot(2,3,4)
%     imagesc(session_vars{2}.Place_cell{1, 3}.dF_lap_map_ROI{ROI})
%     hold on;
%     title(num2str(ROI));
%     ylabel('Lap #'); 
%     xlabel('Spatial bin');
%     caxis([0 2])
%     colormap(gca, 'jet');
%     hold off;
%     
%    
%     
%     subplot(2,3,5)
%     polarplot(x{2},r_scaled{2},'k','Linewidth',1.5)
%     hold on
%     %plot A (2) trial events
%     for ll=1:size(idxMin{2}{1},2)
%         polarscatter(angle(posVectorApprox{2}{1}{ll}{ROI}),r_scaled{2}(idxMin{2}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
%         %place field center
%         %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
%     end
%     %plot tuning specificity vector for all A trials
%     polarplot([0+0i,15*session_vars{2}.Place_cell{4}.Tuning_Specificity.tuning_vector_specificity(ROI)],'b-','LineWidth',2)
% 
%     
%     %plot B (3) trial events
%     for ll=1:size(idxMin{2}{2},2)
%         polarscatter(angle(posVectorApprox{2}{2}{ll}{ROI}),r_scaled{2}(idxMin{2}{2}{ll}{ROI}),'ro','MarkerFaceColor','r')
%         %place field center
%         %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
%     end
%     
%     %plot tuning specificity vector for all A trials
%     polarplot([0+0i,15*session_vars{2}.Place_cell{5}.Tuning_Specificity.tuning_vector_specificity(ROI)],'r-','LineWidth',2)
% 
%     hold off
%     
%     subplot(2,3,6)
%     imagesc(ROI_zooms{ii,2})
%     %imagesc(ROI_zooms{registered.multi.assigned_all(ii,2),2})
%     hold on;
%     colormap(gca, 'gray')
%     xticks([])
%     yticks([])
%     %b = bwboundaries(ROI_outlines{registered.multi.assigned_all(ii,2),2},'noholes');
%     b = bwboundaries(ROI_outlines{ii,2},'noholes');
%     plot(b{1}(:,2),b{1}(:,1),'r')
%     hold off
%     
%     %pause;
%     %clf;
%     
% end


end

