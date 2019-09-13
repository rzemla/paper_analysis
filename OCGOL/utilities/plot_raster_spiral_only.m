function [outputArg1,outputArg2] = plot_raster_spiral_only(plot_raster_vars,session_vars,templates,options)

%% Import variables

start_end_idx = plot_raster_vars.start_end_idx;
lap_idxs = plot_raster_vars.lap_idxs;

x = plot_raster_vars.x;
r_scaled = plot_raster_vars.r_scaled;

idxMin = plot_raster_vars.idxMin;
posVectorApprox = plot_raster_vars.posVectorApprox;

event_norm_time = plot_raster_vars.event_norm_time;
event_norm_pos_run = plot_raster_vars.event_norm_pos_run;

A_keep_sel_com = plot_raster_vars.A_keep_sel_com;
coor_keep_sel = plot_raster_vars.coor_keep_sel;

Imaging_split{1} = session_vars{1}.Imaging_split;
Place_cell{1} = session_vars{1}.Place_cell;
Behavior_split{1} = session_vars{1}.Behavior_split; 

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
            polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI}),r_scaled{1}(idxMin{1}{1}{ll}{ROI}),'bo','MarkerFaceColor','b')
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
    ROI_all = [398	404	455];
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


end

