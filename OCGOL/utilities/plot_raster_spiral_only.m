function [outputArg1,outputArg2] = plot_raster_spiral_only(plot_raster_vars,session_vars,templates,task_remapping_ROIs,options)

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

%% Get time limits and other parameters for plotting

%imaging frame-to-frame time interval (period)
dt = session_vars{1}.Imaging.dt;

%number of laps
nb_laps = size(session_vars{1}.Behavior.lap,2);
%
startT = session_vars{1}.Imaging.time_restricted(1)/60;
%end seconds
endT = size(session_vars{1}.Imaging.time_restricted,1)*dt;
%end minutes
endT = endT/60;

%axis time limits for visualization
time_limits = [startT endT];
%df/f trace axis limits
dff_limits = [-0.2 4.5];

%red and blue colormaps
cmap_blue=cbrewer('seq', 'Blues', 64);
cmap_red=cbrewer('seq', 'Reds', 64);

%set bottom value to red to bottom value of blue
cmap_red(1,:) = [1 1 1];
cmap_blue(1,:) = [1 1 1];

%trial order of laps
trialOrder =session_vars{1}.Behavior.performance.trialOrder;
%all B lap idxs
idx_A = find(trialOrder ==2 | trialOrder ==20);
%all B lap idxs
idx_B = find(trialOrder ==3 | trialOrder ==30);

%% Plot each dF/F raster in a separate figure window

%define each class of ROI to plot
ROI_all = [627 416	269 524 458];
for rr=1:5
    figure('Position',[2951 443 297 313])
    %split into 2 and fill with nans areas where opposite trial is there
    whole_map_temp = session_vars{1}.Place_cell{3}.dF_lap_map_ROI_smooth{ROI_all(rr)};
    whole_map_temp_A = whole_map_temp;
    whole_map_temp_B = whole_map_temp;
    %fill with nans
    whole_map_temp_A(idx_B,:) = nan;
    whole_map_temp_B(idx_A,:) = nan;

    %plot A map
    imAlpha=ones(size(whole_map_temp));
    imAlpha(isnan(whole_map_temp_A))=0;
    ax2 = axes();
    h2 = imagesc(ax2,whole_map_temp_A,'AlphaData',imAlpha);
    %set colormap to
    colormap(ax2,cmap_blue);
    caxis(ax2,[0 2])
        %plot B map
    imAlpha=ones(size(whole_map_temp));
    imAlpha(isnan(whole_map_temp_B))=0;
    ax3 = axes;
    h3 = imagesc(ax3,whole_map_temp_B,'AlphaData',imAlpha);
    colormap(ax3,cmap_red);
    caxis(ax3,[0 2])
    ax2.Visible = 'off';
    ax3.Visible = 'off';

%export each raster
mkdir(fullfile(crossdir,'match_STC'))
disp('Saving match ROIs STC ')
export_fig(f_clip ,fullfile(crossdir,'match_STC','all_matching__nan_d1_clipped_300.png'),'-r300')

end

%make subplots with colorbar figure with colorbar
figure
subplot(1,2,1)
hold on
    ax1 = gca;
    colormap(ax1,cmap_blue);
    caxis(ax1,[0 2])
    c1 = colorbar(ax1);
c1.Ticks = [0 0.5 1 1.5 2];
c1.Label.String = 'dF/F';
c1.Label.FontSize = 16;
c1.FontSize = 16;
axis off

subplot(1,2,2)
hold on
    ax2=gca;
    c2 = colormap(ax2,cmap_red);
    caxis(ax2,[0 2])
c2 = colorbar(ax2)
c2.Ticks = [0 0.5 1 1.5 2];
c2.Label.String = 'dF/F';
c2.Label.FontSize = 16;
c2.FontSize = 16;
axis off

%% Plot figure 2b traces 
if options.plotFigure2 ==1
    

%go through each class of ROIs
%for cc=1:size(task_remapping_ROIs.global_near,2)
    %common, rate, near
    %ROI_all = [627	416	440 524 522];
    %ROI_all = [627	416	269 400 522];
    %ROI_all = [627	416	381 524 169];
    %ROI_all = [627	416	381 524 470];
    ROI_all = [627 416	269 524 458]; % global:78  524 228
% partial select: 640, 591, 534 (nice), 522 (nice), 519 (nice) , 511, 470 (nice), 458 (nice), 437 (nice), 371, 343,
% 244 (nice), 233(nice), 206 (nice), 187 (nice), 176 (nice), 169 (nice) , 

for zz=1
    %ROI_all = [562	404	455];
    %ROI_all(1:3) = task_remapping_ROIs.global_near(cc)
    %ROI_colors = {'b','r','m'};
    %matrix of RGB color values
    %ROI_colors = [65,105,225, 255; 220,20,60, 255; 139, 0, 139, 255]/255;
    %all purple
    ROI_colors = [139, 0, 139, 255; 139, 0, 139, 255; 139, 0, 139, 255; 139, 0, 139, 255; 139, 0, 139, 255]/255;
    %ROI_outline_order = [1,6,11, 16, 21];
    ROI_outline_order = 1:6:36;
    marker_size = 8;
    spiral_line_width = 0.4;
    
    figure('Position',[1920 40 1920 960],'Renderer','painters')
    %Show spatial outline of the component
    for rr=1:5
        subplot(6,6,ROI_outline_order(rr))
        imagesc(templates{1}.template);
        hold on
        %set background axis color to black
        set(gca,'color',0*[1 1 1]);
        title(num2str(ROI_all(rr)))
        axes(gca);
        axis square
        xticks(gca,[])
        yticks(gca,[])
        ax = gca;
        grayMap = brighten(gray,0.4);
        colormap(ax,grayMap)
        %plot componenet outline
        plot(coor_keep_sel{ROI_all(rr)}(1,:),coor_keep_sel{ROI_all(rr)}(2,:),'Color',ROI_colors(rr,:), 'LineWidth',1);
        %zoom into component based on center of mass calculation
        xlim([A_keep_sel_com(ROI_all(rr),2)-15, A_keep_sel_com(ROI_all(rr),2)+15])
        ylim([A_keep_sel_com(ROI_all(rr),1)-15, A_keep_sel_com(ROI_all(rr),1)+15])
hold off
    end
    
    %show color-coded (by trial) calcium traces of the component
    %ROI_trace_order = [2,3,4; 7 8 9; 12 13 14; 17 18 19; 22 23 24];
    ROI_trace_order(:,1) = 2:6:30;
    ROI_trace_order(:,2) = 3:6:30;
    ROI_trace_order(:,3) = 4:6:30;

    for rr=1:5
        subplot(6,6,ROI_trace_order(rr,:))
        hold on
        xlim(time_limits)
        ylim(dff_limits)
        yticks([0 1 2])
        ylabel('dF/F');
        xticks(0:3:15.5);
        set(gca,'FontSize',14)
        set(gca,'LineWidth',1)
        %convert to minuntes and offset to start at 0 min (0s)
        %for each split (cell above) % for A trial traces
        for ii=1:size(start_end_idx.A,1)
            plot(Imaging_split{1}{4}.time_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2))/60,...
                Imaging_split{1}{4}.trace_restricted(start_end_idx.A(ii,1):start_end_idx.A(ii,2),ROI_all(rr)),'Color',[65,105,225, 255]/255,...
                'LineWidth',1)
        end
        %for each split (cell above) % for B trial traces
        for ii=1:size(start_end_idx.B,1)
            plot(Imaging_split{1}{5}.time_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2))/60,...
                Imaging_split{1}{5}.trace_restricted(start_end_idx.B(ii,1):start_end_idx.B(ii,2),ROI_all(rr)),'Color',[220,20,60, 255]/255,...
                'LineWidth',1)
        end
    end
    
    %ROI_polar_order = [5,10,15,20,25];
    ROI_polar_order = 5:6:30;
    for rr=1:5
        %tuning vectors of the component
        pax1 = subplot(6,6,ROI_polar_order(rr),polaraxes);
        hold on
        pax1.FontSize = 14;
        angles = 0;
        pax1.ThetaTick = angles;
        thetaticklabels(pax1,{'lap start'});
        rticks(pax1, []);
        rlim(pax1,[0 nb_laps]);
        
        polarplot(x{1},r_scaled{1},'k','Linewidth',spiral_line_width)
        
        %plot A (2) trial events
        for ll=1:size(idxMin{1}{1},2)
            polarscatter(angle(posVectorApprox{1}{1}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{1}{ll}{ROI_all(rr)}),...
                marker_size,'MarkerEdgeColor',[65,105,225]/255,'Marker','o','MarkerFaceColor',[65,105,225]/255)
            %place field center
            %polarscatter(centerA_angle(ii), 20, 'b*','MarkerFaceColor','b');
        end
        %plot B (3) trial events
        for ll=1:size(idxMin{1}{2},2)
            polarscatter(angle(posVectorApprox{1}{2}{ll}{ROI_all(rr)}),r_scaled{1}(idxMin{1}{2}{ll}{ROI_all(rr)}),...
                marker_size,'MarkerEdgeColor',[220,20,60]/255, 'Marker','o','MarkerFaceColor',[220,20,60]/255)
            %place field center
            %polarscatter(centerB_angle(ii), 20, 'r*','MarkerFaceColor','r');
        end
    end
%end

%dF/F mean rasters for each ROI
ROI_dFF_order = 6:6:36;
for rr=1:5
    subplot(6,6,ROI_dFF_order(rr))
    axis off
    %split into 2 and fill with nans areas where opposite trial is there
    whole_map_temp = session_vars{1}.Place_cell{3}.dF_lap_map_ROI_smooth{ROI_all(rr)};
    whole_map_temp_A = whole_map_temp;
    whole_map_temp_B = whole_map_temp;
    %fill with nans
    whole_map_temp_A(idx_B,:) = nan;
    whole_map_temp_B(idx_A,:) = nan;

    %plot A map
    imAlpha=ones(size(whole_map_temp));
    imAlpha(isnan(whole_map_temp_A))=0;
    ax2 = axes;
    h2 = imagesc(ax2,whole_map_temp_A,'AlphaData',imAlpha);
    %set colormap to
    colormap(ax2,cmap_blue);
    caxis(ax2,[0 2])
        %plot B map
    imAlpha=ones(size(whole_map_temp));
    imAlpha(isnan(whole_map_temp_B))=0;
    ax3 = axes;
    h3 = imagesc(ax3,whole_map_temp_B,'AlphaData',imAlpha);
    colormap(ax3,cmap_red);
    caxis(ax3,[0 2])
    ax2.Visible = 'off';
    ax3.Visible = 'off';


end

%plot normalized position
subplot(6,6,[32 33 34])
hold on
yticks([0 0.5 1])
ylabel('Normalized position')
xlabel('Time [min]');
xticks(0:3:11.2);
xlim(time_limits)
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

end


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


%overlay significant calcium run events
%A
%scatter(event_norm_time.A{ROI},event_norm_pos_run.A{ROI},[],[0 0 1],'*')
%B 
%scatter(event_norm_time.B{ROI},event_norm_pos_run.B{ROI},[],[1 0 0],'*')
%pause;
%clf;

    
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


end

