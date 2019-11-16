function [outputArg1,outputArg2] = cumulative_task_sel_STC(path_dir)


%% Load the STC data

%read relevant data
for aa=1:size(path_dir,2)
    load_data_path{aa} = fullfile(path_dir{aa},'cumul_analysis','task_sel_STC.mat');
    cumulative_STCs{aa} = load(string(load_data_path{aa}));
end

%% Generate a cumulative raster

%combine maps into single cell (animal x rasters)
%preallocate cell
STC_cell = cell(size(path_dir,2), 2);

for aa=1:size(path_dir,2)
    STC_cell(aa,:) = cumulative_STCs{aa}.task_sel_STC.maps{1};
end

%Split into A/B STC matrices
STC_A_sel = cell2mat(STC_cell(:,1));
STC_B_sel = cell2mat(STC_cell(:,2));

%% Global sort of A/B selective STCs

[~,maxBin_all.A] = max(STC_A_sel(:,1:100)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all.A] = sort(maxBin_all.A,'ascend');

[~,maxBin_all.B] = max(STC_B_sel(:,101:200)', [], 1);
%sortIdx - arrangment of ROIs after sorting by max spatial bin acitivity
[~,sortOrder_all.B] = sort(maxBin_all.B,'ascend');

figure
imagesc(STC_A_sel(sortOrder_all.A,:))
hold on

figure
imagesc(STC_B_sel(sortOrder_all.B,:))
hold on


%% Plot color-coded STCs

%reward A and B black marker lines
rewA_bin_pos = 30;
rewB_bin_pos = 70;
odor_bin_pos = 10;

%define colors maps for rasters
cmap_blue=cbrewer('seq', 'Blues', 32);
cmap_red=cbrewer('seq', 'Reds', 32);

%%%%% A selective neurons - STC on A laps and STC on B laps %%%%%

%set 0 value to white for both blue and red
cmap_blue(1,:) = [1 1 1];
cmap_red(1,:) = [1 1 1];

%input raster for display (A and B neighboring)
%for A selective neurons
figure('Position',[2106 101 761 857])
%A
subplot(2,2,1)
%A laps
input_matrix = STC_A_sel(sortOrder_all.A,1:100);
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
title('Asel - A laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_blue);
cbar= colorbar;
cbar.Label.String = 'Normalized activity';
cbar.Ticks = [0 0.5 1];
ax1 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax1.XTick = [1 100];
ax1.XTickLabel = {'0','1'};

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(STC_A_sel,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(STC_A_sel,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(STC_A_sel,1)],'k--')

%make ticks invisible
set(ax1, 'TickLength', [0 0]);

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B laps as input
input_matrix = STC_A_sel(sortOrder_all.A,101:200);

%B 
subplot(2,2,2)
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
title('Asel - B laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_red);
cbar2 = colorbar;
cbar2.Label.String = 'Normalized activity';
cbar2.Ticks = [0 0.5 1];
ax2 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax2.XTick = [1 100];
ax2.XTickLabel = {'0','1'};
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(STC_A_sel,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(STC_A_sel,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(STC_A_sel,1)],'k--')

%make ticks invisible
set(ax2, 'TickLength', [0 0]);


%%%%% B selective neurons - STC on A laps and STC on B laps %%%%%

% %set 0 value to white for both blue and red
% cmap_blue(1,:) = [1 1 1];
% cmap_red(1,:) = [1 1 1];

%input raster for display (A and B neighboring)
%for A selective neurons
%figure('Position',[2182 356 934 559])
%A
subplot(2,2,3)
%A laps
input_matrix = STC_B_sel(sortOrder_all.B,1:100);

%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
title('Bsel - A laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_blue);
cbar= colorbar;
cbar.Label.String = 'Normalized activity';
cbar.Ticks = [0 0.5 1];
ax1 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax1.XTick = [1 100];
ax1.XTickLabel = {'0','1'};

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(STC_B_sel,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(STC_B_sel,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(STC_B_sel,1)],'k--')

%make ticks invisible
set(ax1, 'TickLength', [0 0]);

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%B laps as input
input_matrix = STC_B_sel(sortOrder_all.B,101:200);

%B 
subplot(2,2,4)
%create blank alpha shading matrix where 
imAlpha=ones(size(input_matrix));
imAlpha(isnan(input_matrix))=0;
imagesc(input_matrix,'AlphaData',imAlpha);
hold on
title('Bsel - B laps')
%set background axis color to black
set(gca,'color',0*[1 1 1]);
%set colormap to 
caxis([0 1])
colormap(gca,cmap_red);
cbar2 = colorbar;
cbar2.Label.String = 'Normalized activity';
cbar2.Ticks = [0 0.5 1];
ax2 = gca;
ylabel('Neuron #');
xlabel('Normalized position');
ax2.XTick = [1 100];
ax2.XTickLabel = {'0','1'};
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

%plot reward A and reward B zones
%A
plot([rewA_bin_pos rewA_bin_pos],[1 size(STC_B_sel,1)],'k--')
%B
plot([rewB_bin_pos rewB_bin_pos],[1 size(STC_B_sel,1)],'k--')
%odor pos
plot([odor_bin_pos odor_bin_pos],[1 size(STC_B_sel,1)],'k--')

%make ticks invisible
set(ax2, 'TickLength', [0 0]);



end

