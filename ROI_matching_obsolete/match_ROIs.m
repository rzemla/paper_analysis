
%% Read in directory list
%directories containing the output of CNMF from each exp

% path_dir = {'G:\I52RT_AB_sal_PSEM_113018_120118\1',...
%     'G:\I52RT_AB_sal_PSEM_113018_120118\2',...
%     'G:\I52RT_AB_sal_PSEM_113018_120118\3',...
%     'G:\I52RT_AB_sal_PSEM_113018_120118\4'};
% 
% %cross session directory
% crossdir = 'G:\I52RT_AB_sal_PSEM_113018_120118\crossSession';

path_dir = {'G:\GOL\I55_RTLS_RF_GOL_022219\1',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\2',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\3',...
    'G:\GOL\I55_RTLS_RF_GOL_022219\4'};

%cross session directory
crossdir = 'G:\GOL\I55_RTLS_RF_GOL_022219\crossSession';

%directory path that can be used with fprintf for python config file
%modify this so it works with any path subfolders
%place in for loop and run it across all cells
split_crossdir = split(crossdir,'\');
%crossdir_config = [split_crossdir{1},'\\',split_crossdir{2},'\\',split_crossdir{3}];
crossdir_config = [split_crossdir{1},'\\',split_crossdir{2},'\\',split_crossdir{3},'\\',split_crossdir{4}];
%crossdir_config = 'E:\\I42L_AB_1\\crossDay';

%get the names of the mat files
for ii=1:size(path_dir,2)
    filenam(ii) = dir([path_dir{ii},'\input\*.mat']);
end

%navigate to save directory (cross session directory)
cd(crossdir)
%make python input and output directory
mkdir('python_input');
mkdir('python_output');

%% Read in template and A_keep matrix

%load in removed ROIs
for ii=1:size(path_dir,2)
    %get path name
    session_files_removedROI{ii} = subdir(fullfile(path_dir{ii},'removedROI','*selectedROIs.mat'));
    %load in the logical list of components
    somaticROI{ii} = load(session_files_removedROI{ii}.name,'compSelect');
end

%load all the remaining data for session-by-session comp registration
tic;
for ii=1:size(path_dir,2) %for each session
    session_vars{ii} = load(fullfile(filenam(ii).folder,filenam(ii).name),'A_keep','C_full','options_mc','options','dims','Coor_kp','expDffMedZeroed', 'template');
    
    %load template (for later datasets)
    template{ii} = load(fullfile(path_dir{ii},'template_nr.mat'));
    session_vars{ii}.template = template{ii}.template;
    
    %get keep centers
    session_vars{ii}.centers = com(session_vars{ii}.A_keep,session_vars{ii}.dims(1),session_vars{ii}.dims(2));
end
toc;

%% Set the registration options parameters

% options                 parameter structure with inputs:
%       d1:               number of rows in FOV
options.d1 = session_vars{1}.dims(1);
%       d2:               number of columns in FOV
options.d2 = session_vars{1}.dims(2);
%       options_mc:       motion correction options
options.options_mc = session_vars{1}.options_mc;
%       plot_reg:         create a contour plot of registered ROIs
options.plot_reg = true;

%TINKER WITH THESE PARAMETERS - works well out of the box

%       d3:               number of planes in FOV (default: 1)
%       dist_maxthr:      threshold for turning spatial components into binary masks (default: 0.1)
%       dist_exp:         power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n (default: 1)
%       dist_thr:         threshold for setting a distance to infinity. (default: 0.5)
%       dist_overlap_thr: overlap threshold for detecting if one ROI is a subset of another (default: 0.8)


%       options_mc:       motion correction options
options_mc = session_vars{1}.options_mc;

%correct filename for unix notation
mc_filename = split(options_mc.h5_filename,'/');
options_mc.h5_filename = mc_filename{1};

%% Remove non-somatic componenets

for ii=1:size(path_dir,2)
    compSelect = somaticROI{ii}.compSelect;
    session_vars{ii}.A_keep = session_vars{ii}.A_keep(:,compSelect);
    session_vars{ii}.C_full = session_vars{ii}.C_full(compSelect,:);
    session_vars{ii}.Coor_kp = session_vars{ii}.Coor_kp(compSelect);
    session_vars{ii}.expDffMedZeroed = session_vars{ii}.expDffMedZeroed(compSelect,:);
    session_vars{ii}.centers = session_vars{ii}.centers(compSelect,:);
end

%% Session-by-session registration - create cross match matrix

%define lookup matrix for which combinations to run
lookup_match = triu(ones(size(path_dir,2)),1);

%define equivalent empty cell to store output from component registration
register = cell(size(path_dir,2),size(path_dir,2));

tic;
for ii=1:size(path_dir,2)
    disp([num2str(ii), '/',num2str(size(path_dir,2))]);
    for jj=1:size(path_dir,2)
        
        %check that this combination is to be evaluated
        if lookup_match(ii,jj) == 1
            %define spatial componenets from session 1
            A1 = session_vars{ii}.A_keep;
            %define spatial componenets from session 1
            A2 = session_vars{jj}.A_keep;
            %define the templates
            template1 = session_vars{ii}.template;
            template2 = session_vars{jj}.template;
            
%             [register{ii,jj}.matched_ROIs,register{ii,jj}.nonmatched_1,register{ii,jj}.nonmatched_2,register{ii,jj}.Aout,register{ii,jj}.R,register{ii,jj}.A_union] = ...
%                 register_ROIs(A1,A2,options,template1,template2,options_mc);

            %save memory by excluding some outputs
            [register{ii,jj}.matched_ROIs,register{ii,jj}.nonmatched_1,register{ii,jj}.nonmatched_2,~,~,~] = ...
                register_ROIs(A1,A2,options,template1,template2,options_mc);
        end
    end
end
toc;


%% Export to python for multisession registration
% Read in A matrices, remove discarded components, and save to mat for processing in Python
% use different template read-in for current datasets
A= {};
for ii=1:size(path_dir,2)
    %load in A matrices
    load(fullfile(filenam(ii).folder,filenam(ii).name),'A_keep');
    %convert to full type
    A_keep = full(A_keep);
    
    %remove ROIs that were discarded
    A_keep = A_keep(:,somaticROI{ii}.compSelect);
    
    %save each session A for later visualization (after removed ROIs)
    A{ii} = A_keep;
    
    %load in templates
    load(fullfile(filenam(ii).folder,filenam(ii).name),'template');
    
    %save A_keep and templates for each session in numerical order
    save(fullfile(crossdir,'python_input',[num2str(ii),'.mat']),'A_keep','template')
end

%% Create config file for python to know cross session directory
fileID = fopen(['D:\python_scripts\','crossSesParams.config'], 'w');
%split each absolute location directory into subpart delimited by /
%splitName = strsplit(names{ii},'/');

fprintf(fileID,'[fileNames]\n\n');
fprintf(fileID,['exp_cross_name: ',crossdir_config, '\n']);
%fprintf(fileID,'save_directory: /scratch/rz646/data/stacks/\n');
%fprintf(fileID,'read_directory: /beegfs/rz646/data/stacks/\n\n');

% fprintf(fileID,'[buildParam]\n\n');
% 
% fprintf(fileID,'save_combined_stack = 0\n');
% fprintf(fileID,'split_stack = 0\n');
% fprintf(fileID,'split_every = 2000');

%close text file
fclose(fileID);


%% Execute python scripts from windows command line

%run bat script that load base/caiman environment and executes registration
%script
[stat,cmd_out] = system('D:\python_scripts\launch_python_script.bat','-echo');

%% Read python caiman register function output variables

%navigate to cross session working directory for experiment
cd(crossdir)

%load in python output variables
load(fullfile(crossdir,'python_output','matching_output.mat'))

%account for python 1 index offset (run once)
assignments_adj = assignments + 1;

%components that match across all sessions
all_match_sum = sum(assignments_adj,2);
assign_all_matching = assignments_adj(~isnan(all_match_sum),:);


%% Try plotting the matching A component

%from all matching indices across sessions
%number of sessions
nbSes = size(assign_all_matching,2); 

figure;
for selComp =1:size(assign_all_matching,1)
    for ii=1:nbSes
        subplot(1,nbSes,ii)
        imagesc( reshape(A{ii}(:,assign_all_matching(selComp,ii)),512,512))
        title(num2str(selComp))
    end
    pause;
end

%% Check the component binaries for how they matched

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

%% Plot as one large group of subplots of ROI across sessions

%number of ROIs (rows) by sessions (cols)
rows = 20;
cols = 4;

%matrix corresponding to display
display_matrix = reshape([1:(rows*cols)],cols,rows)';

%total number of matched ROIs to display
displayROInb = size(ROI_zooms,1);

%ROI input list (in chunks of 20 start to finish
start_end_ROI = [(1:20:20*floor(displayROInb/20))',(20:20:20*floor(displayROInb/20))'];
start_end_ROI = [start_end_ROI; [start_end_ROI(end)+1, displayROInb]];

%iterate through each range
for iter=1:size(start_end_ROI,1)
    figure('renderer','painters','Position', [2200 100 300 900])
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

%% Create a matching matrix for session by session registration equivalent to  
%matched ROIs - use lookup_match matrix
%save as that in register struct, but outside
matchedROIs = cell(size(path_dir,2),size(path_dir,2));

%move ROIs that are matched into each cell as matrix
for ii=1:size(path_dir,2)
    for jj=1:size(path_dir,2)
        
        %check that this combination is to be evaluated
        if lookup_match(ii,jj) == 1
            matchedROIs{ii,jj} = register{ii,jj}.matched_ROIs;
        end
    end
end

%for matching across all sessions (will miss some assignments)

%generate all matches
% [C, i2, i3] = intersect(matchedROIs{1,2}(:,1),matchedROIs{1,3}(:,1),'stable');
% 
% %combine into 1 matrix
% combineMatch(:,1) = matchedROIs{1,2}(i2,1);
% combineMatch(:,2) = matchedROIs{1,2}(i2,2);
% combineMatch(:,3) = matchedROIs{1,3}(i3,2);

%sort

%% Process output matrices/cells that assign ROIs to each other across days

%all matches across sessions
multiSessionAlign = assignments_adj;
%neurons that are only matched across all sessions
multiSessionAlign_all = assign_all_matching;

%cell from session by session matching
sesBySesMatchROIs = matchedROIs;

%cell with session-by-session match - should be redundant with
%multiSessionAlign

%preallocate
allSesMatchROIs = cell(size(path_dir,2),size(path_dir,2));
%move ROIs that are matched into each cell as matrix
for ii=1:size(path_dir,2)
    for jj=1:size(path_dir,2)
        
        %check that this combination is to be evaluated
        if lookup_match(ii,jj) == 1
            selectTemp = ~isnan(sum(multiSessionAlign(:,[ii jj]),2));
            
            allSesMatchROIs{ii,jj} = multiSessionAlign(selectTemp,[ii jj]);
        end
    end
end

%% Output data (struct)

registered.multi.assigned = multiSessionAlign;
registered.multi.assigned_all = multiSessionAlign_all;
registered.multi.assign_cell = allSesMatchROIs;
registered.session.assign_cell = sesBySesMatchROIs;



