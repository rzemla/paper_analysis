
function [C_df,CSV,XML,dff_input_vars, F_vars]=combine_sessions_newDff_RZ_V1(directory_name,type)

%navigate to directory of containing .mat, .xml, .csv files - RZ
cd(directory_name);

%list the mat., .xml, and .csv files in the folder - RZ
listmat = dir('*.mat'); %will look for all the mat files
listcsv = dir('*.csv');
listxml = dir('*.xml');

%sample ROI - change later to be able to select
ROI = 2;

%for each .mat file (containing the C_df matrix) - count based on C_df data
for i = 1:length(listmat)
    %take in the name of the CSV file (not whole directory)
    csvfile{i}=listcsv(i).name;
    
    disp(['Reading CSV ' ,num2str(i),'/', num2str(length(listmat))]);
    %read the CSV file as a matrix in a separate cell for each session
    %uses relative directory, so must navigate to it
    CSV{i} = csvread(csvfile{i},1,0);
    
    disp(['Reading XML ' ,num2str(i),'/', num2str(length(listmat))]);
    
    %must provide absolute directory name by using fullfile
    %each XMl file in read in into a struct
    XML{i} = xml2structV2((fullfile(listxml(i).folder,listxml(i).name)));
    
    %select which type of calcium signal to read in
    switch type
        %exponentially smoothed and median subtracted
        case 'expdff'
            %load in all variables - can make this faster by loading only
            %select variable of interest -RZ implement in FUTURE
            session_imaging{i} = load((listmat(i).name));
            C_df{i}=session_imaging{i}.expDffMedZeroed;
            %each column is an ROI
            C_df{i}=(full(C_df{i})');
            
        %C_df without smoothing and baseline shifting
        case 'Cdf'
            session_imaging{i} = load((listmat(i).name));
            C_df{i}=session_imaging{i}.C_df;
            C_df{i}=(full(C_df{i})');
            
        %deconvolved calcium signal
        case 'spikes'
            session_imaging{i} = load((listmat(i).name));
            C_df{i}=session_imaging{i}.S;
            C_df{i}=(full(C_df{i})');
        
        %calcium signal - intial dff - refined later in event detection
        %section
        case 'Jia/Danielson'
            
            ROI = 54;
            %load in necessary variables to calculate dFFs
            dff_input_vars = load((listmat(i).name), 'F');
            
            [F_vars] = dFF_Fiji_V1(dff_input_vars,ROI);
            
            %load the dFF types into cells
            %dFF_types{1} = F_vars.F_dff_exp; %Rolling median (with Jia/Danielson constant)
            %dFF_types{2} = F_vars.F_df_exp;  %Jia/Danielson
            C_df{i} = F_vars.F_df_exp';
            
        case 'Rolling median'
                        %load in necessary variables to calculate dFFs
            dff_input_vars = load((listmat(i).name), 'A_keep','b','FOV','C_full',...
                'f_full','F_dark','T','R_full', 'options');
            
            [F_vars] = dFF(dff_input_vars,ROI);
            
            %load the dFF types into cells
            dFF_types{1} = F_vars.F_dff_exp; %Rolling median (with Jia/Danielson constant)
            dFF_types{2} = F_vars.F_df_exp;  %Jia/Danielson
            C_df{i} = dFF_types{1}';
    end
end
end

