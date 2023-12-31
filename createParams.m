%% Create Parameters
%   Editor:    Daniel Elbich
%   Created:   3/14/19
%
%   Creates parameters .mat file for use with MVPA scripts
%
% 
%   Updated:   5/22/20

%   Major Updates:
%   -Removed ROI processing. Only directories are handled now
%   -fMRIPrep output data is not integrated into the pipeline
%   -Reorganized script to handle variable preprocessing.
%
%   Updated:   9/19/19
%
%   Major updates:
%   -Added flag for regressing reaction time (RT) from betas prior to
%   classification
%   -Added flag and variables for bootstrap MVPA classification
%   -Added ERS functionality
%+
%
%   Minor updates:
%   -Cleaned up code to get subject list/directories
%   -Cleaned up code saving mat files with ROI names to preven overwriting
%   between ROI & Searchlight analyses
%
%
%   Updated:   4/22/19
%
%   Major updates:
%   -Added section to create separate params file for use with Specify
%   Model script. This will be deleted during creation of final param mat
%   file.
%   -Rescling now completed through AFNI. Must loaded/added to system PATH
%
%   Minor updates:
%   -Searches common folder for ROI masks instead of project folder
%   -List dialog to select 1 or more regions for reslicing
%   -User input to select searchlight/no searchlight anlaysis
%   -User input to select conditions of interest
%   -Saves subjects, rois, and conds variables to var directory in project
%   folder


%% Set Pipeline Parameters

%  Set Path Variables

% Parent/Project Path - directory containing project data and analyses
directory.Project = inputdlg(['Enter path to parent directory: '],...
    'File Path',[1 35],{pwd});
directory.Project = directory.Project{:};

% Directories for raw functional data and behavioral data
rawData.funcDir = inputdlg(['Enter Path to raw functional data starting from: ' directory.Project],...
    'File Path',[1 35],{'e.g. CPM/subj/func'});
rawData.funcDir = [directory.Project filesep rawData.funcDir{:}];

% Directory of the behavioral data file
rawData.behavDir = inputdlg(['Enter Path to behavioral file starting from: ' directory.Project],...
    'File Path',[1 35],{'e.g. CPM/subj/behav'});
rawData.behavDir = rawData.behavDir{:};

% Extension of the behavioral file to help search
rawData.behavFile = inputdlg('Enter Trial Tag File Extension:',...
    'File extension',[1 35],{'e.g. xlsx, csv'});
rawData.behavFile = rawData.behavFile{:};

% Functional Data Preprocessing Program
preprocPipeline = questdlg('Select Preprocessing Pipeline',...
    'Confirm Pipeline',...
    'spm12','fmriprep','Cancel','fmriprep');

%% Task and Data Information

try
    % Task Name
    taskInfo.Name = inputdlg('Enter Task Name: [as written in functional filename]',...
        'Task Name',[1 35],{'e.g. encoding, nback, rest'});
    taskInfo.Name = taskInfo.Name{:};
    
    % Conditions of Interest
    taskInfo.Conditions = inputdlg('Enter Conditions of Interest: [separated by commas]',...
        'Conditions',[1 35],{'e.g. face, object, place'});
    
    % Conditions of Interest
    taskInfo.accuracyFlag = questdlg('Is accuracy possible? (i.e. interactive vs. static task)',...
        'Confirm Accuracy Potential',...
        'Yes','No','Cancel','No');
    
    % Length of task (in volumes)
    taskInfo.Datapoints = inputdlg('How many volumes [Length of task]?',...
        'Volume Number',[1 35],{'150'});
    taskInfo.Datapoints = str2double(taskInfo.Datapoints{:});
    
    taskInfo.Runs = inputdlg('How many runs [# of functional runs]?',...
        'Number of Runs',[1 35],{'1'});
    taskInfo.Runs = str2double(taskInfo.Runs{:});
    
    taskInfo.Trials = inputdlg('How many trials per run [sum trials for all conditions]?',...
        'Trial Number',[1 35],{'40'});
    taskInfo.Trials = str2double(taskInfo.Trials{:});
    
    
    %% Set Project Specific Variables
    
    % Please specify:
    % -The units used (i.e., 'scans' or 'secs')
    % -The TR or repition time
    Model.units = 'secs';
    %Model.TR    = 2.5;
    Model.TR  = inputdlg('Enter TR: [in seconds]',...
        'Repetition Time',[1 35],{'2'});
    Model.TR  = str2double(Model.TR{:});
    
    % Please specify if a mask is used during model estimation
    Mask.on   = 0; %Default - no mask used
    Mask.dir  = '/path/to/mask/directory';
    Mask.name = 'name_of_mask.nii';
    
    switch preprocPipeline
        case 'spm12'
            Func.dir         = [directory.Project filesep rawData.funcDir];
            Func.wildcard    = '^ar.*\.nii'; % File
            taskInfo.Slices = inputdlg('How many slices [# of slices in volume]?',...
                'Number of Slices',[1 35],{'50'});
            taskInfo.Slices = str2double(taskInfo.Slices{:});
            
        case 'fmriprep'
            Func.dir         = [directory.Project filesep rawData.funcDir];
            Func.wildcard    = ['^*' taskInfo.Name '_run-']; % File
            
    end
    
catch
    warning('Unable to set project information!.')
end

%% Project Paths

try
    % General path setup
    directory.Analysis = [directory.Project filesep 'multivariate'];
    directory.Model = [directory.Analysis filesep 'models' filesep ...
        'SingleTrialModel' taskInfo.Name];
    setenv('analysisPath',directory.Model);
    !mkdir -p $analysisPath
    
    analysisList = {'MVPA','RSA','ERS'};
    
    [index,tf] = listdlg('PromptString','Select Multivariate Analysis:',...
        'SelectionMode','Single','ListString',analysisList);
    
    classType = analysisList{index};
    
    % Select analysis type. No searchlight is default
    analysisType = questdlg('Select Analysis Level for Multivariate Test:',...
        'Confirm Searchlight',...
        'Searchlight','ROI','Cancel','ROI');

    % Account for RT/regress out. No is default
    regressRT.flag = questdlg('Regress out Reaction Time (RT)?',...
        'Confirm Regression',...
        'Yes','No','Cancel','No');
catch
    error('Unable to set project path information!');
end

%% Classification Flags

try
    % Bootstrapping Setup
    
    bootstrap.flag  = questdlg('Perform Bootstrap?',...
        'Confirm Bootstrap',...
        'Yes','No','Cancel','No');
    
    if strcmpi(bootstrap.flag,'Yes')==1
        bootstrap.numRuns     = taskInfo.Runs;
        bootstrap.numTrials   = taskInfo.Trials;
        bootstrap.perm        = inputdlg('How many permutations?',...
            'Permutation Number',[1 35],{'1000'});
        bootstrap.perm = str2double(bootstrap.perm{:});
        
        bootstrap.trialsPerRun = randperm(length...
            (1:bootstrap.numTrials));
        
        for i=1:bootstrap.numRuns
            bootstrap.structNames{i} = ['perm' num2str(i)];
        end
    end
    switch classType
        % MVPA Flags
        case 'MVPA'
            % Compute MVPA with entire run to train/test or individual
            % trials to train/test
            try
                trialAnalysis = 'Run';
                
                leaveRuns = NaN;
                
                if strcmpi(analysisType, 'Searchlight')==1
                    searchlight.Metric = questdlg('Searchlight Metric',...
                        'SearchlightSize','count','radius','Cancel','count');
                    
                    searchlight.Size = inputdlg('Size of the searchlight',...
                        'Searchlight Size',[1 35],{'50'});
                    searchlight.Size = str2double(searchlight.Size{:});
                end
                
            catch
                warning(['Error in setting MVPA analysis flags. '...
                    'Set to debug mode.']);
            end
            
            % RSA Flags
        case 'RSA'
            % Compute RSA with mean activation pattern (average all trials)
            % or individual trial pattern (all trials are separate)
            try
                trialAnalysis = 'Individual';
                
                leaveRuns = NaN;
             
                if strcmpi(analysisType, 'Searchlight')==1
                    searchlight.Metric = questdlg('Searchlight Metric',...
                        'SearchlightSize','count','radius','Cancel','count');
                    
                    searchlight.Size = inputdlg('Size of the searchlight',...
                        'Searchlight Size',[1 35],{'50'});
                    searchlight.Size = str2double(searchlight.Size{:});
                end
            
            catch
                warning(['Error in setting RSA analysis flags. '...
                    'Set to debug mode.']);
            end
            
        case 'ERS'
            try
                [index,tf] = listdlg('Name','Possible Tasks',...
                    'PromptString','Ctrl+Click to select conditions:',...
                    'ListString',tasks,...
                    'ListSize',[280,300]);
                tasks=tasks(index);
                
                save([Project 'vars/tasks.mat'],'tasks');
                
            catch
                warning(['Error in setting ERS analysis flags. '...
                    'Set to debug mode.']);
            end
    end
catch
    warning('Error in setting classification flags. Set to debug mode.');
end

%% Create Subject List

try
    % List subject directory
    subjDir=dir(rawData.funcDir);
    
    % Remove any files in directory
    for i=1:length(subjDir)
        if subjDir(i).isdir==0
            fileFlag(i) = 0;
        else
            fileFlag(i) = 1;
        end
    end
    subjDir = subjDir(logical(fileFlag));
    
    % Remove non-subject directories
    % change subject tags according to project** 
    for i=1:length(subjDir)
        subFlag(i) = ~isempty(strfind(subjDir(i).name, 'sub-')); %change here to reflect you subject naming
    end
    subjDir = subjDir(subFlag);
    
    % Search directory to get list of subjects
    for i=1:length(subjDir)
        subjects{i,1}=subjDir(i).name;
    end
    
    clear subjCount i subjDir;
catch
    warning('Unable to create subject list. Set to debug mode.');
end

%% Assign Conditions

try
    taskInfo.Conditions = strsplit(taskInfo.Conditions{:},',');
    
    subconditionFlag = questdlg('Are there subconditions? (If unsure select No)',...
        'Confirm Subconditions','Yes','No','Cancel','No');
        
    switch subconditionFlag
        case 'Yes'
            [index,tf] = listdlg('Name','Possible Conditions',...
                'PromptString','Ctrl+Click to select conditions:',...
                'ListString',conds,...
                'ListSize',[280,300]);
            
            subconds=subconds(index);
            subconditionFlag = 'TRUE';
            
        case 'No'
            subconditionFlag = 'FALSE';
    end
    
catch
    warning('Unable to assign conditions. Set to debug mode.');
end

%% Set filename for output

try
    % Set filename for parameter .mat file
    for i=1:length(taskInfo.Conditions)
        if i==1
            conds=taskInfo.Conditions{i};
        else
            conds=[conds '-' taskInfo.Conditions{i}];
        end
        
    end
    switch analysisType
        case 'Searchlight'
            filename=strcat(directory.Analysis,'/params_',preprocPipeline,...
                '_',taskInfo.Name,'_',classType,'_',analysisType,'_Conds_',...
                conds,'_subConds_',subconditionFlag,'_Bootstrap_',...
                bootstrap.flag,'_Searchlight',searchlight.Metric,'_',...
                num2str(searchlight.Size),'.mat');
        otherwise
            filename=strcat(directory.Analysis,'/params_',preprocPipeline,...
                '_',taskInfo.Name,'_',classType,'_',analysisType,'_Conds_',...
                conds,'_subConds_',subconditionFlag,'_Bootstrap_',...
                bootstrap.flag,'_Searchlight_No.mat');
    end
    
catch
    warning('Unable to set params filename. Set to debug mode.');
end

%% Save Params File
switch analysisType
    case 'Searchlight'
        save(filename,'directory','rawData','preprocPipeline',...
            'taskInfo','Model','Mask', 'Func','classType',...
            'subjects','analysisType','regressRT','searchlight',...
            'bootstrap','trialAnalysis','leaveRuns');
    otherwise
        save(filename,'directory','rawData','preprocPipeline',...
            'taskInfo','Model','Mask','Func','classType',...
            'subjects','analysisType','regressRT',...
            'bootstrap'); %'trialAnalysis' ,'leaveRuns'
end

%% Create PBS Script for Command Line Processing

% Identify shell script
cmd=which('multi-FRAME/modeling/createPBSScript.sh');
setenv('cmd',cmd);

% Break filename into parts
[filepath,name]=fileparts(filename);

% Set Variables
setenv('fileName',name);
setenv('filePath',filepath);
setenv('project',directory.Project);

% Create PBS script
!sh $cmd --paramsFile $fileName --paramsDir $filePath --projectDir $project

%% Cleanup

clear;
clc;
disp('All finished!!');
            
