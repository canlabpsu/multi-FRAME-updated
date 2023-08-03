%% Run Representational Similarity Analysis (RSA)
%   Editor:    Daniel Elbich
%   Updated:   4/29/19
%
% Representational similarity analysis (RSA) for a single subject or Encoding-Retrieval Similarity (ERS) for single subjects.
% Flagged for either single ROI RSA or searchlight analysis, but only
% currently supports ROI RSA. For RSA Searchlight, see runRSASearchlight.m
%
% Load single-trial beta images from each subject, apply ROI mask, calculate
% correlations between trial patterns, take the mean across trial types
%
%
% Updates:
%
% 3/28/19 - Loads in parameter file created by createParams.m subscript. Mat
% file should contain paths to roi and data folders, lists of rois,
% subjects, & conditions, and analysis name. See createParams.m for full
% list.

%Outputs within-category similarity, between-category similarity and
%distinctiveness described by Haxby 2001

% Only works for 2 conditions

% Requires zipped betas in nifti format and co-registered ROIs 

%% Set Analysis Parameters & Paths
%CG added 8/25A
addpath(genpath...
    ('/gpfs/group/nad12/default/nad12/toolbox/CoSMoMVPA-master/CoSMoMVPA-master'));

% Load all relevent project information
if exist('commandFlag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

% turn cosmo warnings off
cosmo_warning('off');

% Filepath for results folder
parentDir = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

%Debug
%subjects(2)=[];

%% Main Body
for iteration=1:length(subjects)
    %% Subject-Specific Directories
    
    % Current subject data paths:
    %  dataPath = fullpath to this subject's Single Trial Model directory
    %  spmFile  = fullpath to this subject
    %%%% For Just model folders w/Masks in's SPM.mat file. Note: the
    %                 :beta appended to the end tells cosmo to pull the beta
    %                 information from the SPM.mat file.
    
    dataPath   = fullfile(directory.Model, subjects{iteration});
    outputPath = fullfile(analysis, subjects{iteration});   
    spmFile = [dataPath '/SPM_gz.mat']; %For fMRIprep   
    
    % create the output path if it doesn't already exist
    if ~exist(outputPath, 'dir')
        mkdir(outputPath)
    end
    
    %% Load Mask Data
    %For SPM preprocessed
    %    masks = dir([directory.Analysis filesep 'masks' filesep file(1:end-4)...
    %    filesep subjects{iteration} filesep '*.nii.gz']);
    
    %For fMRIprep preprocessed
  
    %%%% For Individual Subject Masks folder setup
         %masks = dir([directory.Analysis filesep 'masks' filename(50:end-4)... 
         %filesep subjects{iteration} filesep '*.nii']);

    %%%% For Just model folders w/Masks inside
        %masks = dir([directory.Analysis filesep 'masks' ...
        %filesep subjects{iteration} filesep '*.nii']);  
   
   masks = dir([directory.Analysis filesep 'masks' filesep file(1:end-4)...
   filesep subjects{iteration} filesep '*.nii']); 
          
    for curMask = 1:length(masks)
        
        switch classType
            case 'RSA'
                % Path to current region mask
                curROI = fullfile(masks(curMask).folder, masks(curMask).name);
                
                % Current region name
                regionName=erase(masks(curMask).name,'.nii'); 
                %regionName=erase(masks(curMask).name,'.nii.gz');
                
                %%% Loading ROI data into CosmoMVPA - use SPM betas
                % Note that loading data through the SPM.mat file will automatically
                % "chunk" by runs, which is what we want
                fprintf('Current subject: %s\n',subjects{iteration});
                fprintf('Loading data from ROI: %s\n',masks(curMask).name);
                currDataset=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                
                % Clear errant Not a Number (NaN) values
                % Remove constant features
                currDataset=cosmo_remove_useless_data(currDataset);
                
                switch regressRT.flag
                    case 'Yes'
                        files = dir([study_path filesep subject filesep 'Run*']);
                        for i=1:length(files)
                            curMat(i) = load([files(i).folder filesep files(i).name]);
                            if i==length(files)
                                rtCell = [curMat.RT];
                                
                                % Convert from cell to double for regression
                                for ii=1:length(rtCell)
                                    
                                    % Flag outlier RT greater than 4 seconds
                                    if double(rtCell{ii}) >= regressRT.trialSec
                                        rtDouble(ii,1) = regressRT.trialSec;
                                    else
                                        rtDouble(ii,1) = double(rtCell{ii});
                                    end
                                end
                                
                                % Replace with trial duration (set in params)
                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                            end
                        end
                        
                        for i=1:length(currDataset.samples)
                            model = LinearModel.fit(rtDouble,currDataset.samples(:,i));
                            if i==1
                                allResiduals = model.Residuals.Raw;
                            else
                                allResiduals = [allResiduals model.Residuals.Raw];
                            end
                        end
                        
                        zscoreResid = zscore(allResiduals);
                        currDataset.samples = zscoreResid;
                        
                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                end
            case 'ERS'
                disp('Skipping to ERS...');
        end
        
        %% Define trial information
        switch classType
            case 'RSA'
                try
                    % Identify trials of interest for each condition
                    if exist('subConds','var')
                        
                        for ii=1:length(taskInfo.Conditions)
                            
                            subCond.(taskInfo.Conditions{1,ii})=contains(currDataset.sa.labels, taskInfo.Conditions{1,ii});
                            counter=1;
                            
                            for iii=1:length(subConds)
                                subCond.(taskInfo.Conditions{1,ii})(:,counter+1)=contains...
                                    (currDataset.sa.labels, subConds{1,iii});
                                counter=counter+1;
                            end
                            
                            subCond.(taskInfo.Conditions{1,ii})=double(subCond.(taskInfo.Conditions{1,ii}));
                            subCond.(taskInfo.Conditions{1,ii})(:,counter+1)=sum(subCond.(taskInfo.Conditions{1,ii})(:,1:counter),2);
                            
                        end
                        
                        Cond(1).idx = find(subCond.(taskInfo.Conditions{1,1})(:,counter+1) == counter);
                        Cond(2).idx = find(subCond.(taskInfo.Conditions{1,2})(:,counter+1) == counter);
                        
                    else
                        
                        CondList = zeros(size(currDataset.samples,1),1);
                        for ii=1:length(taskInfo.Conditions)
                            Cond(ii).labels = ~cellfun(@isempty, strfind...
                                (currDataset.sa.labels, taskInfo.Conditions{ii}));
                            Cond(ii).idx = find(Cond(ii).labels == 1);
                            
                            CondList(Cond(ii).idx) = ii;
                            
                        end
                        
                    end
                    
                    currDataset.sa.targets = CondList;
                    
                    % Codes trials/conditions of no interest as 0 (see SpecifyModel script
                    % for trial tag information)
                    Zeroidx = find(CondList == 0);
                    
                    % Removes all trials of no interest from analysis
                    if isempty(Zeroidx)==0
                        currDataset.samples(Zeroidx,:)=[];
                        currDataset.sa.targets(Zeroidx)=[];
                        currDataset.sa.beta_index(Zeroidx)=[];
                        currDataset.sa.chunks(Zeroidx)=[];
                        currDataset.sa.fname(Zeroidx)=[];
                        currDataset.sa.labels(Zeroidx)=[];
                    end
                    
                    fprintf('Number of possible targets: %i\n',length...
                        (unique(currDataset.sa.targets)));
                    
                    % Print dataset
                    fprintf('Dataset input:\n');
                    
                    % Print dataset
                    fprintf('Number of samples: %i\n',...
                        size(currDataset.samples,1));
                    fprintf('Number of features (voxels): %i\n',...
                        size(currDataset.samples,2));
                    fprintf('Number of chunks (runs): %i\n',...
                        length(unique(currDataset.sa.chunks)));
                    
                    % RSA ROI analysis
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    switch analysisType
                        case 'Searchlight' %%%CC added 3/30 not functional yet     
                       %%% Loading ROI data into CosmoMVPA - use SPM betas
                        fprintf('Loading data from ROI: %s\n',ROI);
          
                            spmFile = [parentDir '/SingleTrialModelenc' '/' ...
                                subject '/SPM_gz.mat']; %For fMRIprep
                            currDataset(i)=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                            
                        % Calculate separate matrices for each condition & compare
                            % those
                            
                            % Split dataset into separate condition
                            % variables
                            for i=1:length(Cond)
                                index = find(currDataset(i).sa.targets == i);
                                Conditions(i) = currDataset(i);
                                Conditions(i).samples = Conditions(i).samples(index,:);
                                Conditions(i).sa.beta_index = Conditions(i).sa.beta_index(index);
                                Conditions(i).sa.chunks = Conditions(i).sa.chunks(index);
                                Conditions(i).sa.fname = Conditions(i).sa.fname(index);
                                Conditions(i).sa.labels = Conditions(i).sa.labels(index);
                                
                                % Re-index trials within condition into
                                % separate targets
                                Conditions(i).sa.targets = [1:length(index)]';
                                Conditions(i).samples = Conditions(i).samples';   
                            end
                                
                   % Setup DSM and searchlight arguments
                        %DSM
                        dsmArgs.metric = 'correlation';
                        dsmArgs.center_data = 1;
                        
                        %Searchlight
                        searchlightSize = 20;
                        metric = 'count';
                        measure = @cosmo_target_dsm_corr_measure;
                        measure_args = struct();
                        measure_args.type = 'Spearman';
                        measure_args.center_data = 1;
                        
                        for i=1:length(Cond)
                            
                            % Create Target DSM from Encoding Run
                            targetDSM(i) = cosmo_dissimilarity_matrix_measure(Conditions(i), dsmArgs);
                            
                            if i==1
                                % Create Searchlight
                                nbrhood=cosmo_spherical_neighborhood(Conditions(i),metric,searchlightSize);
                            end
                            
                            % Set target DSM
                            measure_args.target_dsm = targetDSM(i).samples;
                            
                            % Searchlight RSA for each condition separately
                            results(i) = cosmo_searchlight(Conditions(i),nbrhood,measure,measure_args);
                            
                        end
                        
                        save([outputPath '/searchlightResults_' metric '_' ...
                            num2str(searchlightSize) '.mat'],'results');
                        
                        for i=1:length(conds)
                            % Define output location
                                output{i}=strcat(outputPath,'/',subject,...
                                    '_RSA_Searchlight_',metric,'_',...
                                    num2str(searchlightSize),'_',...
                                    conds{i},'.nii');
                            
                            % Store results to disk
                            SearchlightRSAoutput=cosmo_map2fmri(results(i), output{i});
                            
                        end
                            
                        case 'ROI'
                            
                            % Calculate separate matrices for each condition & compare
                            % those
                            
                            % Split dataset into separate condition
                            % variables
                            for i=1:length(Cond)
                                index = find(currDataset.sa.targets == i);
                                Conditions(i) = currDataset;
                                Conditions(i).samples = Conditions(i).samples(index,:);
                                Conditions(i).sa.beta_index = Conditions(i).sa.beta_index(index);
                                Conditions(i).sa.chunks = Conditions(i).sa.chunks(index);
                                Conditions(i).sa.fname = Conditions(i).sa.fname(index);
                                Conditions(i).sa.labels = Conditions(i).sa.labels(index);
                                
                                % Re-index trials within condition into
                                % separate targets
                                Conditions(i).sa.targets = [1:length(index)]';
                                Conditions(i).samples = Conditions(i).samples';
                               
                            end
                    end
                    
                    %Compute within-condition correlations the Haxby way
                    withinAvgCorr = {};
                    
                    for t=1:length(taskInfo.Conditions)
                        currCorrWithin=corrcoef(Conditions(t).samples);
                        withinAvgCorr{t} = mean(atanh(nonzeros(tril(currCorrWithin,-1))));
                    end
                    
                    averageRSAWithinForROI = mean(cell2mat(withinAvgCorr)),;
   
                    %Compute between-condition correlations the Haxby way
                    betweenAvgCorr = {};
                    
                    for j = 1:size(Conditions(1).samples,2)
                        for c = 1:size(Conditions(2).samples,2)
                            betweenAvgCorr{c,j} = corr(Conditions(1).samples(:,j),Conditions(2).samples(:,c),'Type','Pearson');
                        end
                    end
                    
                    averageRSABetweenForROI = mean2(atanh(cell2mat(betweenAvgCorr)));
                    
                    %Compute distinctiveness score
                    %distinctivenessRSAForROI = averageRSAWithinForROI - averageRSABetweenForROI;
             
                    % Obtain number of combinations with 2
                    % conditions selected
                    combinations=factorial(length(Cond))/...
                        ((factorial(2)*(factorial(length(Cond)-2))));
                    
                    
                    % Searchlight ERS for each condition separately
                    iter=1;
                    for i=1:combinations
                        if i<length(Cond)
                            condCount=1;
                        elseif i<combinations
                            condCount=2;
                            if i==length(Cond)
                                iter=condCount;
                            end
                        else
                            condCount=3;
                            iter=condCount;
                        end
                        
                        TrialTypeCombo(1,i) = {strcat(taskInfo.Conditions{condCount},'_v_',taskInfo.Conditions{iter+1})};
                        
                        iter=iter+1;
                    end
                    
                    clear targetDSM Conditions;
                end
                
                %% ERS analysis
            case 'ERS'
                switch analysisType
                    case 'Searchlight'
                        %%% Loading ROI data into CosmoMVPA - use SPM betas
                        fprintf('Loading data from ROI: %s\n',ROI);
                        for i=1:length(tasks)
                            
                            spmFile = [parentDir '/SingleTrialModel' tasks{i} '/' ...
                                %subject '/SPM_gz.mat']; %For fMRIprep
                                subject '/SPM.mat'];
                            %spmFile = [dataPath '/SPM_gz.mat']; %For fMRIprep
                            spmFile = [dataPath '/SPM.mat'];
                            currDataset(i)=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                            
                        end
                        
                        % Remove NaN voxels in both datasets
                        for i=1:size(currDataset(1).samples,2)
                            for j=1:length(tasks)
                                remove(j) = isnan(currDataset(j).samples(1,i));
                            end
                            
                            if sum(remove)>=1
                                removeVoxels(i) = 1;
                            else
                                removeVoxels(i) = 0;
                            end
                        end
                        
                        for i=1:length(tasks)
                            currDataset(i).samples(:,logical(removeVoxels))=[];
                            currDataset(i).fa.i(:,logical(removeVoxels))=[];
                            currDataset(i).fa.j(:,logical(removeVoxels))=[];
                            currDataset(i).fa.k(:,logical(removeVoxels))=[];
                            
                            if strcmpi(task{i},'Retrieval')==1
                                switch regressRT.flag
                                    case 'Yes'
                                        files = dir([study_path filesep subject filesep 'Run*']);
                                        for j=1:length(files)
                                            curMat(j) = load([files(j).folder filesep files(j).name]);
                                            if j==length(files)
                                                rtCell = [curMat.RT];
                                                
                                                % Convert from cell to double for regression
                                                for k=1:length(rtCell)
                                                    
                                                    % Flag outlier RT greater than 4 seconds
                                                    if double(rtCell{k}) >= regressRT.trialSec
                                                        rtDouble(k,1) = regressRT.trialSec;
                                                    else
                                                        rtDouble(k,1) = double(rtCell{k});
                                                    end
                                                end
                                                
                                                % Replace with trial duration (set in params)
                                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                                            end
                                        end
                                        
                                        % Z-Score RT values
                                        rtDouble = zscore(rtDouble);
                                        
                                        for jj=1:length(currDataset.samples)
                                            model = LinearModel.fit(rtDouble,currDataset.samples(:,j));
                                            if j==1
                                                allResiduals = model.Residuals.Raw;
                                            else
                                                allResiduals = [allResiduals model.Residuals.Raw];
                                            end
                                        end
                                        
                                        zscoreResid = zscore(allResiduals);
                                        currDataset.samples = zscoreResid;
                                        
                                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                                end
                                
                            end
                            
                            try
                                if exist('subConds','var')
                                    
                                    subConds{1,2}=[];
                                    subConds={subConds{1,1}};
                                    
                                    for ii=1:length(conds)
                                        
                                        subCond(ii).idx=contains(currDataset(i).sa.labels, ['-' conds{ii}]);
                                        counter=1;
                                        
                                        for iii=1:length(subConds)
                                            subCond(ii).idx(:,counter+1)=contains(currDataset(i).sa.labels, subConds{1,iii});
                                            counter=counter+1;
                                        end
                                        
                                        subCond(ii).idx=double(subCond(ii).idx);
                                        subCond(ii).idx(:,counter+1)=sum(subCond(ii).idx(:,1:counter),2);
                                        
                                    end
                                    
                                    Cond(1).idx = find(subCond(1).idx(:,counter+1) == counter);
                                    Cond(2).idx = find(subCond(2).idx(:,counter+1) == counter);
                                    
                                else
                                    
                                    for ii=1:length(conds)
                                        Cond(ii).labels = ~cellfun(@isempty, strfind...
                                            (currDataset(i).sa.labels, ['-' conds{ii}]));
                                        Cond(ii).idx = find(Cond(ii).labels == 1);
                                    end
                                    
                                end
                                
                            catch
                                warning('Failure to define trial information. Set to debug mode.');
                            end
                        end
                        
                        % Parse data by task and condition
                        counter=1;
                        for i=1:length(currDataset)
                            for ii=1:length(Cond)
                                if strcmpi(task{i},'Retrieval')==1
                                    retData(ii)=currDataset(i);
                                    retData(ii).samples=retData(ii).samples(Cond(ii).idx,:);
                                    retData(ii).sa.beta_index=retData(ii).sa.beta_index(Cond(ii).idx);
                                    retData(ii).sa.chunks=retData(ii).sa.chunks(Cond(ii).idx);
                                    retData(ii).sa.fname=retData(ii).sa.fname(Cond(ii).idx);
                                    retData(ii).sa.labels=retData(ii).sa.labels(Cond(ii).idx);
                                    retData(ii).sa.targets=retData(ii).sa.targets(Cond(ii).idx);
                                    %retData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                else
                                    encData(ii)=currDataset(i);
                                    encData(ii).samples=encData(ii).samples(Cond(ii).idx,:);
                                    encData(ii).sa.beta_index=encData(ii).sa.beta_index(Cond(ii).idx);
                                    encData(ii).sa.chunks=encData(ii).sa.chunks(Cond(ii).idx);
                                    encData(ii).sa.fname=encData(ii).sa.fname(Cond(ii).idx);
                                    encData(ii).sa.labels=encData(ii).sa.labels(Cond(ii).idx);
                                    encData(ii).sa.targets=encData(ii).sa.targets(Cond(ii).idx);
                                    %encData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                end
                            end
                        end
                        
                        % Setup DSM and searchlight arguments
                        %DSM
                        dsmArgs.metric = 'correlation';
                        dsmArgs.center_data = 1;
                        
                        %Searchlight
                        searchlightSize = 20;
                        metric = 'count';
                        measure = @cosmo_target_dsm_corr_measure;
                        measure_args = struct();
                        measure_args.type = 'Spearman';
                        measure_args.center_data = 1;
                        
                        for i=1:length(Cond)
                            
                            % Create Target DSM from Encoding Run
                            targetDSM(i) = cosmo_dissimilarity_matrix_measure(encData(i), dsmArgs);
                            
                            if i==1
                                % Create Searchlight
                                nbrhood=cosmo_spherical_neighborhood(retData(i),metric,searchlightSize);
                            end
                            
                            % Set target DSM
                            measure_args.target_dsm = targetDSM(i).samples;
                            
                            % Searchlight ERS for each condition separately
                            results(i) = cosmo_searchlight(retData(i),nbrhood,measure,measure_args);
                            
                        end
                        
                        save([outputPath '/searchlightResults_' metric '_' ...
                            num2str(searchlightSize) '.mat'],'results');
                        
                        for i=1:length(conds)
                            % Define output location
                            if ~exist('subConds','var')
                                output{i}=strcat(outputPath,'/',subject,...
                                    '_ERS_Searchlight_',metric,'_',...
                                    num2str(searchlightSize),'_',...
                                    conds{i},'.nii');
                            else
                                output=strcat(outputPath,'/',subject,'_',...
                                    classifier.name,'_Searchlight_',metric,'_',...
                                    num2str(searchlightSize),'_',...
                                    conds{1,1},'_vs_',conds{1,2},'_',subConds{1,1},'.nii');
                            end
                            
                            % Store results to disk
                            cosmo_map2fmri(results(i), output{i});
                            
                        end
                        
                    case 'ROI'
                        curROI = fullfile(masks(curMask).folder, masks(curMask).name);
                        %%% Loading ROI data into CosmoMVPA - use SPM betas
                        fprintf('Current subject: %s\n',subjects{iteration});
                        fprintf('Loading data from ROI: %s\n',masks(curMask).name);
                        for i=1:length(tasks)
                            
                            spmFile = ['/gpfs/group/nad12/default/nad12/ICEE/multivariate/models/SingleTrialModel' ...
                                tasks{i} '/' subjects{iteration} '/SPM.mat'];
                            
                            currDataset(i)=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                        end
                        
                        % Remove errant voxels in both datasets
                        for i=1:size(currDataset(1).samples,2)
                            for j=1:length(tasks)
                                remove(j) = isnan(currDataset(j).samples(1,i));
                            end
                            
                            if sum(remove)>=1
                                removeVoxels(i) = 1;
                            else
                                removeVoxels(i) = 0;
                            end
                        end
                        
                        %%% Indentify trials of interest for each condition
                        for i=1:length(tasks)
                            currDataset(i).samples(:,logical(removeVoxels))=[];
                            currDataset(i).fa.i(:,logical(removeVoxels))=[];
                            currDataset(i).fa.j(:,logical(removeVoxels))=[];
                            currDataset(i).fa.k(:,logical(removeVoxels))=[];
                            
                            % Mean center data
                            currDataset(i).samples = bsxfun...
                                (@minus,currDataset(i).samples,mean(currDataset(i).samples,1));
                            
                            if strcmpi(tasks{i},'Retrieval')==1
                                switch regressRT.flag
                                    case 'Yes'
                                        files = dir([study_path filesep subject filesep 'Run*']);
                                        for j=1:length(files)
                                            curMat(j) = load([files(j).folder filesep files(j).name]);
                                            if j==length(files)
                                                rtCell = [curMat.RT];
                                                
                                                % Convert from cell to double for regression
                                                for k=1:length(rtCell)
                                                    
                                                    % Flag outlier RT greater than 4 seconds
                                                    if double(rtCell{k}) >= regressRT.trialSec
                                                        rtDouble(k,1) = regressRT.trialSec;
                                                    else
                                                        rtDouble(k,1) = double(rtCell{k});
                                                    end
                                                end
                                                
                                                % Replace with trial duration (set in params)
                                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                                            end
                                        end
                                        
                                        % Z-Score RT values
                                        rtDouble = zscore(rtDouble);
                                        
                                        for jj=1:length(currDataset.samples)
                                            model = LinearModel.fit(rtDouble,currDataset.samples(:,j));
                                            if j==1
                                                allResiduals = model.Residuals.Raw;
                                            else
                                                allResiduals = [allResiduals model.Residuals.Raw];
                                            end
                                        end
                                        
                                        zscoreResid = zscore(allResiduals);
                                        currDataset.samples = zscoreResid;
                                        
                                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                                end
                                
                            end
                            
                            try
                                if exist('subConds','var')
                                    
                                    subConds{1,2}=[];
                                    subConds={subConds{1,1}};
                                    
                                    for ii=1:length(conds)
                                        
                                        subCond(ii).idx=contains(currDataset(i).sa.labels, ['-' conds{ii}]);
                                        counter=1;
                                        
                                        for iii=1:length(subConds)
                                            subCond(ii).idx(:,counter+1)=contains(currDataset(i).sa.labels, subConds{1,iii});
                                            counter=counter+1;
                                        end
                                        
                                        subCond(ii).idx=double(subCond(ii).idx);
                                        subCond(ii).idx(:,counter+1)=sum(subCond(ii).idx(:,1:counter),2);
                                        
                                    end
                                    
                                    Cond(1).idx = find(subCond(1).idx(:,counter+1) == counter);
                                    Cond(2).idx = find(subCond(2).idx(:,counter+1) == counter);
                                    
                                else
                                    
                                    currDataset(i).sa.targets = zeros(size(currDataset(i).samples,1),1);
                                    for ii=1:length(conds)
                                        Cond(ii).(tasks{i}).labels = ~cellfun(@isempty, strfind...
                                            (currDataset(i).sa.labels, conds{ii}));
                                        Cond(ii).(tasks{i}).idx = find(Cond(ii).(tasks{i}).labels == 1);
                                        currDataset(i).sa.targets(Cond(ii).(tasks{i}).idx) = ii;
                                    end
                                    
                                end
                                
                            catch
                                warning('Failure to define trial information. Set to debug mode.');
                            end
                        end
                        
                        % Parse data by task and condition
                        counter=1;
                        for i=1:length(currDataset)
                            for ii=1:length(Cond)
                                switch tasks{i}
                                    case 'Encoding'
                                        encData(ii) = currDataset(i);
                                        %encData(ii).samples=encData(ii).samples(Cond(ii).(tasks{i}).idx,:);
                                        encData(ii).samples=encData(ii).samples;
                                        encData(ii).sa.beta_index=encData(ii).sa.beta_index;
                                        encData(ii).sa.chunks=encData(ii).sa.chunks;
                                        encData(ii).sa.fname=encData(ii).sa.fname;
                                        encData(ii).sa.labels=encData(ii).sa.labels;
                                        encData(ii).sa.targets=encData(ii).sa.targets;
                                    case 'Retrieval'
                                        retData(ii) = currDataset(i);
                                        retData(ii).samples = retData(ii).samples;
                                        retData(ii).sa.beta_index = retData(ii).sa.beta_index;
                                        retData(ii).sa.chunks = retData(ii).sa.chunks;
                                        retData(ii).sa.fname = retData(ii).sa.fname;
                                        retData(ii).sa.labels = retData(ii).sa.labels;
                                        retData(ii).sa.targets = retData(ii).sa.targets;
                                end
                            end
                        end
                        
                        for i=1:length(Cond)
                            
                            for j=1:length(Cond(i).Encoding.idx)
                                for k=1:length(Cond(i).Retrieval.idx)
                                    corrVal(j,k) = corr(encData(i).samples(j,:)',retData(i).samples(k,:)','Type','Pearson');
                                end
                            end
                            
                            rho.Conds{i} = mean(mean(corrVal));
                            clear corrVal;
                        end
                        
                        clear currDataset removeVoxels encData retData
                end
        end
        
        %% Save text output of SVM Classification
        if strcmpi(analysisType,'Searchlight')==0
            % Create a tidyverse formatted table for final statistical analysis
            
            switch classType
                case 'RSA'
                    
                    % Create outptu stats table
                    %stats_table = table...
                    %    (subjectid, roiid, TrialTypeCombo, rho(1,2));
                    
                    % create subjectid and roiid columns
                    subjectid   = repmat(subjects(iteration), length(TrialTypeCombo(1)), 1);
                    roiid       = repmat({regionName}, length(TrialTypeCombo(1)), 1);
                    
                    for i=1:combinations
                        
                        conditions   = repmat({TrialTypeCombo(i)}, length(TrialTypeCombo(1)), 1);
                        
                        %Make table for within
                        if exist('averageRSAWithinForROI','var')
                            averageRSAWithinForROI = averageRSAWithinForROI;
                        else
                            averageRSAWithinForROI = str2double('NaN');
                        end
                        corrValWithin = averageRSAWithinForROI;
                        
                        if i==1
                            statsTableWithin = table(subjectid, roiid, conditions, corrValWithin, withinAvgCorr{i}, withinAvgCorr{i+1});
                            statsTableWithin.Properties.VariableNames{3}=['conditions'];
                            statsTableWithin.Properties.VariableNames{4}=['corrValWithin'];
                            statsTableWithin.Properties.VariableNames{5}=['corrValCond1'];
                            statsTableWithin.Properties.VariableNames{6}=['corrValCond2'];

                        else
                            tempTable = table(conditions, corrValWithin);
                            tempTable.Properties.VariableNames{1}=['conditions' num2str(i)];
                            tempTable.Properties.VariableNames{2}=['corrValWithin' num2str(i)];
                            statsTableWithin = [statsTableWithin tempTable];
                        end
                        
                        %Make table for between
                        if exist('averageRSABetweenForROI','var')
                            averageRSABetweenForROI = averageRSABetweenForROI;
                        else
                            averageRSABetweenForROI = str2double('NaN');
                        end
                        corrValBetween = averageRSABetweenForROI;
                        
                        if i==1
                            statsTableBetween = table(subjectid, roiid, conditions, corrValBetween);
                            statsTableBetween.Properties.VariableNames{3}=['conditions'];
                            statsTableBetween.Properties.VariableNames{4}=['corrValBetween'];
                        else
                            tempTable = table(conditions, corrValBetween);
                            tempTable.Properties.VariableNames{1}=['conditions'];
                            tempTable.Properties.VariableNames{2}=['corrValBetween'];
                            statsTableBetween = [statsTableBetween tempTable];
                        end
                        
                        %Make table for distinctiveness
                        %Compute distinctiveness score
                            distinctivenessRSAForROI1 = withinAvgCorr{i} - averageRSABetweenForROI;
                            distinctivenessRSAForROI2 = withinAvgCorr{i+1} - averageRSABetweenForROI;
                        
                        if exist('distinctivenessRSAForROI1','var')
                            distinctivenessRSAForROI1 = distinctivenessRSAForROI1;
                        else
                            distinctivenessRSAForROI1 = str2double('NaN');
                        end
                        
                        if exist('distinctivenessRSAForROI2','var')
                            distinctivenessRSAForROI2 = distinctivenessRSAForROI2;
                        else
                           distinctivenessRSAForROI2 = str2double('NaN');
                        end
                        
                        
                        distinctivenessRSAForROI1 = distinctivenessRSAForROI1;
                        distinctivenessRSAForROI2 = distinctivenessRSAForROI2;
                        
                        
                        if i==1
                            statsTableDistinctiveness = table(subjectid, roiid, conditions, distinctivenessRSAForROI1, distinctivenessRSAForROI2);
                            statsTableDistinctiveness.Properties.VariableNames{3}=['conditions' num2str(i)];
                            statsTableDistinctiveness.Properties.VariableNames{4}=['distinctivenessScore1' num2str(i)];
                            statsTableDistinctiveness.Properties.VariableNames{5}=['distinctivenessScore2' num2str(i)];

                            
                      
                        else
                            tempTable = table(conditions, distinctivenessRSAForROI);
                            tempTable.Properties.VariableNames{1}=['conditions'];
                            tempTable.Properties.VariableNames{2}=['distinctivenessScore'];
                            statsTableDistinctiveness = [statsTableDistinctiveness tempTable];
                        end
                        
                    end
                    
                    % within stats table
                    filenameWithin = sprintf('sub-%s_roiid-%s_RSAwithin_statistics-table.csv', subjectid{:}, roiid{:});
                    writetable(statsTableWithin, fullfile(outputPath, filenameWithin));
                    
                    %between stats table
                    filenameBetween = sprintf('sub-%s_roiid-%s_RSAbetween_statistics-table.csv', subjectid{:}, roiid{:});
                    writetable(statsTableBetween, fullfile(outputPath, filenameBetween));
                    
                    %distinctiveness stats table
                    filenameDistinctiveness = sprintf('sub-%s_roiid-%s_RSAdistinctiveness_statistics-table.csv', subjectid{:}, roiid{:});
                    writetable(statsTableDistinctiveness, fullfile(outputPath, filenameDistinctiveness));
                    
                case 'ERS'
                    % Create subjectid and roiid columns
                    TrialTypeCombo = {strcat(tasks{1},'_v_',tasks{2})};
                    regionName  = erase(masks(curMask).name,{'reslice_','_bilat.nii'});
                    
                    % Create final stats table for region
                    statsTable = cell(2,3+length(Cond));
                    for i=1:length(Cond)
                        if i==1
                            statsTable(1,1:3) = {'subjectid','roiid','TrialTypeCombo'};
                            statsTable(2,1:3) = {subjects{iteration}, regionName, char(TrialTypeCombo)};
                        end
                        statsTable{1,i+3} = conds{i};
                        statsTable{2,i+3} = rho.Conds{i};
                    end
                    
                    % Write output summary file
                    file = fopen([sprintf([outputPath '/'...
                        'sub-%s_roiid-%s_statistics-table.csv'], subjects{iteration}, regionName)], 'w');
                    
                    for a=1:size(statsTable,1)
                        for b=1:size(statsTable,2)
                            var = eval('statsTable{a,b}');
                            try
                                fprintf(file, '%s', var);
                            end
                            fprintf(file, ',');
                        end
                        fprintf(file, '\n');
                    end
                    fclose(file);
            end
            
            %% Create aggregate table for easy viewing
            % Create headers
            if iteration==1 && curMask==1
                
                switch classType
                    case 'RSA'
                        
                        %within
                        summaryWithin=cell(length(subjects)+1,length(masks)*2);
                        summaryWithin{1,1}='subjectid';
                        summaryWithin{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summaryWithin{1,tmpCnt}=[masks(header).name(1:end-7)...
                                '_Similarity'];
                            summaryWithin{1,tmpCnt+1}=[masks(header).name(1:end-7)...
                                '_Similarity1'];
                            summaryWithin{1,tmpCnt+2}=[masks(header).name(1:end-7)...
                                '_Similarity2'];
                            tmpCnt=tmpCnt+3;

                        end
                        
                        %between
                        summaryBetween=cell(length(subjects)+1,length(masks)*2);
                        summaryBetween{1,1}='subjectid';
                        summaryBetween{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summaryBetween{1,tmpCnt}=[masks(header).name(1:end-7)...
                                '_SimilarityBetween'];
                            tmpCnt=tmpCnt+3;
                        end
                        
                        %distinctiveness
                        summaryDistinctiveness=cell(length(subjects)+1,length(masks)*2);
                        summaryDistinctiveness{1,1}='subjectid';
                        summaryDistinctiveness{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summaryDistinctiveness{1,tmpCnt}=[masks(header).name(1:end-7)...
                                '_Distinctiveness1'];
                            summaryDistinctiveness{1,tmpCnt+1}=[masks(header).name(1:end-7)...
                                '_Distinctiveness2'];
                            tmpCnt=tmpCnt+3;
                        end
                    
                        
                    case 'ERS'
                        summary=cell(length(subjects)+1,length(masks)*2);
                        summary{1,1}='subjectid';
                        summary{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summary{1,tmpCnt}=strcat(masks(header).name(9:end-7),...
                                '_',taskInfo.Conditions{1},'_Similarity');
                            summary{1,tmpCnt+1}=strcat(masks(header).name(9:end-7),...
                                '_',taskInfo.Conditions{2},'_Similarity');
                            summary{1,tmpCnt+2}=strcat(masks(header).name(9:end-7),...
                                '_',taskInfo.Conditions{3},'_Similarity');
                            summary{1,tmpCnt+3}=strcat(masks(header).name(9:end-7),...
                                '_',taskInfo.Conditions{4},'_Similarity');
                            tmpCnt=tmpCnt+4;
                        end
                end
                
                  switch classType
                      case 'RSA'
                          
                          %within
                          summaryWithin{1,tmpCnt}=['num_' taskInfo.Conditions{1}];
                          summaryWithin{1,tmpCnt+1}=['num_' taskInfo.Conditions{2}];

                          
                          %between
                          summaryBetween{1,tmpCnt}=['num_' taskInfo.Conditions{1}];
                          summaryBetween{1,tmpCnt+1}=['num_' taskInfo.Conditions{2}];

                          
                          %distinctiveness
                          summaryDistinctiveness{1,tmpCnt}=['num_' taskInfo.Conditions{1}];
                          summaryDistinctiveness{1,tmpCnt+1}=['num_' taskInfo.Conditions{2}];

                          
                      case 'ERS'
                          summary{1,tmpCnt}=['num_' taskInfo.Conditions{1}];
                          summary{1,tmpCnt+1}=['num_' taskInfo.Conditions{2}];
                          summary{1,tmpCnt+2}=['num_' taskInfo.Conditions{3}];
                          summary{1,tmpCnt+3}=['num_' taskInfo.Conditions{4}];
                          
                  end
                  
                row=2;
                header=3;
                clear tmpCnt;
                
            end
            
            % Counter for resetting to next row. Uses remainder from divison of
            % region counter over total (e.g. 1/14) to check data should be
            % read into next subject line/row.
            iterCheck=mod(curMask,length(masks));
            
            % Add subject, trial type, and accuracy to table
            switch classType
                case 'RSA'
                    %within
                    summaryWithin{row,1}=subjects{iteration};
                    summaryWithin{row,2}=TrialTypeCombo{1,1};
              
                    
                    %between
                    summaryBetween{row,1}=subjects{iteration};
                    summaryBetween{row,2}=TrialTypeCombo{1,1};
                    
                    %distinctiveness
                    summaryDistinctiveness{row,1}=subjects{iteration};
                    summaryDistinctiveness{row,2}=TrialTypeCombo{1,1};
                    
                case 'ERS'
                    summary{row,1}=subjects{iteration};
                    summary{row,2}=TrialTypeCombo{1,1};
            end
            
            switch classType
                case 'RSA'
                    summaryWithin{row,header}=averageRSAWithinForROI;
                    summaryWithin{row,header+1}=withinAvgCorr{i};
                    summaryWithin{row,header+2}=withinAvgCorr{i+1};
                    
                    
                    summaryBetween{row,header}=averageRSABetweenForROI;
                    
                    summaryDistinctiveness{row,header}=distinctivenessRSAForROI1;
                    summaryDistinctiveness{row,header+1}=distinctivenessRSAForROI2;

                    header=header+3;
                    
                    
                case 'ERS'
                    summary{row,header}=rho.Conds{1};
                    summary{row,header+1}=rho.Conds{2};
                    summary{row,header+2}=rho.Conds{3};
                    summary{row,header+3}=rho.Conds{4};
                    header=header+4;
            end

            % Drops to next row if remainder is 0 (e.g. all regions have been
            % entered for a given subject)
            if iterCheck == 0   
                switch classType
                    case 'RSA'
                        %within
                        summaryWithin{row,header}=num2str(length(Cond(1).idx));
                        summaryWithin{row,header+1}=num2str(length(Cond(2).idx));
                        
                        %between
                        summaryBetween{row,header}=num2str(length(Cond(1).idx));
                        summaryBetween{row,header+1}=num2str(length(Cond(2).idx));
                        
                        %distinctiveness
                        summaryDistinctiveness{row,header}=num2str(length(Cond(1).idx));
                        summaryDistinctiveness{row,header+1}=num2str(length(Cond(2).idx));
                        
                    case 'ERS'
                        summary{row,header}=num2str(length(Cond(1).Retrieval.idx));
                        summary{row,header+1}=num2str(length(Cond(2).Retrieval.idx));
                        summary{row,header+2}=num2str(length(Cond(3).Retrieval.idx));
                        summary{row,header+3}=num2str(length(Cond(4).Retrieval.idx));
                end
                row=row+1;
                header=3;
            end
            
        end
    end
end

switch classType
    case 'RSA'
        save([fileparts(outputPath) filesep 'summaryWithin_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.mat'],'summaryWithin');
        save([fileparts(outputPath) filesep 'summaryBetween_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.mat'],'summaryBetween');
        save([fileparts(outputPath) filesep 'summaryDistinctiveness_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.mat'],'summaryDistinctiveness');
    case 'ERS'
        % Save mat file with statistics separately
        save([fileparts(outputPath) filesep 'summary.mat'],'summary');
end
%% Save summary files of RSA.

switch classType
    case 'RSA'
        
        if strcmpi(analysisType,'Searchlight')==0
            
            % Write output summary file
            file1 = fopen([fileparts(outputPath) filesep 'allSimilaritiesWithinSummary_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.csv'], 'w');
            
            for a=1:size(summaryWithin,1)
                for b=1:size(summaryWithin,2)
                    var = eval('summaryWithin{a,b}');
                    try
                        fprintf(file1, '%s', var);
                    end
                    fprintf(file1, ',');
                end
                fprintf(file1, '\n');
            end
        end
        
        fclose(file1);
        
        % Write output summary file
        file2 = fopen([fileparts(outputPath) filesep 'allSimilaritiesBetweenSummary_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.csv'], 'w');
        
        for a=1:size(summaryBetween,1)
            for b=1:size(summaryBetween,2)
                var = eval('summaryBetween{a,b}');
                try
                    fprintf(file2, '%s', var);
                end
                fprintf(file2, ',');
            end
            fprintf(file2, '\n');
        end
        
        fclose(file2);
        
         % Write output summary file
        file3 = fopen([fileparts(outputPath) filesep 'distinctivenessSummary_' taskInfo.Conditions{1} '_' taskInfo.Conditions{2} '_' taskInfo.Name '.csv'], 'w');
        
        for a=1:size(summaryDistinctiveness,1)
            for b=1:size(summaryDistinctiveness,2)
                var = eval('summaryDistinctiveness{a,b}');
                try
                    fprintf(file3, '%s', var);
                end
                fprintf(file3, ',');
            end
            fprintf(file3, '\n');
        end
        
        fclose(file3);
        
        clc;
        clear;
        
    case 'ERS'
        if strcmpi(analysisType,'Searchlight')==0
            
            % Write output summary file
            file = fopen([fileparts(outputPath) filesep 'allSimilaritiesSummary.csv'], 'w');
            
            for a=1:size(summary,1)
                for b=1:size(summary,2)
                    var = eval('summary{a,b}');
                    try
                        fprintf(file, '%s', var);
                    end
                    fprintf(file, ',');
                end
                fprintf(file, '\n');
            end
            
            fclose(file);
            clc;
            clear;
        end
end