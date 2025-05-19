%% Run Representational Similarity Analysis (RSA)
%   Editor:    CAN Lab
%   Updated:   4/9/2025
%
% Representational similarity analysis (RSA) for a single subject for single subject.
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

% 4/9/25 - removed all mention of ERS -- separate script now 

%Outputs within-category similarity, between-category similarity and
%distinctiveness described by Haxby 2001

% Only works for 2 conditions

% Requires zipped betas in nifti format and co-registered ROIs 

%% Set Analysis Parameters & Paths
addpath(genpath...
    ('/storage/group/nad12/default/nad12/toolbox/CoSMoMVPA-master/CoSMoMVPA-master'));

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

    %%%% For Just model folders w/Masks inside
        masks = dir([directory.Analysis filesep 'masks' ...
        filesep file(1:end-4) filesep '*.nii']); %  subjects{iteration} 
   
   %masks = dir([directory.Analysis filesep 'masks' filesep file(1:end-4)...
   %filesep subjects{iteration} filesep '*.nii']); 
          
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
                
                   
                % Trial by Trial Correlation Matrix - All Conditions --
                % Script Condition: Working (CAT CARPENTER 5/9/2025)
                %Flip Dataset
                    currDataset_flipped = currDataset.samples.';

                %Run Correlations
                    EntireMatrix=(atanh(corrcoef(currDataset_flipped)));
                        
                    %%%% You can save this out if you'd like too! %%%%
                        %Save out correlation matrices by subject and ROI in csv:
                        %csvwrite([analysis filesep subjects{iteration} filesep subjects{iteration} '_' regionName '_EntireCorrMatrix.csv'],
                        %EntireMatrix);

                %Plot 
                     % figure;
                     % imagesc(atanh(corrcoef(currDataset_flipped))); 
                     % colorbar; 
                     % colormap(jet); 
                     % caxis([-1 1]);
                     % numVars = size(currDataset_flipped, 2);  
                     % set(gca, 'XTick', 1:numVars, 'YTick', 1:numVars);
                     % set(gca, 'XTickLabel', currDataset.sa.labels, 'YTickLabel', currDataset.sa.labels);


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
                        %case 'Searchlight'  
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


                %Plot Matrix: Within Correlation matrix -- 
                    %ADDITION FOR DATA VISUALIZATION PURPOSES
                        % figure; 
                        % imagesc(atanh(corrcoef(Conditions(1).samples)));
                        % colorbar;
                        % colormap(jet);
                        % caxis([-1 1]);
                        % numVars = size(Conditions(1).samples, 2);  
                        % set(gca, 'XTick', 1:numVars, 'YTick', 1:numVars);
                        % set(gca, 'XTickLabel', Conditions(1).sa.labels, 'YTickLabel', Conditions(1).sa.labels);

                       
                      %Save out correlation matrices by subject and ROI in csv:
                        % Within1 = atanh(corrcoef(Conditions(1).samples));
                        % csvwrite([analysis filesep subjects{iteration} filesep subjects{iteration} '_' regionName '_WithinCond1Sim.csv'], Within1);


                     %ADDITION FOR DATA VISUALIZATION PURPOSES  
                        % figure;
                        % imagesc(atanh(corrcoef(Conditions(2).samples))); 
                        % colorbar; 
                        % colormap(jet); 
                        % caxis([-1 1]);
                        % numVars = size(Conditions(2).samples, 2);  
                        % set(gca, 'XTick', 1:numVars, 'YTick', 1:numVars);
                        % set(gca, 'XTickLabel', Conditions(2).sa.labels, 'YTickLabel', Conditions(2).sa.labels);
                   
                    %Save out correlation matrices by subject and ROI in csv:
                        %Within2 = atanh(corrcoef(Conditions(2).samples));
                        %csvwrite([analysis filesep subjects{iteration} filesep subjects{iteration} '_'regionName '_WithinCond2Sim.csv'], Within2);



                    %Compute between-condition correlations the Haxby way
                    betweenAvgCorr = {};
                    
                    for j = 1:size(Conditions(1).samples,2)
                        for c = 1:size(Conditions(2).samples,2)
                            betweenAvgCorr{c,j} = corr(Conditions(1).samples(:,j),Conditions(2).samples(:,c),'Type','Pearson');
                        end
                    end
                    
                    averageRSABetweenForROI = mean2(atanh(cell2mat(betweenAvgCorr)));


                    %Plot Matrix: Between Correlation matrix
                        % figure; 
                        % imagesc(atanh(cell2mat(betweenAvgCorr)));
                        % colorbar; 
                        % colormap(jet); 
                        % caxis([-1 1]);
                        % numVars = size(betweenAvgCorr, 2);  
                        % set(gca, 'XTick', 1:numVars, 'YTick', 1:numVars);
                        % set(gca, 'XTickLabel', Conditions(1).sa.labels, 'YTickLabel', Conditions(2).sa.labels);


                    %Save out correlation matrices by subject and ROI in csv:
                        %csvwrite([analysis filesep subjects{iteration} filesep subjects{iteration} '_' regionName '_BetweenCondSim.csv'],
                        %betweenAvgCorr);
                    

                    
                   % Obtain number of combinations with 2
                    % conditions selected
                    combinations=factorial(length(Cond))/...
                        ((factorial(2)*(factorial(length(Cond)-2))));
                    
                    
                    % Compute for each condition separately
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
             
        end
        
        %% Save text output of SVM Classification
        if strcmpi(analysisType,'Searchlight')==0
            % Create a tidyverse formatted table for final statistical analysis
            
            switch classType
                case 'RSA'
                    
                    % Create output stats table
                    %stats_table = table...
                    %    (subjectid, roiid, TrialTypeCombo, rho(1,2));
                    
                    % create subjectid and roiid columns
                    subjectid   = repmat(subjects(iteration), length(TrialTypeCombo(1)), 1);
                    roiid       = repmat({regionName}, length(TrialTypeCombo(1)), 1);
                    
                    for i=1:combinations
                        
                        conditions = repmat({TrialTypeCombo(i)}, length(TrialTypeCombo(1)), 1);
                        
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
        
end
