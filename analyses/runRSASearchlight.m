%% RSA Searchlight

% This script was created by Cat Carpenter - contact cmc84@psu.edu for any
% questions
% Last updated 7/25/23

% This script implements the cosmoMVPA toolbox and uses a dissimilarity
% matrix within the searchlight

%This is multi-FRAME compatible as of 7/26/23

% This script is currently set up to run a binary similarity searchlight, 
% wherein the similarity between conditions are weighted to 1, 
% and similarity within conditions is weighted to  0
% from here, similarity is computed on each voxel within the specified searchlight.

%%  demo from http://www.cosmomvpa.org/_static/publish/demo_fmri_searchlight_rsm.html 

%path to cosmo-mvpa
addpath(genpath...
    ('/gpfs/group/nad12/default/nad12/toolbox/CoSMoMVPA-master/CoSMoMVPA-master'));

%% Set data paths
% The function cosmo_config() returns a struct containing paths to tutorial
% data. (Alternatively the paths can be set manually without using
% cosmo_config.)
config=cosmo_config();

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
study_path = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

% reset citation list
cosmo_check_external('-tic');

%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iteration = 1:length(subjects)

  % define data filenames & load data   
  data_path=fullfile(study_path, subjects{iteration}); 
  output_path = fullfile(analysis, subjects{iteration});  
  masks = dir([directory.Analysis filesep 'masks' filesep file(1:end-4)...
   filesep subjects{iteration} filesep '*.nii']); 
  data_fn=[data_path '/SPM_gz.mat']; %For fMRIprep  

  % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end

for curMask = 1:length(masks)

curROI = fullfile(masks(curMask).folder, masks(curMask).name); 
regionName=erase(masks(curMask).name,'.nii'); 

% Set conditions
Conds = taskInfo.Conditions;

ds=cosmo_fmri_dataset([data_fn ':beta'],'mask',curROI);

%ds=cosmo_fmri_dataset([data_fn ':beta'],'mask',curROI,'targets', 1:144, 'chunks', 1);

% Clear errant Not a Number (NaN) values
    % Remove constant features
   ds=cosmo_remove_useless_data(ds);  

%% Simple RSM searchlight
% define 'simple' target structure where trials from the same condition are
% equal to 0. (distance = 0)
% all other pairs are equally dissimilar (distance=1).
curROI

for ii=1:length(Conds)                    
     Cond{1,ii}=contains(ds.sa.labels, Conds{1,ii});
     counter=1;                          
end
   
Condition(1).idx = find(Conds{1,1}(:,counter+1) == counter);
Condition(2).idx = find(Conds{1,2}(:,counter+1) == counter);

CondList = zeros(size(ds.samples,1),1);
      for ii=1:length(Conds)
          Condition(ii).labels = ~cellfun(@isempty, strfind...
              (ds.sa.labels, Conds{1,ii}));
          Condition(ii).idx = find(Condition(ii).labels == 1);
                            
         CondList(Condition(ii).idx) = ii;
                            
      end
      
     ds.sa.targets = CondList;
     
% Codes trials/conditions of no interest as 0 (see SpecifyModel script
     % for trial tag information)
     Zeroidx = find(CondList == 0);
                    
       % Removes all trials of no interest from analysis
       if isempty(Zeroidx)==0
           ds.samples(Zeroidx,:)=[];
           ds.sa.targets(Zeroidx)=[];
           ds.sa.beta_index(Zeroidx)=[];
           ds.sa.chunks(Zeroidx)=[];
           ds.sa.fname(Zeroidx)=[];
           ds.sa.labels(Zeroidx)=[];
      end     
                                      
% Split dataset into separate condition variables
      for i=1:length(Cond)
          index = find(ds.sa.targets == i);
          Conditions(i) = ds;
          Conditions(i).samples = Conditions(i).samples(index,:);
          Conditions(i).sa.beta_index = Conditions(i).sa.beta_index(index);
          Conditions(i).sa.chunks = Conditions(i).sa.chunks(index);
          Conditions(i).sa.fname = Conditions(i).sa.fname(index);
          Conditions(i).sa.labels = Conditions(i).sa.labels(index);
                                
 %Re-index trials within condition into separate targets
         Conditions(i).sa.targets = [1:length(index)]';
         Conditions(i).samples = Conditions(i).samples';                    
      end
   
nsamples=size(ds.samples,1);
target_dsm=zeros(nsamples);
target_class=ds.sa.targets;

%% Linear RSA Searchlight

% simple sanity check to ensure all attributes are set properly
cosmo_check_dataset(ds);

% print dataset
fprintf('Dataset input:\n');
cosmo_disp(ds);

%% Define feature neighorhoods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Important note! 
       %%%Your searchlight size will be
       %%%influenced by your voxel dimensions. So if your voxel is
       %%%3mm isotropic and your searchlight radius is 2 voxels, your
       %%%searchlight radius will be 6mm. But if you have 2mm isotropic
       %%%voxels, then your radius would be 4mm. Searchlight
       %%%size will need to be study specific

%   'radius', r         } either use a radius of r, or select
%   'count', c          } approximately c voxels per searchlight
%                       Notes:
%                       - These two options are mutually exclusive
%                       - When using this option for an fmri dataset, the
%                         radius r is expressed in voxel units; for an meeg
%                         source dataset, the radius r is in whatever units
%                         the source dataset uses for the positions
                    
            %can determine searchlight size by count or radius in voxels,
            %radius is more common
            
            %searchlightSize=1.5;  %Average 17.5 voxels per light
            %searchlightSize=2.5;  %Average 70.9 voxels per light
            %searchlightSize=3;  %Average 120.3 voxels per light
            
  searchlightSize = searchlight.Size;
  searchlightMetric = searchlight.Metric; 
  
% The neighborhood defined here is used three times (one for each target
% similarity matrix), so it is not recomputed for every searchlight call.
fprintf('Defining neighborhood for each feature\n');
%nbrhood=cosmo_spherical_neighborhood(ds,'radius',searchlightSize);

nbrhood=cosmo_spherical_neighborhood(ds,searchlightMetric,searchlightSize);

% print neighborhood

fprintf('Searchlight neighborhood definition:\n');
cosmo_disp(nbrhood);

%% For linear RSA Searchlight
%computes absolute difference between each pair of samples
%target_dsm=abs(bsxfun(@minus,target_class,target_class'));

%% For binary RSA Searchlight 
for row=1:nsamples
    for col=1:nsamples
        same_target_class=target_class(row)==target_class(col);

      if same_target_class
            target_dsm(row,col)= 0; %% weighted as dissimilar
        else
            target_dsm(row,col)= 1; %% weighted as similar
        end
    end
end

fprintf('Using the following target dsm\n');
disp(target_dsm);
imagesc(target_dsm)
set(gca,'XTick',1:nsamples,'XTickLabel',ds.sa.labels,...
        'YTick',1:nsamples,'YTickLabel',ds.sa.labels)

% set measure
measure=@cosmo_target_dsm_corr_measure; 
measure_args=struct();
measure_args.target_dsm=target_dsm;

% print measure and arguments
fprintf('Searchlight measure:\n');
cosmo_disp(measure);
fprintf('Searchlight measure arguments:\n');
cosmo_disp(measure_args);     

% run searchlight
ds_rsm=cosmo_searchlight(ds,nbrhood,measure,measure_args);

% Note: when these results are used for group analysis across multiple
% participants, it may be good to Fisher-transform the correlation values,
% so that they are more normally distributed. This can be done by:
ds_rsm.samples=atanh(ds_rsm.samples);

% show results
cosmo_plot_slices(ds_rsm);

% store results
output_fn=strcat(output_path,'/', 'rsm','_',...
                                    num2str(searchlightSize),'_',...
                                    regionName,'.nii');
cosmo_map2fmri(ds_rsm, output_fn)

save([sprintf([output_path '/', 'rsm','_',...
                                    num2str(searchlightSize),'_',...
                                    regionName,'.mat'])]);
 end

end