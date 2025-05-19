# multi-FRAME-updated

This set of scripts was initially created by user dbe5007 as a part of our lab's multivariate analysis pipeline. It has since been updated and additional analyses have been added to the pipeline.

The package requires:
* MATLAB R2017b (only for runMVPAClassification.m, all other scripts may use an updated version of MATLAB).
* SPM12
* CoSMoMVPA toolbox (Also available on Github)


# Main Package Functions

```createParams.m```
Creates a parameter MAT file with all relevant information for preprocessing model specification/estimation, and multivariate analyses. Also includes querying location and registration of regions of interest (ROIs) to be used for multivariate analyses. Parameter file is named using processing choices made within the script to aid in file identification.

```preprocessData.m```
Performs 1st level preprocessing on functional data, including registration of functional images, slice scan time correction, and motion calculation.

```specifyModel.m```
Creates single trial models for functional data across all runs. Requires behavioral file to be sourced for onsets and trial tag information

```estimateModel.m```
Estimates a single trial model for every trial of interest using SPM12. The resulting model assigns a beta value for each trial to be used for multivariate analysis. Models are gzipped to save space and will not impede classification analysis. Requires SPM.

```runMVPAClassification.m```
Runs MVPA within singular ROIs or a Searchlight within mask provided by user. Currently employs SVM classifier for ROI level, and LDA for searchlight. Output is saved as tidyverse-formatted CSV file and as MAT files. Searchlight results are saved as weights NiFTi files. Summary CSV file is generated as well, concatenating all subject accuracies. Requires CoSMoMVPA

```runRSA.m```
Runs RSA within singular ROIs, or ERS analyses. Output is saved as tidyverse-formatted CSV file and as MAT files. Summary CSV file is generated as well, concatenating all subject accuracies. This version now allows for Trial-Level or Category-Level analysis. Requires CoSMoMVPA

```runRSASearchlight.m```
Runs RSA Searchlight in user provided masks. Output is saved as a summary .mat file, and a .nii file containing subject specific correlation maps. Requires CoSMoMVPA


See https://github.com/dbe5007/multi-FRAME/blob/master/README.md for additional information. 
