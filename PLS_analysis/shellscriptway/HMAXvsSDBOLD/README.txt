===============================================================================
                       = LNDG PLS file creation = 
                   = Setup guide for our PLS scripts =
===============================================================================

-------------------------------------------------------------------------------

Table of Contents:

I. Introduction 
II. Requirements
III. SD-BOLD PLS file creation
IV. Known Issues and TODOs 
V. Contribution and feedback

-------------------------------------------------------------------------------

I. Introduction

This will be a quick guide for the creation of SD-BOLD PLS session data files necessary for PLS analysis. The standard way which we implement requires the production of several intermediary files before we have our SD-BOLD files. These SD-BOLD PLS files serve as the basis for PLS analysis using the PLS toolbox by Baycrest.

-------------------------------------------------------------------------------

II. Requirements

Some external tools are necessary for using our scripts:

- PLS toolbox from Baycrest (https://www.rotman-baycrest.on.ca/index.php?section=84)
- NIfTI toolbox (https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) 
- Also necessary is the spm_squeeze function from spm12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- Additionally, we use a custom script to resize NIfTI images, but it is not essential, as the script which calls on it also contains the equivalent code, it simply doesn't do all the checks that our code does.

-------------------------------------------------------------------------------

III. SD-BOLD PLS file creation

The file creation process happens in four steps, with an optional fifth if the user will perform mean-BOLD PLS analysis. Generally speaking, these are the steps:

- Step 0: Create folder structure
- Step 1: Create text files to create mean-BOLD files
- Step 2: Create mean-BOLD PLS files using plsgui
- Step 3: Create gray matter (GM) masked common coordiantes(CC)
- Step 4: Create SD-BOLD PLS files (using GM masked CC)
- Step 5: [Optional] Create mean-BOLD PLS files which have been GM/CC masked

Steps 2, 3, 4 & 5 can be run with via a compiled version of the .m scripts which get called by the .sh handler or by running an interactive matlab job. The matlab processes involved in PLS file creation aren't hugely resource intensive and, if the user so desires, can be called without going through the intermediary process of compilation, but this means running everything sequentially. The gain in compilation comes from running jobs in parallel on the grid (steps 2 & 4), but if the files are small enough then it's an unnecessary usage of resources. Step 3 is explicitly encouraged to be done via an interactive matlab window on the grid.

When running matlab interactively, this should be done either on the interactive node or the testing node. This avoids slowing down the master node for other users. The interactive and testing nodes are accessible via the following commands:

- qsub -I (job lasting 24 hours)
- qsub -I -q testing (job lasting 1 hour)

The bash scripts depend on a configuration file which is creatively named pls_config.sh. Within this file the user must set several parameters, which are used by the handler scripts. The file directories are set in the configuration and, if the user is following the general preproc pipeline structure, do not need to be changed. The PreprocPipe variable exists here to keep track of the preprocessing pipeline being used as an input. The following parameters must be set by the user:

- ProjectName
- PreprocPipe
- SubjectID
- SessionID
- RunID
- VoxelSize

== Step 0 ==

This is a simple script called step_0_setup_folders.sh which creates the directories necessary for creating PLS files. If the user is using alternate directories than those recommended by this pipeline, then this can be skipped.

== Step 1 ==

Step 1 uses a template file (called make_batch_template.txt) which contains the appropriate information for each condition & their respective onset/duration, as well as dummy code for both the NIfTI data file directory and for a prefix containing the subject name/identifier. The step_1_create_batchfiles.sh script loops over each participant, copies the template file and then uses the sed command to change the dummy code for the individual text files. The sed command changes the dummyDirectory and dummyID from the template to the appropriate individual subject information. This text file will serve as the batch input file to the pls toolbox for the purpose of creating the mean-BOLD PLS files. The actual format of the text file depends on the amount of images/runs/timings/conditions used for the study. As such, an explenation of how to set this up for the user's particular study is beyond the scope of this text. The user is referred to the Baycrest page and the demo files in the PLS toolbox for information pertaining to the setup of the template file. Just know that the user will have to add as new dummy code (and it's replacement in the .sh script) for each NIfTI.

== Step 2 ==

Step 2 operates through the make_meanbold_datamat.m script which feeds the previously produced text batch files to the batch_plsgui function in order to create the mean-BOLD PLS files. The script calls on the function 'batch_plsgui' with the batch text file as input. This is essentially all the script does. The script takes the SubjectID file name as an input.

== Step 3 ==

Step 3 loads all the NIfTI files for the study and creates the common coordinates of activated voxels, called common_coords in the make_commoncoords.m script. A binary MNI GM image of the appropriate voxel size is loaded and saved as GM_coords and used to mask the common_coords, which results in the final GM masked common coordinates. These masked coordinates are called final_coords in the script. This final_coords represents the voxels for which we have measured values in every NIfTI image. The user should manually set the SubjectID list in the script, even if running a compiled version.

== Step 4 ==

Step 4 uses the make_sdbold_datamat.m script to load the data structure of the mean-BOLD PLS files and replaces the mean data with the SD data. It also masks the data matrices with the final_coords created in step 3. This script is pretty well documented so going through the script while reading the comments should suffice to understand whats going on. In short: We change several parts of our original mean-BOLD PLS files to create the SD-BOLD PLS files. First we remove the st_datamat field (which we will alter to hold SD values rather than mean values) and the st_coords field (which will be replaced with our final_coords). We also change the session_info.datamat_prefix and the pls_data_path so that they now represent the SD file which we are creating. We then determine the correct number of scans, load the relevant NIfTIs, divide the data into the appropriate conditions, normalize data and calculate the SD values. These SD values get placed into the st_datamat field and then the file is saved. If the user hasn't changed the basic directories or file names, then the script runs simply by passing it the SubjectID as an input.

== Step 5 [Optional] ==

Step 5 is completely optional and only relevant if the user will be performing a mean-BOLD PLS analysis. It depends on a script called make_commoncoord_meanbold.m which preforms the masking process on the mean-BOLD PLS files using the same the final_coords from Step 3. The script changes the relevant fields, similar to how it is done in Step 4. It operates by taking the SubjectID as an input.

-------------------------------------------------------------------------------

IV. Known Issues and TODOs 

Optimally, all the text file creation gets phased out in favor of directly creating matfiles. This might be the next job on the agenda in terms of PLS scripts. There already exist some template files for this purpose but they scripts used for them still need to be generalized.

-------------------------------------------------------------------------------

V. Contribution and feedback

This is a work in progress

-------------------------------------------------------------------------------
