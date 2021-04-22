# fMRI Data Preprocessing of the EyeMem study

> Check out the wiki for the extensive documentation.
> The preprocessing of EyeMem2 was run on SLURM.  

## Requirements

#### E_preproc_detrend_filter_excecute Executable

Verify that version 9.3 (R2017b) of the MATLAB Runtime is installed.   

If the MATLAB Runtime is not installed, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.

Alternatively, download and install the Linux version of the MATLAB Runtime for R2017b 
from the following link on the MathWorks website:

    http://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
Package and Distribute in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.    




## Files to Package for Standalone 

-E_preproc_detrend_filter_excecute 
-run_E_preproc_detrend_filter_excecute.sh (shell script for temporarily setting 
                                           environment variables and executing the 
                                           application)
   -to run the shell script, type
   
       ./run_E_preproc_detrend_filter_excecute.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where version 9.3 of the MATLAB Runtime is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have version 9.3 of the MATLAB Runtime installed in 
    /mathworks/home/application/v93, run the shell script as:
    
       ./run_E_preproc_detrend_filter_excecute.sh /mathworks/home/application/v93
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_E_preproc_detrend_filter_excecute.sh /mathworks/devel/application/matlab
-MCRInstaller.zip
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime downloaded from web" link in the
    Deployment Tool.
-This readme file 

## Definitions

For information on deployment terminology, go to
http://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.

## Appendix 

A. Linux systems:
In the following directions, replace MR by the directory where MATLAB or the MATLAB 
   Runtime is installed on the target machine.

(1) Set the environment variable XAPPLRESDIR to this value:

MR/v93/X11/app-defaults


(2) If the environment variable LD_LIBRARY_PATH is undefined, set it to the following:

MR/v93/runtime/glnxa64:MR/v93/bin/glnxa64:MR/v93/sys/os/glnxa64:MR/v93/sys/opengl/lib/glnxa64

If it is defined, set it to the following:

${LD_LIBRARY_PATH}:MR/v93/runtime/glnxa64:MR/v93/bin/glnxa64:MR/v93/sys/os/glnxa64:MR/v93/sys/opengl/lib/glnxa64

    For more detailed information about setting the MATLAB Runtime paths, see Package and 
   Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, you can
              run the shell script file run_E_preproc_detrend_filter_excecute.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    



5. Documentation of Preprocessing Steps

    0.  Anatomical Image: /EYEMEM###_t1.nii.gz 
        Functional Images: EYEMEM###/mri/<condition> EYEMEM###_<condition>.nii.gz 
        Preprocessing Images: /EYEMEM###/preproc2/<condition>/EYEMEM###_<condition>.nii.gz 
        Configuration Script: preproc2_config.sh
            Important: When reviewing the IC's it was noticed that several participants had large spikes in the beginning of the timecourse, which led us to consider deleting volumes. 
            We chose to delete the first 4 volumes for resting state images and the first 12 volumes for task  images. 
    1. BET
        Script: 01_preproc2_BET.sh
        Input: EYEMEM###/mri/t1/EYEMEM###_t1.nii.gz
        Output: EYEMEM###/mri/t1/EYEMEM###_t1_brain.nii.gz
        Parameters: ?
    
    2. FEAT 
        Script: 02_preproc2_FEAT.sh 
        Input: /EYEMEM###/mri/<condition>/EYEMEM###_<condition>.nii.gz 
        Output: EYEMEM###/preproc2/<condition>/FEAT.feat 
        Parameters
            TR: 1.0 
            Volumes: Resting State: 600; Test Conditions: 474
            Deleted Volumes: 12 (task) or 4 (rest)
            Smoothing: 7
            Bandpass: OFF
            BET Input: YES
            MNI registration: 3mm 
    
    3. Detrend and Highpasss Filter
        Script: 03_preproc2_detrend_highpass_call.sh
        Detrending
            Input: ./FEAT.feat/filtered_func_data.nii.gz
            Output: ./EYEMEM###_<condition>_feat_detrended.nii.gz
            Parameters
                Polynomical Order:3 
        Highpass Filter
            Input:  ./EYEMEM###_<condition>_feat_detrended.nii.gz
            Output: ./EYEMEM###_<condition>_feat_detrended_highpassed.nii.gz 
            Parameters
                TR: 1.0
                LowCutoff: 0.01
                No High Cut Off
                Filter Order: 8 
    
    4. ICA 
        Script: 04_preproc2_ICA.sh   
        Input: ./EYEMEM###_<condition>_feat_detrended_highpassed.nii.gz 
        Output: ./FEAT.feat/filtered_func_data.ica 
        Parameters
            Dimest: mdl 
            Bgthreshold: 3 
            Report: ./FEAT.feat/filtered_func_data.ica/report.html 
            TR: 1.0 
            Mmthresh: 0.5 
            Background: ./anat2func.nii.gz 







    Re-preprocessing of Eyemem data to include unwarping in FSL, ANTS BET, new quality checking tools and filtering (only high-pass for Task)
