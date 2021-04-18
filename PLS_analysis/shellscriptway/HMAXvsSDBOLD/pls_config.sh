#!/bin/bash

## Configuration file
# This file will set environments and variables which will be used throughout the preprocessing procedures.

#TODO: Change study parameters for your study
## Study parameters

# Project name. This should be the name of the folder in which the study data is saved.
	# EX: ProjectName="Mock_Study"
ProjectName="EyeMem"									# No default, must be set by user

# Name of the specific pipeline which is used as input.
	# EX: "new_preproc"
PreprocPipe="variability" 							# Default is preproc
PreprocSuffix="sd" 								# Suffix of input images, the preproc steps which have been applied to them. EX: FEAT_detrended_filtered_denoised_2mm_MNI
# Set subject ID list. Use an explicit list. No commas between subjects.

SubjectID="sub-09 sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28 sub-29 sub-30 sub-31 sub-32 sub-33 sub-34 sub-35 sub-36 sub-37 sub-38 sub-39 sub-40 sub-41 sub-42 sub-43 sub-44 sub-45 sub-46 sub-47 sub-48 sub-49 sub-50 sub-51 sub-52 sub-53 sub-54 sub-55 sub-56 sub-57 sub-58 sub-59 sub-61 sub-62 sub-64 sub-65 sub-66 sub-67 sub-68 sub-69 sub-70 sub-71 sub-72 sub-73 sub-74 sub-75 sub-76 sub-77 sub-78 sub-79 sub-80 sub-81 sub-82 sub-83 sub-84 sub-85 sub-86 sub-87 sub-88 sub-89 sub-90 sub-91 sub-92 sub-93 sub-94 sub-95 sub-96 sub-97 sub-98 sub-99 sub-100 sub-101"

# No default, must be set by user

# Set session ID list. Leave as an empty string if no sessions in data path. No commas if containing any session information.
	# EX: SessionID="ses1 ses2"
SessionID="" 									# Default is empty string

# Name of experimental conditions, runs or task data to be analyzed. No commas between runs/conditions.
	# EX: RunID="task1_run1 task1_run2 task2 restingstate"
RunID="HMAXvsSDBOLD" 					# No default, must be set by user

# Voxel size for MNI image
VoxelSize="3" 									# No default, must be set by user

#TODO: These next values do not need modifying
## Set directories

# Base directory.
# CurrentDirectory=`pwd`
# UserName=`whoami`
# if [[ ${CurrentDirectory} == "/home/mpib/LNDG"* ]]; then
BaseDirectory="/Users/kloosterman/gridmaster2012/LNDG/"
# elif [[ ${CurrentDirectory} == "/home/mpib/${UserName}"* ]]; then
# 	BaseDirectory="/home/mpib/${UserName}"
# fi

## Project directories
ProjectDirectory="${BaseDirectory}/${ProjectName}/data/PLS/HMAXvsSDBOLD"  				# Base project directory
ScriptsPath="/Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/PLS_analysis/HMAXvsSDBOLD" 	# Pipe specific scripts
PLSDirectory="${BaseDirectory}/${ProjectName}/data/PLS/HMAXvsSDBOLD"					# PLS directory
LogPath="${BaseDirectory}/${ProjectName}/data/PLS/HMAXvsSDBOLD/logs"					# Common log paths for PLS processing
MeanPLS="/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/PLS/HMAXvsSDBOLD/mean_datamats" 							# Location of Mean PLS files &  batch text files
SDPLS="/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/PLS/HMAXvsSDBOLD/sd_datamats" 								# Location of SD PLS files