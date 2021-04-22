#!/bin/bash

# Configuration file

# This file will set environments and variables which will be used throughout the preprocessing procedures.

#################################################################
#################################################################
#TODO: REMEMBER TO VERIFY IF THE IMAGES HAVE TO BE RESLICED!!!!!#
#################################################################
#################################################################

## Study variables

# Project name. This should be the name of the folder in which the study data is saved.
ProjectName="eyemem_tardis"

## Set the base directory. This will be the place where we can access common files which are not particular to this project, for example, MNI images, gray matter masks and any shared toolboxes.

BaseDirectory="/home/mpib/kamp"

## Set the working directory for the current project & it's script/data paths

WorkingDirectory="${BaseDirectory}/${ProjectName}"  		# Common project directory
DataPath="${WorkingDirectory}/BIDS/B_data"	 						# Root Data
ScriptsPath="${WorkingDirectory}/study_information/A_scripts/G_Git" 	# Pipe specific scripts
LogPath="${WorkingDirectory}/study_information/B_logs"			# Common log paths
SharedFilesPath_standards="${BaseDirectory}/Standards" 				# Toolboxes, Standards, etc
SharedFilesPath_toolboxes="${BaseDirectory}/toolboxes" 

if [ ! -d ${LogPath} ]; then
	mkdir -p ${LogPath}
	chmod 770 ${LogPath}
fi

# Set subject ID list. Use an explicit list. No commas.

# All IDs
#sub-11009 sub-11012 sub-11018 sub-11019 sub-11021 sub-11022 sub-11024 sub-11029 sub-11031 sub-11032 sub-11038 sub-11046 sub-11048 sub-11049 sub-11051 sub-11060 
#sub-11067 sub-11074 sub-11075 sub-11083 sub-11092 sub-11102 sub-12005 sub-12014 sub-12016 sub-12030 sub-12033 sub-12034 sub-12036 sub-12037 sub-12039 sub-12050 
#sub-12053 sub-12054 sub-12070 sub-12084 sub-12090 sub-12091 sub-12100 sub-12106 sub-21006 sub-21008 sub-21010 sub-21026 sub-21035 sub-21044 sub-21047 sub-21052 
#sub-21059 sub-21065 sub-21069 sub-21071 sub-21072 sub-21079 sub-21081 sub-21086 sub-21088 sub-21093 sub-21094 sub-21095 sub-21104 sub-21105 sub-21114 sub-21115 
#sub-21116 sub-22004 sub-22007 sub-22011 sub-22020 sub-22021 sub-22023 sub-22025 sub-22040 sub-22041 sub-22043 sub-22045 sub-22056 sub-22058 sub-22061 sub-22062 
#sub-22064 sub-22066 sub-22068 sub-22076 sub-22080 sub-22099 sub-22103

SubjectID="sub-12050 sub-22007 sub-22020 sub-22021"

# Set task ID list (e.g. rest, eyemem, eyemem2). Task names must adhere to BIDS naming!
TaskID="eyemem2" # "rest" or "eyemem" or "rest eyemem" or "eyemem2"

# Name of experimental conditions, runs or task data to be analyzed. No commas.
RunID="01 02 03 04 05 06"

# Voxel sizes & TR:
VoxelSize="3"
#TR="1"

# FEAT standard variables
ToggleMCFLIRT="1"								# 0=No, 1=Yes: Default is 1
BETFunc="1" 									# 0=No, 1=Yes: Default is 1
#TotalVolumes="600"
DeleteVolumes_rest="12"
DeleteVolumes_task="12" 			
HighpassFEAT="0"			# 0=No, 1=Yes‚
SmoothingKernel="7"
RegisterStructDOF="BBR" 						# Default is BBR, other DOF options: 12, 9, 7, 6, 3
MNIImage="${SharedFilesPath_standards}/MNI152_T1_3mm_brain"

# Secondary FEAT variables (normally unused)
NonLinearReg="0"								# 0=No, 1=Yes: Default is 0
NonLinearWarp="10"								# Default is 10, applied only if NonLinearReg=1
IntensityNormalization="0" 						# 0=No, 1=Yes: Default is 0
SliceTimingCorrection="0"						# Default is 0; 0:None,1:Regular up,2:Regular down,3:Use slice order file,4:Use slice timings file,5:Interleaved

# Other FIELDmap variables
Unwarping="1"
UnwarpDir="-y"
EpiSpacing="0.28499967" # in ms
EpiTE="30" # in ms
SignalLossThresh="10"

# Motion outlier detection metric
MoutMetric="dvars"

# Detrend variables
PolyOrder="3" 									# Default is 3

# Filter variables
HighpassFilterLowCutoff="0.01"					# Default is 0.01, can be set to "off"" if not perforing Highpass
LowpassFilterHighCutoff="off" 					# Default is 0.1, can be set to "off" if not performing Lowpass
FilterOrder="8" 								# Default is 8

# ICA variables
dimestVALUE="mdl" 								# Default is mdl
bgthresholdVALUE="3" 							# Default is 3
mmthreshVALUE="0.5" 							# Default is 0.5
dimensionalityVALUE="0" 						# Default is 0
AdditionalParameters="-v --Ostats" 				# Default are '-v --Ostats', verbose and, output thresholded maps and probability maps

# FIX
## Test Set for FIX
### Need to double check number of younger and old subjects
TestSetID="sub-11 sub-18 sub-39 sub-56 sub-62 sub-88 sub-97 sub-17 sub-44 sub-51 sub-70 sub-73 sub-78 sub-90 sub-50 sub-15 sub-29 sub-40 sub-42 sub-46 sub-75 sub-85 sub-91 sub-14 sub-28 sub-37 sub-38 sub-86 sub-89 sub-99 sub-24 sub-54 sub-10 sub-36 sub-81 sub-57 sub-74 sub-82"

## Accepted FIX Threshold
FixThreshold=""

# Additional parameters
StandardsAndMasks="${ProjectName}_Standards" #For template
MeanValue="10000" #For re-adding mean