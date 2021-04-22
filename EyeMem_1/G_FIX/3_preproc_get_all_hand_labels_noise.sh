#!/bin/bash

## Noise Component Rejection

# This will remove the components which have been selected for rejection by FIX.

source ../preproc2_config.sh

# PBS Log Info
CurrentPreproc="Denoise"
CurrentLog="${LogPath}/${CurrentPreproc}"

## Create string with all FEAT directories
# Loop over participants, sessions (if they exist) & runs
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source ${ScriptsPath}/preproc2_config.sh; fi
		
		for RUN in ${RunID}; do
			
			# Name of functional image.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold_feat_detrended_highpassed"
			elif [ ${TASK} == "eyemem" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed"
			fi
			
			# functional image path
			if [ ${TASK} == "rest" ]; then
				FuncPath="${WorkingDirectory}/data/mri/resting_state/preproc/${SUB}"	# Path for run specific functional image
			else
				FuncPath="${WorkingDirectory}/data/mri/task/preproc/${SUB}/run-${RUN}"
			fi
			
			if [ ! -f ${FuncPath}/${FuncImage}.nii.gz ]; then
				continue
			elif [ -f ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz ]; then
				continue
			fi	
			
			## Create rejcomps.txt file for func/condition
			cd ${FuncPath}/FEAT.feat
			Training="hand_labels_noise.txt"
			if [ -f ${Training} ]; then
				compList=`cat ${Training} | tail -1`
				
				echo "${SUB}_task-${TASK}_run-${RUN} >> ${compList}"  >> ${ScriptsPath}/G_FIX/rejcomps/rejcomps_all.txt
				
			else
				echo "${SUB}_task-${TASK}_run-${RUN} ${Training} file does not exist, please verify" >> ${ScriptsPath}/G_FIX/rejcomps/rejcomps_all.txt
				continue
			fi
	
			
		done
	done
done
