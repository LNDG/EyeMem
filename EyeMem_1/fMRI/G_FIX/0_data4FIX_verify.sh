#!/bin/bash

## FIX: File Preparation

# Rename the original filtered_func_data.nii.gz image and replace it with the preprocessed data. Will also attempt to take reg folder from FEAT+.feat if the default reg folder doesn't contain a highres image (necessary input for FIX).
source ../preproc2_config.sh

# Create rejcomps folder for 3_create_hand_labels_noise
if [ ! -d ${ScriptsPath}/G_FIX/rejcomps ]; then mkdir ${ScriptsPath}/G_FIX/rejcomps; fi; chmod 770 ${ScriptsPath}/G_FIX/rejcomps

# Error Log
CurrentPreproc="data4FIX_verify"
CurrentLog="${LogPath}/FIX"; if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; chmod 770 ${CurrentLog}; fi
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

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
			fi

			FIX_Files="filtered_func_data.nii.gz filtered_func_data.ica/melodic_mix mc/prefiltered_func_data_mcf.par mask.nii.gz mean_func.nii.gz reg/example_func.nii.gz reg/highres.nii.gz reg/highres2example_func.mat design.fsf" 
			
			for FID in ${FIX_Files}; do
				if [ ! -f ${FuncPath}/FEAT.feat/${FID} ]; then
					echo "${FuncImage} missing ${FID}" >> ${Error_Log}	
					echo "${FuncImage} missing ${FID}"
				fi
			done
			
		done
	done
done