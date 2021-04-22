#!/bin/bash

## FIX: Feature Extraction

source ../preproc2_config.sh

# PBS Log Info
CurrentPreproc="Extract_Features"
CurrentLog="${LogPath}/FIX/${CurrentPreproc}"

# Create rejcomps folder for 3_create_hand_labels_noise
if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; fi; chmod 770 ${CurrentLog}

# Error Log
Error_Log="${LogPath}/FIX/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

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

			# !!!!! Deleting fix to test if our local FSL copy works !!!!!!!!!!
			
			rm -rf ${FuncPath}/FEAT.feat/fix
			cd ${FuncPath}/FEAT.feat
			delfix=`ls | grep "fix4*"`
			for DelID in ${delfix}l; do
				rm -rf ${DelID}
			done
			rm -rf fix4melview_Training_EyeMem_N101_LOO_thr50.txt
		done
	done
done