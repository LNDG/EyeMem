#!/bin/bash

## FIX: hand_labels_file.txt creation
# FIX expects the ICA folders to be called 'filtered_func_data.ica', so either do this to rename all your folders appropriately or change all the scripts in FIX...

#TODO: User must first manually create the <subjID>_<run/condition>_rejcomps.txt files for all subjects in the Test Set. These text files must be locatedf in the ${ScriptsPath}/FIX/rejcomps directory.

source ../preproc2_config.sh

# Error Log
CurrentPreproc="Create_Hand_Labels_Noise"
CurrentLog="${LogPath}/FIX"
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

# Loop over participants, sessions (if they exist) & runs
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source ${ScriptsPath}/preproc2_config.sh; fi
				
		for RUN in ${RunID}; do
			
			# Name of functional image.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold_feat_detrended_highpassed"
				FuncName="${SUB}_task-${TASK}"
			elif [ ${TASK} == "eyemem" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed"
				FuncName="${SUB}_task-${TASK}_run-${RUN}"
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

			# Training set does not have to be explicitly indicated if verifying for the rejcomps file
			cd ${ScriptsPath}/G_FIX/rejcomps
			if [ ! -f ${FuncName}_rejcomps.txt ]; then
				echo "${FuncName}_rejcomps.txt does not exist" >> ${Error_Log}
				continue
			fi
			
			# Create hand_labels_noise.txt file
			if [ -f ${FuncName}_rejcomps.txt ]; then
				echo  "${SUB} ${RUN}: creating hand_labels_noise.txt"
				cp ${FuncName}_rejcomps.txt ${FuncPath}/FEAT.feat/hand_labels_noise.txt
			elif [ -f ${FuncPath}/FEAT.feat/filetered_func_data.ica/rejcomps ]; then 
				echo  "${SUB} ${RUN}: creating hand_labels_noise.txt"
				cat ${FuncPath}/FEAT.feat/rejcomps | tail -1  > ${FuncPath}/FEAT.feat/hand_labels_noise.txt
			fi

			# Error Log
			if [ ! -f ${FuncPath}/FEAT.feat/hand_labels_noise.txt ]; then
				echo "${FuncPath}/FEAT.feat/hand_labels_noise.txt was not created" >> ${Error_Log}
			fi
		
			#chgrp -R lip-lndg ${FuncPath}/FEAT.feat/hand_labels_noise.txt
			chmod -R 770 ${FuncPath}/FEAT.feat/hand_labels_noise.txt
		
		done
	done
done