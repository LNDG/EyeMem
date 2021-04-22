#!/bin/bash

## FIX: Feature Extraction

source ../preproc2_config.sh

# Error Log
CurrentPreproc="Log_Matlab"
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
			
			# Verify Matlab log file
			
			cd ${FuncPath}/FEAT.feat/fix
			Log=`grep "End of Matlab Script" logMatlab.txt | tail -1` # Get line containing our desired text output
			
			if [ ! -d ${FuncPath}/FEAT.feat/fix ]; then
				echo "${SUB} ${RUN}: missing fix folder" >> ${Error_Log}
				continue
			elif [ ! "$Log" == "End of Matlab Script" ]; then
				echo "${SUB} ${RUN}: fix did not terminate properly" >> ${Error_Log}
			fi

		done
	done
done