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
			elif [ -d ${FuncPath}/FEAT.feat/fix ]; then
				cd ${FuncPath}/FEAT.feat/fix
				Log=`grep "End of Matlab Script" logMatlab.txt | tail -1` # Get line containing our desired text output
				if [ ! "$Log" == "End of Matlab Script" ]; then
					echo "${SUB} ${RUN}: fix was incomplete, deleting and re-running"
					rm -rf ${FuncPath}/FEAT.feat/fix
				else
					continue
				fi
			fi

			cd ${FuncPath}/FEAT.feat
			if [ ! -f filtered_func_data.nii.gz ]; then
				echo "${SUB} ${RUN}: filtered_func_data.nii.gz missing" >> ${Error_Log}
				continue
			elif [ ! -d filtered_func_data.ica ]; then
				echo "${SUB} ${RUN}: filtered_func_data.ica missing" >> ${Error_Log}
				continue	
			elif [ ! -f mc/prefiltered_func_data_mcf.par ]; then
				echo "${SUB} ${RUN}: prefiltered_func_data_mcf.par missing" >> ${Error_Log}
				continue
			elif [ ! -f mask.nii.gz ]; then
				echo "${SUB} ${RUN}: mask.nii.gz missing" >> ${Error_Log}
				continue
			elif [ ! -f mean_func.nii.gz ]; then
				echo "${SUB} ${RUN}: mean_func.nii.gz missing" >> ${Error_Log}
				continue
			elif [ ! -f reg/example_func.nii.gz ]; then
				echo "${SUB} ${RUN}: reg/example_func.nii.gz missing" >> ${Error_Log}
				continue
			elif [ ! -f reg/highres.nii.gz ]; then
				echo "${SUB} ${RUN}: reg/highres.nii.gz missing" >> ${Error_Log}
				continue
			elif [ ! -f reg/highres2example_func.mat ]; then
				echo "${SUB} ${RUN}: reg/highres2example_func.mat missing" >> ${Error_Log}
				continue
			elif [ ! -f design.fsf ]; then
				echo "${SUB} ${RUN}: design.fsf missing" >> ${Error_Log}
				continue																								
			fi

			echo "#PBS -N ${CurrentPreproc}_${SUB}_${RUN}" 						>> job # job name 
			echo "#PBS -l walltime=3:00:00" 									>> job # time until job is killed 
			echo "#PBS -l mem=6gb" 												>> job # books 4gb RAM for the job 
			#echo "#PBS -m ae" 													>> job # email notification on abort/end   -n no notification 
			echo "#PBS -o ${CurrentLog}" 										>> job # write (error) log to group log folder 
			echo "#PBS -e ${CurrentLog}" 										>> job 
			
			
			# Initialize FSL
			#echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
			FSLDIR="/home/mpib/LNDG/FSL/fsl-5.0.11"
			echo "FSLDIR=${FSLDIR}" 								>> job
			echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
			echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
			echo "export FSLDIR PATH"                               >> job
	
			echo "module load fsl_fix" 											>> job
			# Extract features for FIX; FEAT directories with all required files and folders have to be provided
			echo "cd ${FuncPath}"												>> job
			echo "fix -f ${FuncPath}/FEAT.feat" 								>> job 

			# Change group permissions for all documents created using this script
			echo "cd ${FuncPath}/FEAT.feat/" 									>> job
			#echo "chgrp -R lip-lndg ."  										>> job
			echo "chmod -R 770 ."  												>> job

			qsub job  
			rm job
			
		done
	done
done