#!/bin/bash

# Register to MNI

source preproc2_config.sh

# PBS Log Info
CurrentPreproc="MNI"
CurrentLog="${LogPath}/${CurrentPreproc}"
if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; chmod 770 ${CurrentLog}; fi
	
# Error log
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${CurrentLog}

FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11

# Loop over participants, sessions (if they exist) & runs/conditions/tasks/etc
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source ${ScriptsPath}/preproc2_config.sh; fi
				
		for RUN in ${RunID}; do
			
			# Name of functional image.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold_feat_detrended_highpassed_denoised"
			elif [ ${TASK} == "eyemem" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed_denoised"
			fi
			
			# Path to the functional image
			if [ ${TASK} == "rest" ]; then
				FuncPath="${WorkingDirectory}/data/mri/resting_state/preproc/${SUB}"	# Path for run specific functional image
			else
				FuncPath="${WorkingDirectory}/data/mri/task/preproc/${SUB}/run-${RUN}"
			fi
			
			# Name of the output image
			if [ ${TASK} == "rest" ]; then
				Func2MNI="${SUB}_task-${TASK}_bold_feat_detrended_highpassed_denoised_${CurrentPreproc}"
			elif [ ${TASK} == "eyemem" ]; then
				Func2MNI="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed_denoised_${CurrentPreproc}"
			fi
			
			if [ -f ${FuncPath}/${Func2MNI}.nii.gz ]; then
				echo "${Func2MNI} already produced. Skipping..."
				continue
			elif [ ! -f ${FuncPath}/${FuncImage}.nii.gz ]; then
				echo "${FuncImage} not found" >> ${Error_Log}
				continue
			fi
			
			# Gridwise
			echo "#PBS -N ${CurrentPreproc}_${FuncImage}" 				>> job # Job name 
			echo "#PBS -l walltime=1:00:00" 							>> job # Time until job is killed 
			echo "#PBS -l mem=4gb" 										>> job # Books 4gb RAM for the job 
			echo "#PBS -m n" 											>> job # Email notification on abort/end, use 'n' for no notification 
			echo "#PBS -o ${CurrentLog}" 								>> job # Write (output) log to group log folder 
			echo "#PBS -e ${CurrentLog}" 								>> job # Write (error) log to group log folder 

			# Initialize FSL
			# echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
			echo "FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11"  >> job
			echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
			echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
			echo "export FSLDIR PATH"                               >> job
			
			echo "flirt -in ${FuncPath}/${FuncImage} -ref ${MNIImage} -applyxfm -init ${FuncPath}/FEAT.feat/reg/example_func2standard.mat -out ${FuncPath}/${Func2MNI}" >> job
			
			qsub job
			rm job
			
		done
	done
done