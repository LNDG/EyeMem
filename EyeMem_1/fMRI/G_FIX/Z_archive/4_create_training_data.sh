#!/bin/bash

## FIX: Training Data creation
# Create Training.Rdata from Training Set (subjects must already have hand_labels_noise.txt in their FEAT directories)

source ../preproc2_config.sh

# String with all FEAT folders

TrainingSet=""

# Error Log
CurrentPreproc="Create_Training_Data"
CurrentLog="${LogPath}/FIX"
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

## Create string with all FEAT directories
# Loop over participants, sessions (if they exist) & runs
for SUB in ${TestSetID} ; do
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
			
			if [ ! -f ${FuncPath}/${FuncImage}.nii.gz ]; then
				echo "${FuncImage}.nii.gz does not exist. Process will fail when attempting to create Rdata." >> ${Error_Log}
				continue
			fi
			
			# Training set does not have to be explicitly indicated if verifying for the hand_labels_noise file
			if [ ! -f ${FuncPath}/FEAT.feat/hand_labels_noise.txt ]; then
				echo "${FuncImage}_rejcomps.txt does not exist" >> ${Error_Log}
				continue
			fi
			
			FeatPath="${FuncPath}/FEAT.feat"
			
			TrainingSet="${TrainingSet}${FeatPath} "
			
		done
	done
done

# Number of subjects
Nsubjects=`echo ${SubjectID} | wc -w`

# Location where Training.Rdata file will be saved
cd ${ScriptsPath}/G_FIX

#. /etc/fsl/5.0/fsl.sh

#FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11
#. ${FSLDIR}/etc/fslconf/fsl.sh      
#PATH=${FSLDIR}/bin:${PATH}         
#export FSLDIR PATH
#
#module load fsl_fix
#
#fix -t Training_${ProjectName}_N${Nsubjects} -l ${TrainingSet}

if [ ! -f Training_${ProjectName}_N${Nsubjects} ]; then
	# Create Training file
	echo "Creating Training_${ProjectName}_N${Nsubjects}"
	
	echo "#PBS -N ${CurrentPreproc}_TrainingSet" 						>> job # job name 
	echo "#PBS -l walltime=48:00:00" 									>> job # time until job is killed 
	echo "#PBS -l mem=32gb" 												>> job # books 4gb RAM for the job 
	#echo "#PBS -m ae" 													>> job # email notification on abort/end   -n no notification 
	echo "#PBS -o ${CurrentLog}" 										>> job # write (error) log to group log folder 
	echo "#PBS -e ${CurrentLog}" 										>> job 
	
	
	# Initialize FSL
	#echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
	FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11
	echo "FSLDIR=${FSLDIR}"  >> job
	echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
	echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
	echo "export FSLDIR PATH"                               >> job
	
	echo "module load fsl_fix" 											>> job
	# Run Training Set
	echo "cd ${ScriptsPath}/G_FIX" 									>> job
	echo "fix -t Training_${ProjectName}_N${Nsubjects} -l ${TrainingSet}" 								>> job 
	
	# Change group permissions for all documents created using this script
	#echo "chgrp -R lip-lndg ."  										>> job
	echo "chmod -R 770 ."  												>> job
	
	qsub job  
	rm job
else
	echo "Training_${ProjectName}_N${Nsubjects} already exists. Delete existing files if you want to run again."
fi