#!/bin/bash

## Noise Component Rejection

# This will remove the components which have been selected for rejection by FIX.

source ../preproc2_config.sh

FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11

OutputStage="denoised"

# Number of subjects in test set
Nsubjects=`echo ${TestSetID} | wc -w`

# Test
#SubjectID="EYEMEMtest"; FixThreshold="70"; Nsubjects="1"

# PBS Log Info
CurrentPreproc="Denoise"
CurrentLog="${LogPath}/${CurrentPreproc}"
if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; chmod 770 ${CurrentLog}; fi

# Error Log
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

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
			
			if [ ! -f hand_labels_noise.txt ]; then
				echo "no hand_labels_noise file for ${SUB} run-${RUN}"
				continue
			fi
			
			Training="hand_labels_noise.txt"
			Rejected=`cat ${Training} | tail -1`
			#if [ -f ${Training} ]
			#	cat ${Training} | tail -1  > rejcomps.txt
			#else
			#	echo "${Training} does not exist, cannot extract component list." >> ${Error_Log}
			#	continue
			#fi
			
			# Gridwise
			echo "#PBS -N ${CurrentPreproc}_${FuncImage}" 			>> job # Job name 
			echo "#PBS -l walltime=1:00:00" 						>> job # Time until job is killed 
			echo "#PBS -l mem=2gb" 									>> job # Books 4gb RAM for the job 
			#echo "#PBS -m n" 										>> job # Email notification on abort/end, use 'n' for no notification 
			echo "#PBS -o ${CurrentLog}" 							>> job # Write (output) log to group log folder 
			echo "#PBS -e ${CurrentLog}" 							>> job # Write (error) log to group log folder 

			# Initialize FSL
			# echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
			echo "FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11"  >> job
			echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
			echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
			echo "export FSLDIR PATH"                               >> job	
					
			# Variables for denoising
	
			Preproc="${FuncPath}/${FuncImage}.nii.gz"									# Preprocessed data image
			Denoised="${FuncPath}/${FuncImage}_${OutputStage}.nii.gz"								# Denoised image
			Melodic="${FuncPath}/FEAT.feat/filtered_func_data.ica/melodic_mix"						# Location of ICA generated Melodic Mix
			
																	# List of components to be removed
			
			# Run fsl_regfilt command
			echo "fsl_regfilt -i ${Preproc} -o ${Denoised} -d ${Melodic} -f \"${Rejected}\""  		>> job
			
			# Change permissions
			echo "chmod 770 ${Denoised}"  															>> job
			
			# Error Log
			echo "Difference=\`cmp ${Preproc} ${Denoised}\`" >> job
			echo "if [ -z \${Difference} ]; then echo 'Denoising did not change the preprocessing image: ${FuncImage}' >> ${Error_Log}; fi" >> job
			
			
			
			qsub job
			rm job
			
		done
	done
done
