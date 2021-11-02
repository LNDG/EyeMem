#!/bin/bash

## Noise Component Rejection

# This will remove the ICA components which have been selected manually for rejection.

source preproc2_config.sh

OutputStage="denoised"

# SLURM Log Info
CurrentLog="${LogPath}/denoise"
if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; chmod 770 ${CurrentLog}; fi

# Error Log
Error_Log="${CurrentLog}/denoise_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

# Loop over participants, runs/conditions/tasks/etc
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source ${ScriptsPath}/preproc2_config.sh; fi
				
		for RUN in ${RunID}; do
			
			# Name of functional image.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold_feat_detrended_highpassed"
			elif [ ${TASK} == "eyemem2" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed"
			fi

			if [ ${TASK} == "rest" ]; then
				FuncPath="${WorkingDirectory}/data/mri/resting_state/preproc/${SUB}"	# Path for run specific functional image
			else
				FuncPath="${WorkingDirectory}/data/mri/task/preproc/${SUB}/run-${RUN}"
			fi

			FuncName="${SUB}_task-${TASK}_run-${RUN}"

			if [ ! -f ${FuncPath}/${FuncImage}.nii.gz ]; then
				echo "No functional image found for ${SUB} run ${RUN}"
				continue
			elif [ -f ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz ]; then
				rm ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz
				echo "Deleting ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz"
				#continue
			fi

			# verify manual rejcomp file
			cd ${ScriptsPath}/rejcomps
			if [ ! -f ${FuncName}_rejcomps.txt ]; then
				echo "${FuncName}_rejcomps.txt does not exist" >> ${Error_Log}
				continue
			fi

			# Create hand_labels_noise.txt from manual rejcomp file required for filtering
			if [ -f ${FuncName}_rejcomps.txt ]; then
				echo  "${SUB} ${RUN}: creating hand_labels_noise.txt"
				if [ -f ${FuncPath}/FEAT.feat/hand_labels_noise.txt ]; then
					rm ${FuncPath}/FEAT.feat/hand_labels_noise.txt
				fi 
				cp ${FuncName}_rejcomps.txt ${FuncPath}/FEAT.feat/hand_labels_noise.txt
			fi
			
			chmod -R 770 ${FuncPath}/FEAT.feat/hand_labels_noise.txt
			
			## Remove rejected components
			cd ${FuncPath}/FEAT.feat
			Training="hand_labels_noise.txt"
			Rejected=`cat ${Training} | tail -1`
		
			# Gridwise (SLURM)
			echo "#!/bin/bash"											> job.slurm
			echo "#SBATCH --job-name Denoising_${FuncImage}"		 	>> job.slurm # Job name 
			echo "#SBATCH --time 4:00:00" 								>> job.slurm # Time until job is killed 
			echo "#SBATCH --mem 2gb" 									>> job.slurm # Books 4gb RAM for the job 
			echo "#SBATCH --output ${CurrentLog}/slurm-%j.out" 			>> job.slurm # Write (output) log to group log folder 

			# Initialize FSL
			echo "module load fsl"										>> job.slurm # Loads FSL environment
				
			# Variables for denoising
		
			Preproc="${FuncPath}/${FuncImage}.nii.gz"									# Preprocessed data image
			Denoised="${FuncPath}/${FuncImage}_${OutputStage}.nii.gz"								# Denoised image
			Melodic="${FuncPath}/FEAT.feat/filtered_func_data.ica/melodic_mix"						# Location of ICA generated Melodic Mix
						
			# Run fsl_regfilt command
			echo "fsl_regfilt -i ${Preproc} -o ${Denoised} -d ${Melodic} -f \"${Rejected}\""  		>> job.slurm
			
			# Change permissions
			echo "chmod 770 ${Denoised}"  															>> job.slurm
			
			# Error Log
			echo "Difference=$(cmp ${Preproc} ${Denoised})" >> job.slurm
			echo "if [ -z ${Difference} ]; then echo 'Denoising did not change the preprocessing image: ${FuncImage}' >> ${Error_Log}; fi" >> job.slurm
			
			
			sbatch job.slurm
			rm job.slurm
				
		done
	done
done 
