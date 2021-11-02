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
			elif [ ${TASK} == "eyemem2" ]; then
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
			elif [ ${TASK} == "eyemem2" ]; then
				Func2MNI="${SUB}_task-${TASK}_run-${RUN}_bold_feat_detrended_highpassed_denoised_${CurrentPreproc}"
			fi
			
			if [ -f ${FuncPath}/${Func2MNI}.nii.gz ]; then
				echo "${Func2MNI} already produced. Skipping..."
				continue
			elif [ ! -f ${FuncPath}/${FuncImage}.nii.gz ]; then
				echo "${FuncImage} not found" >> ${Error_Log}
				continue
			fi
			
			# Gridwise (SLURM)
			echo "#!/bin/bash"													> job.slurm
			echo "#SBATCH --job-name ${CurrentPreproc}_${FuncImage}" 			>> job.slurm # Job name 
			echo "#SBATCH --time 1:00:00" 										>> job.slurm # Time until job is killed 
			echo "#SBATCH --mem 4GB" 											>> job.slurm # Books 4gb RAM for the job 
			echo "#SBATCH --output ${CurrentLog}/slurm-%j.out" 					>> job.slurm # Write (output) log to group log folder 

			# Initialize FSL
			echo "module load fsl"											 	>> job.slurm # Loads FSL environment

			# echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
			#echo "FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11"  >> job
			#echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
			#echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
			#echo "export FSLDIR PATH"                               >> job
			
			echo "flirt -in ${FuncPath}/${FuncImage} -ref ${MNIImage} -applyxfm -init ${FuncPath}/FEAT.feat/reg/example_func2standard.mat -out ${FuncPath}/${Func2MNI}" >> job.slurm
			
			sbatch job.slurm
			rm job.slurm
			
		done
	done
done