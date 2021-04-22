#!/bin/bash

## Compiled Matlab matlab_master_do

source preproc2_config.sh

# load FSL

module load fsl 
#FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11
#. ${FSLDIR}/etc/fslconf/fsl.sh      
#PATH=${FSLDIR}/bin:${PATH}         
#export FSLDIR PATH

# Preprocessing suffix. Denotes the preprocessing stage of the data.
InputStage="feat" 		# Before doing steps on matlab
OutputStage="feat_detrended_highpassed" 		# After running Detrend & Filtering

# PBS Log Info
CurrentPreproc="Detrend_Filt"
CurrentLog="${LogPath}/${CurrentPreproc}"
if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; chmod 770 ${CurrentLog}; fi

# Error log
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${CurrentLog}

# Loop over participants, sessions (if they exist) & runs/conditions/tasks/etc
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source preproc2_config.sh; fi
				
		for RUN in ${RunID}; do
			
			# Path to the original functional image folder.
			OriginalPath="${DataPath}/${SUB}/func"
			
			# Name of original functional image.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold"
			elif [ ${TASK} == "eyemem2" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold"
			fi
			
			# Path to the pipeline specific folder.
			if [ ${TASK} == "rest" ]; then
				FuncPath="${WorkingDirectory}/data/mri/resting_state/preproc/${SUB}"
			elif [ ${TASK} == "eyemem2" ]; then
				FuncPath="${WorkingDirectory}/data/mri/task/preproc/${SUB}/run-${RUN}"
			fi
			
			if [ ! -f ${FuncPath}/FEAT.feat/filtered_func_data.nii.gz ]; then
				echo "no filtered_func_data file found for ${FuncImage}"
				continue
			elif [ -f ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz ] && [ ! -f ${FuncPath}/${FuncImage}_${InputStage}_detrended.nii.gz ] ; then
				continue
			fi
				
			# get TR from original image
			TR=`fslinfo ${OriginalPath}/${FuncImage}.nii.gz | grep pixdim4`; TR=${TR:15}

			# Gridwise (SLURM)
			echo "#!/bin/bash"										> job.slurm
			echo "#SBATCH --job-name ${CurrentPreproc}_${SUB}_${TASK}_${RUN}" 			>> job.slurm # Job name 
			echo "#SBATCH --time 4:00:00" 							>> job.slurm # Time until job is killed 
			echo "#SBATCH --mem 4gb" 								>> job.slurm # Books 4gb RAM for the job 
			echo "#SBATCH --output ${CurrentLog}/slurm-%j.out" 		>> job.slurm # Write (output) log to group log folder 
			
			echo "sleep $(( RANDOM % 120 ))"						>> job.slurm # Sleep for a random period between 1-60 seconds, used to avoid interferences            
			echo "cd ${ScriptsPath}"	 							>> job.slurm
			
			# Inputs must be given in the correct order: 
			# FuncPath, FuncImage, PolyOrder, TR, HighpassFilterLowCutoff, LowpassFilterHighCutoff, FilterOrder
			echo -n "./run_E_preproc_detrend_filter_excecute.sh /opt/matlab/R2017b/ " 																>> job.slurm 
			echo "${FuncPath} ${FuncImage}_${InputStage} ${PolyOrder} ${TR} ${HighpassFilterLowCutoff} ${LowpassFilterHighCutoff} ${FilterOrder}" 	>> job.slurm  
			
			# Cleanup
			echo "if [ -f ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz ];" 					>> job.slurm
			echo "then rm -rf ${FuncPath}/${FuncImage}_${InputStage}_detrended.nii.gz; fi" 		>> job.slurm
			
			# Error Log
			echo "if [ ! -f ${FuncPath}/${FuncImage}_${OutputStage}.nii.gz ];" 					>> job.slurm
			echo "then echo 'Error in ${FuncImage}' >> ${Error_Log}; fi"						>> job.slurm 
			
			sbatch job.slurm
			rm job.slurm
			
		done
	done
done