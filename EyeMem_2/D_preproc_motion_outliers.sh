source preproc2_config.sh

# PBS Log Info
CurrentPreproc="mo_out"
CurrentLog="${LogPath}/${CurrentPreproc}"
if [ ! -d ${CurrentLog} ]; then mkdir -p ${CurrentLog}; chmod 770 ${CurrentLog}; fi

# Error log
Error_Log="${CurrentLog}/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${CurrentLog}

# Loop over participants, sessions (if they exist) & runs/conditions/tasks/etc
for SUB in ${SubjectID} ; do
	for TASK in ${TaskID}; do
		
		if [ $TASK == "rest" ]; then RunID="NoRun"; else source preproc2_config.sh; fi
				
		for RUN in ${RunID}; do
			
			# Name of functional image to be used.
			if [ ${TASK} == "rest" ]; then
				FuncImage="${SUB}_task-${TASK}_bold"
			elif [ ${TASK} == "eyemem2" ]; then
				FuncImage="${SUB}_task-${TASK}_run-${RUN}_bold"
			fi
			
			# Path to the original functional image folder.
			OriginalPath="${DataPath}/${SUB}/func"
			
			# Path to the pipeline specific folder.
			if [ ${TASK} == "rest" ]; then
				FuncPath="${WorkingDirectory}/data/mri/resting_state/preproc/${SUB}"
			elif [ ${TASK} == "eyemem2" ]; then
				FuncPath="${WorkingDirectory}/data/mri/task/preproc/${SUB}/run-${RUN}"
			fi
				
			if [ ! -f ${OriginalPath}/${FuncImage}.nii.gz ]; then
				continue
			elif [ -f ${FuncPath}/motionout/${SUB}_motionout.txt ]; then
				continue
			fi
			
			if [ ${TASK} == "rest" ]; then
				DeleteVolumes=${DeleteVolumes_rest}
			elif [ ${TASK} == "eyemem2" ]; then
				DeleteVolumes=${DeleteVolumes_rest}
			else
				echo "Task not found"; continue
			fi
			
			# Create output path and temp dir for motion outlier detection
			if [ ! -d ${FuncPath}/motionout ]; then mkdir -p ${FuncPath}/motionout; fi
			if [ ! -d ${FuncPath}/motionout/tmp ]; then mkdir -p ${FuncPath}/motionout/tmp; fi

			# Gridwise (SLURM)
			echo "#!/bin/bash"													> job.slurm
			echo "#SBATCH --job-name ${CurrentPreproc}_${SUB}_${TASK}_${RUN}" 	>> job.slurm # Job name 
			echo "#SBATCH --time 12:00:00" 										>> job.slurm # Time until job is killed 
			echo "#SBATCH --mem 4GB" 											>> job.slurm # Books 8gb RAM for the job 
			echo "#SBATCH --output ${CurrentLog}/slurm-%j.out" 					>> job.slurm # Write (output) log to group log folder 

			echo "module load fsl"											 	>> job.slurm # Loads FSL environment
			#echo ". /etc/fsl/5.0/fsl.sh"										>> job # Set fsl environment 	
			#echo "FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11"  >> job
			#echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
			#echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
			#echo "export FSLDIR PATH"                               >> job
			
			echo "cd ${FuncPath}"           									>> job.slurm

			# Run motion outlier detection
			echo -n "fsl_motion_outliers -i ${OriginalPath}/${FuncImage} -o ${FuncPath}/motionout/${SUB}_motionout.txt -s ${FuncPath}/motionout/${SUB}_${MoutMetric}.txt " >> job.slurm  
			echo "-p ${FuncPath}/motionout/${SUB}_${MoutMetric}_plot.png --${MoutMetric} -t ${FuncPath}/motionout/tmp --dummy=${DeleteVolumes} -v >> ${FuncPath}/motionout/report.txt" >> job.slurm
			
			# Error log
			echo "if [ ! -f ${FuncPath}/motionout/${SUB}_motionout.txt ];"  	>> job.slurm
			echo "then echo 'Error in ${FuncImage}' >> ${Error_Log}; fi"		>> job.slurm
			
			sbatch job.slurm
			rm job.slurm
			
		done
	done
done