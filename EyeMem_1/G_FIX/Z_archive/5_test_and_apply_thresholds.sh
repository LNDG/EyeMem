#!/bin/bash

## FIX: Test & Apply Thresholding Values
# Apply R.data with certain threshold on your Test Set; can be run locally (doesn't need much time)
# Threshold: the higher, the more components will be rejected; recommended: 20; we have used 50-60 

source ../preproc2_config.sh

# Evaluation Subjects, comment out when threshold has been determined and you wish to create labels for Test Data
	# EX: SubjectID="SUB001 SUB002"
#SubjectID="sub-15 sub-50 sub-24 sub-54"; Evaluation="Y"

SubjectID="sub-94"

Evaluation="Y"

# Number of subjects in test set
Nsubjects=`echo ${TestSetID} | wc -w`

# PBS Log Info
CurrentPreproc="Apply_Threshold"
CurrentLog="${LogPath}/FIX/${CurrentPreproc}"

if [ ! -d ${CurrentLog} ]; then mkdir ${CurrentLog}; fi; chmod 770 ${CurrentLog}

# Error Log
Error_Log="${LogPath}/FIX/${CurrentPreproc}_error_summary.txt"; echo "" >> ${Error_Log}; chmod 770 ${Error_Log}

# loop over old participants
for SUB in ${SubjectID} ; do
	
	# Gridwise
	echo "#PBS -N ${CurrentPreproc}_${SUB}" 				>> job # Job name 
	echo "#PBS -l walltime=1:00:00" 						>> job # Time until job is killed 
	echo "#PBS -l mem=1gb" 									>> job # Books 4gb RAM for the job 
	#echo "#PBS -m n" 										>> job # Email notification on abort/end, use 'n' for no notification 
	echo "#PBS -o ${CurrentLog}" 							>> job # Write (output) log to group log folder 
	echo "#PBS -e ${CurrentLog}" 							>> job # Write (error) log to group log folder 
	
	#echo "sleep $(( RANDOM % 60 ))"						>> job
	
	#echo ". /etc/fsl/5.0/fsl.sh"								>> job # Set fsl environment 	
	FSLDIR="/home/mpib/LNDG/FSL/fsl-5.0.11"
	echo "FSLDIR=${FSLDIR}" 								>> job
	echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
	echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
	echo "export FSLDIR PATH"                               >> job

	echo "module load fsl_fix" 											>> job
	
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
			
			# Number of subjects
			Nsubjects=`echo ${SubjectID} | wc -w`
			
			# Test Thresholds
			Training_Object="${ScriptsPath}/G_FIX/Training_${ProjectName}_N${Nsubjects}.RData"
			Training_File="fix4melview_Training_${ProjectName}_N${Nsubjects}_thr${FixThreshold}.txt"
			if [ ${Evaluation} = "Y" ]; then 
				
				# Output: txt-file (fix4melview_Training_thr##.txt) with names of noise components
				# TODO: Once a satisfactory threshold has been determined, set the FixThreshold variable in the configuration file.
				# TODO: Make sure to comment out the Evaluation Subjects when applying threhold.
				
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 20" >> job
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 30" >> job
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 40" >> job
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 50" >> job
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 60" >> job
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} 70" >> job
				
			# Apply Threshold
			else 
				echo "fix -c ${FuncPath}/FEAT.feat ${Training_Object} ${FixThreshold}" 					>> job
				# Error Log
				echo "if [ ! -f ${FuncPath}/FEAT.feat/${Training_File} ];" 								>> job
				echo "then echo '${FuncImage}: ${Training_File} does not exist' >> ${Error_Log}; fi " 	>> job

			fi
			
		done
	done
	
	qsub job
	rm job
	
done
