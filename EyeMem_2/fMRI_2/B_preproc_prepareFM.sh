#!/bin/bash

source preproc2_config.sh

# PBS Log Info
CurrentPreproc="prepareFM"
CurrentLog="${LogPath}/${CurrentPreproc}"

if [ ! -d ${CurrentLog} ]; then mkdir -p ${CurrentLog}; chmod 770 ${CurrentLog}; fi

#Loop over participants & sessions (if they exist)

for SUB in ${SubjectID} ; do
		
		# Path to the gre image folder.
		FieldOrigPath="${DataPath}/${SUB}/fmap"		# Path for field images
		
		# Name of phase image to be used.
		PhaseImage="${SUB}_phasediff" 	
		
		# Name of magnitudue images to be used
		MagImage1="${SUB}_magnitude1" # magnitude image corresponding to shorter TE
		MagImage2="${SUB}_magnitude2" # magnitude image corresponding to longer TE
				
		# Name of output preprocessing dir		
		FMpath="${WorkingDirectory}/data/mri/fmap/preproc/${SUB}"

		# Output path for field map prep procedure with ANTs (i.e. brain extraction)
		ANTsPath="${FMpath}/ANTs_temp"
		ANTsName="ANTs_"
		
		# Start log
		StartLog="${ANTsPath}/${ANTsName}started.txt"
		# Error message if ANTs did not produce the expected output
		CrashLog="${ANTsPath}/${ANTsName}failed.txt"
		
		# ANTs settings
		KeepTemporaryFiles="0" # don't save temporary files
		ImageDimension="3" # 3d
		
		# ANTs-specific file paths
		TemplatePath="${SharedFilesPath_standards}/ANTS/MICCAI2012-Multi-Atlas-Challenge-Data" 		# Directory for ANTs template to be used.
		TemplateImage="${TemplatePath}/T_template0.nii.gz" 											# ANTs bet template image (e.g. averaged anatomical image) - mandatory
		ProbabilityImage="${TemplatePath}/T_template0_BrainCerebellumProbabilityMask.nii.gz" 		# ANTs bet brain probability image of the template image - mandatory
		RegistrationMask="${TemplatePath}/T_template0_BrainCerebellumRegistrationMask.nii.gz" 		# ANTs bet brain mask of the template image (i.e. rough binary mask of brain location) - optional (recommended)
		
		# check if FM image exists
		if [ ! -f ${FieldOrigPath}/${PhaseImage}.nii.gz ]; then
			echo "No Fieldmap image found for ${SUB}"
			continue
		fi
		
		# check if ANTs results already produced
		if [ -f ${ANTsPath}/${ANTsName}BrainExtractionBrain.nii.gz ] || [ -f ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz ]; then 	# Verify if ANTs output was already created
			continue
		else
			if [ -f ${CrashLog} ]; then 							# Verify if crash log exists, if so, delete intermediary ANTs files and re-run ANTs.
				rm -rf ${AntsPath}				
			elif [ -f ${StartLog} ]; then 							# Verify if ANTs job started. Could be problematic if job did not finish.
				continue
			fi
		fi
		
		# Gridwise
		echo "#!/bin/bash"                    				>  job.slurm
		echo "#SBATCH --job-name ${CurrentPreproc}_${SUB}" 	>> job.slurm # Job name 
		echo "#SBATCH --time 12:00:00" 						>> job.slurm # Time until job is killed 
		echo "#SBATCH --mem 5GB" 							>> job.slurm # Books 10gb RAM for the job 
		echo "#SBATCH --output ${CurrentLog}/slurm-%j.out" 	>> job.slurm # Write output to currentlog

		echo "sleep $(( RANDOM % 120 ))"						>> job.slurm # Sleep for a random period between 1-60 seconds, used to avoid interfierence when running antsBrainExtraction.sh
		
		echo "module load fsl" 									>> job.slurm
		#echo "FSLDIR=/home/mpib/LNDG/FSL/fsl-5.0.11"  >> job
		#echo ". ${FSLDIR}/etc/fslconf/fsl.sh"                   >> job
		#echo "PATH=${FSLDIR}/bin:${PATH}"                       >> job
		#echo "export FSLDIR PATH"                               >> job
		
		echo "module load ants" 							>> job.slurm
		
		# Create temporary folder.
		echo "mkdir -p ${ANTsPath}"							>> job.slurm
		echo "chmod 770 ${ANTsPath}"						>> job.slurm
		echo "chmod 770 ${FMpath}"							>> job.slurm
		echo "echo 'ANTs will start now' >> ${StartLog}"	>> job.slurm
		
		#### From HCP Pipeline ###
		
		echo "cd ${FMpath}"									>> job.slurm
		
		#divide by two

		echo "fslmaths ${FieldOrigPath}/${MagImage1}.nii.gz -add ${FieldOrigPath}/${MagImage2}.nii.gz ${FMpath}/${SUB}_MagnitudeAdd.nii.gz" >> job.slurm
		echo "fslmaths ${FMpath}/${SUB}_MagnitudeAdd.nii.gz -div 2 ${FMpath}/${SUB}_fmap_MeanMagnitude.nii.gz" 								>> job.slurm
		echo "rm ${FMpath}/${SUB}_MagnitudeAdd.nii.gz" 		>> job.slurm

		### They say BET, I'm going to use ANTs bra

		#${FSLDIR}/bin/bet ${WD}/Magnitude ${WD}/Magnitude_brain -f 0.35 -m #Brain extract the magnitude image
		
		# Perform Brain Extraction
		
		echo -n "antsBrainExtraction.sh -d ${ImageDimension} -a ${FMpath}/${SUB}_fmap_MeanMagnitude.nii.gz -e ${TemplateImage} " 	>> job.slurm
		echo  "-m ${ProbabilityImage} -f ${RegistrationMask} -k ${KeepTemporaryFiles} -o ${ANTsPath}/${ANTsName}" 					>> job.slurm

		# If the final ANTs output isn't created, write a text file to be used as a verification of the output outcome.
		echo "if [ ! -f ${ANTsPath}/${ANTsName}BrainExtractionBrain.nii.gz ]; then echo 'BrainExtractionBrain file was not produced.' >> ${CrashLog}; exit ;fi" >> job.slurm
		
		echo "mv ${ANTsPath}/${ANTsName}BrainExtractionBrain.nii.gz ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz" >> job.slurm
		echo "rm -r ${ANTsPath}" >> job.slurm
		
		# perform bet extraction
		
		echo "bet ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz" >> job.slurm
		
		# erode boundary voxels (recommended in FSL documentation for field map preparation)
		
		echo "fslmaths ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz -ero ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz" >> job.slurm
		
		# create binary mask
		
		echo "fslmaths ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz -bin ${FMpath}/${SUB}_fmap_MeanMagnitude_brain_mask.nii.gz" >> job.slurm

# AP - binaries erode brain mask for QC purposes

#${FSLDIR}/bin/imcp ${PhaseInputName} ${WD}/Phase

		echo "fsl_prepare_fieldmap SIEMENS ${FieldOrigPath}/${PhaseImage}.nii.gz ${FMpath}/${SUB}_fmap_MeanMagnitude_brain.nii.gz ${FMpath}/${SUB}_fmap_rads 2.46"		>> job.slurm

		sbatch job.slurm
		rm job.slurm
done