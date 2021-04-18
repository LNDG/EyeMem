#!/bin/bash

## Make PLS text batch files

source pls_config.sh

for SID in $SubjectID; do
    # Make individual batch files       	
	cp ${ScriptsPath}/make_batch_template.txt ${MeanPLS}/${SID}_meanbold_batch.txt
	
	# Modify file using sed command to edit the subject ID and file location/path within the batch text file.
	
	#TODO: Note that the rest of this script will have to me modified for however many sessions/runs/conditions must be specified for the project.
	
	old_1="dummyID"
	new_1="${SID}"
	old_2="dummyDirectory1"
	new_2="${ProjectDirectory}/data/mri/task/${PreprocPipe}/B_data/${SID}/fractals/SDNIFTI/${SID}_sd_fractals_${PreprocSuffix}.nii"
	old_3="dummyDirectory2"
	new_3="${ProjectDirectory}/data/mri/task/${PreprocPipe}/B_data/${SID}/landscapes/SDNIFTI/${SID}_sd_landscapes_${PreprocSuffix}.nii"
	old_4="dummyDirectory3"
	new_4="${ProjectDirectory}/data/mri/task/${PreprocPipe}/B_data/${SID}/naturals/SDNIFTI/${SID}_sd_naturals_${PreprocSuffix}.nii"
	old_5="dummyDirectory4"
	new_5="${ProjectDirectory}/data/mri/task/${PreprocPipe}/B_data/${SID}/streets1/SDNIFTI/${SID}_sd_streets1_${PreprocSuffix}.nii"
	old_6="dummyDirectory5"
	new_6="${ProjectDirectory}/data/mri/task/${PreprocPipe}/B_data/${SID}/streets2/SDNIFTI/${SID}_sd_streets2_${PreprocSuffix}.nii"
	
	# The '' is neccessary when performing on a local machine, but not on the grid 
	#sed -i '' "s|${old_1}|${new_1}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	#sed -i '' "s|${old_2}|${new_2}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	#sed -i '' "s|${old_3}|${new_3}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	#sed -i '' "s|${old_4}|${new_4}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	#sed -i '' "s|${old_5}|${new_5}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	#sed -i '' "s|${old_6}|${new_6}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	
	# The '' is neccessary when performing on a local machine, but not on the grid 
	sed -i "s|${old_1}|${new_1}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i "s|${old_2}|${new_2}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i "s|${old_3}|${new_3}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i "s|${old_4}|${new_4}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i "s|${old_5}|${new_5}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i "s|${old_6}|${new_6}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	
	# Unzip all NIfTI files
	# gunzip ${new_2}.gz
				
done