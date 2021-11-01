#!/bin/bash

## Make PLS text batch files

source pls_config.sh

for SID in $SubjectID; do
    # Make individual batch files       	
	cp ${ScriptsPath}/make_batch_template.txt ${MeanPLS}/${SID}_meanbold_batch.txt
	
	# Modify file using sed command to edit the subject ID and file location/path within the batch text file.
	
	#TODO: Note that the rest of this script will have to me modified for however many sessions/runs/conditions must be specified for the project.
	
	old="dummyID"
	old_1="dummyDirectory1"
	old_2="dummyDirectory2"
	old_3="dummyDirectory3"
	old_4="dummyDirectory4"
	old_5="dummyDirectory5"
	old_6="dummyDirectory6"
	old_7="dummyDirectory7"
	old_8="dummyDirectory8"
	old_9="dummyDirectory9"
	old_10="dummyDirectory10"
	old_11="dummyDirectory11"
	old_12="dummyDirectory12"
	old_13="dummyDirectory13"
	old_14="dummyDirectory14"
	old_15="dummyDirectory15"
	old_16="dummyDirectory16"
	old_17="dummyDirectory17"
	old_18="dummyDirectory18"
	old_19="dummyDirectory19"
	old_20="dummyDirectory20"
	old_21="dummyDirectory21"
	old_22="dummyDirectory22"
	old_23="dummyDirectory23"
	old_24="dummyDirectory24"
	old_25="dummyDirectory25"
	old_26="dummyDirectory26"
	old_27="dummyDirectory27"
	old_28="dummyDirectory28"
	old_29="dummyDirectory29"
	old_30="dummyDirectory30"
	
	new="${SID}"
	new_1="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial5_beta_series_SD.nii"
	new_2="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial10_beta_series_SD.nii"
	new_3="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial15_beta_series_SD.nii"
	new_4="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial20_beta_series_SD.nii"
	new_5="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial25_beta_series_SD.nii"
	new_6="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial30_beta_series_SD.nii"
	new_7="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial35_beta_series_SD.nii"
	new_8="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial40_beta_series_SD.nii"
	new_9="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial45_beta_series_SD.nii"
	new_10="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial50_beta_series_SD.nii"
	new_11="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial55_beta_series_SD.nii"
	new_12="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial60_beta_series_SD.nii"
	new_13="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial65_beta_series_SD.nii"
	new_14="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial70_beta_series_SD.nii"
	new_15="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial75_beta_series_SD.nii"
	new_16="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial80_beta_series_SD.nii"
	new_17="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial85_beta_series_SD.nii"
	new_18="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial90_beta_series_SD.nii"
	new_19="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial95_beta_series_SD.nii"
	new_20="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial100_beta_series_SD.nii"
	new_21="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial105_beta_series_SD.nii"
	new_22="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial110_beta_series_SD.nii"
	new_23="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial115_beta_series_SD.nii"
	new_24="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial120_beta_series_SD.nii"
	new_25="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial125_beta_series_SD.nii"
	new_26="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial130_beta_series_SD.nii"
	new_27="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial135_beta_series_SD.nii"
	new_28="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial140_beta_series_SD.nii"
	new_29="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial145_beta_series_SD.nii"
	new_30="${ProjectDirectory}/SDNIFTI/${SID}_sd_trial150_beta_series_SD.nii"
	
	
	# The '' is neccessary when performing on a local machine, but not on the grid 
	sed -i '' "s|${old}|${new}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_1}|${new_1}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_2}|${new_2}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_3}|${new_3}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_4}|${new_4}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_5}|${new_5}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_6}|${new_6}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_7}|${new_7}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_8}|${new_8}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_9}|${new_9}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_10}|${new_10}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_11}|${new_11}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_12}|${new_12}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_13}|${new_13}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_14}|${new_14}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_15}|${new_15}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_16}|${new_16}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_17}|${new_17}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_18}|${new_18}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_19}|${new_19}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_20}|${new_20}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_21}|${new_21}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_22}|${new_22}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_23}|${new_23}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_24}|${new_24}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_25}|${new_25}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_26}|${new_26}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_27}|${new_27}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_28}|${new_28}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_29}|${new_29}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	sed -i '' "s|${old_30}|${new_30}|g" ${MeanPLS}/${SID}_meanbold_batch.txt

	
	# The '' is neccessary when performing on a local machine, but not on the grid 
	# sed -i "s|${old_1}|${new_1}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	# sed -i "s|${old_2}|${new_2}|g" ${MeanPLS}/${SID}_meanbold_batch.txt
	# sed -i "s|${old_3}|${new_3}|g" ${MeanPLS}/${SID}_meanbold_batch.txt

	
	# Unzip all NIfTI files
	# gunzip ${new_2}.gz
				
done