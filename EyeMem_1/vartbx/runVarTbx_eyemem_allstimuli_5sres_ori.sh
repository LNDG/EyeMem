#!/bin/bash

RootPath=/home/mpib/LNDG
DATADIR=${RootPath}/EyeMem/data/mri/task/variability/B_data

subjList="sub-09 sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28 sub-29 sub-30 sub-31 sub-32 sub-33 sub-34 sub-35 sub-36 sub-37 sub-38 sub-39 sub-40 sub-41 sub-42 sub-43 sub-44 sub-45 sub-46 sub-47 sub-48 sub-49 sub-50 sub-51 sub-52 sub-53 sub-54 sub-55 sub-56 sub-57 sub-58 sub-59 sub-61 sub-62 sub-64 sub-65 sub-66 sub-67 sub-68 sub-69 sub-70 sub-71 sub-72 sub-73 sub-74 sub-75 sub-76 sub-77 sub-78 sub-79 sub-80 sub-81 sub-82 sub-83 sub-84 sub-85 sub-86 sub-87 sub-88 sub-89 sub-90 sub-91 sub-92 sub-93 sub-94 sub-95 sub-96 sub-97 sub-98 sub-99 sub-100 sub-101"


for subj in $subjList; do
	echo "#PBS -N Eyemem_VarTbx_stimlui_${subj}" 					> job
	echo "#PBS -l nodes=1:ppn=1"									>> job
	echo "#PBS -l mem=8gb" 											>> job
	echo "#PBS -l walltime=120:00:00" 								>> job
	echo "#PBS -j oe" 												>> job
	echo "#PBS -o ${RootPath}/Y_logs/variability/stimuli_5sres/"	>> job
	echo "#PBS -m n"												>> job
	echo "#PBS -d ." 												>> job
	
	echo "/home/mpib/LNDG/EyeMem/tools/SPM_compiled_Marija/standalone/run_spm12.sh /opt/matlab/R2016b run ${DATADIR}/${subj}/all_stimuli_5sres/jobs/${subj}_all_stimuli_5sres_sdmodel.mat" 		>> job


qsub job
rm job
done
