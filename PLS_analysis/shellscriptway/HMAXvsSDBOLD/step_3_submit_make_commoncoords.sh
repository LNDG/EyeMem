#!/bin/bash

## Create Common Coordinates file

source pls_config.sh

echo "#PBS -N ${ProjectName}_common_coords" 	>> job # job name 
echo "#PBS -l walltime=4:00:00" 				>> job # time until job is killed
echo "#PBS -l mem=4gb" 							>> job # books 10gb RAM for the job --> 1 CPU hat normalerweise 4gb
#echo "#PBS -l nodes=1:ppn=16" 					>> job 
#echo "#PBS -m ae" 								>> job # email notification on abort/end   -n no notification 
echo "#PBS -o ${LogPath}"						>> job # write (error) log to log folder 
echo "#PBS -e ${LogPath}" 						>> job 

echo "cd ${ScriptsPath}" 						>> job

echo "${ScriptsPath}/run_make_commoncoords.sh /opt/matlab/R2014b ${VoxelSize}" 	>> job 

echo "chmod -R 770 coords_EVAL.mat" 								>> job
 
# submit job
qsub job  
rm job
