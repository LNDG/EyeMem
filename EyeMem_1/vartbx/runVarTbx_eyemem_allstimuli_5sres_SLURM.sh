#!/bin/bash

RootPath=/home/mpib/kloosterman
DATADIR=${RootPath}/projectdata/eyemem/

subjList="sub-09 sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28 sub-29 sub-30 sub-31 sub-32 sub-33 sub-34 sub-35 sub-36 sub-37 sub-38 sub-39 sub-40 sub-41 sub-42 sub-43 sub-44 sub-45 sub-46 sub-47 sub-48 sub-49 sub-50 sub-51 sub-52 sub-53 sub-54 sub-55 sub-56 sub-57 sub-58 sub-59 sub-61 sub-62 sub-64 sub-65 sub-66 sub-67 sub-68 sub-69 sub-70 sub-71 sub-72 sub-73 sub-74 sub-75 sub-76 sub-77 sub-78 sub-79 sub-80 sub-81 sub-82 sub-83 sub-84 sub-85 sub-86 sub-87 sub-88 sub-89 sub-90 sub-91 sub-92 sub-93 sub-94 sub-95 sub-96 sub-97 sub-98 sub-99 sub-100 sub-101"
#subjList="sub-09"


for subj in $subjList; do
#     echo "#PBS -N Eyemem_VarTbx_${subj}"                     > job
#     echo "#PBS -l nodes=1:ppn=1"                                    >> job
#     echo "#PBS -l mem=20gb"                                             >> job
#     echo "#PBS -l walltime=130:00:00"                                 >> job
#     echo "#PBS -j oe"                                                 >> job
#     echo "#PBS -o /home/mpib/kloosterman/qsub"    >> job
#     echo "#PBS -m n"                                                >> job
#     echo "#PBS -d ."                                                 >> job
#
#     echo "/home/mpib/kloosterman/MATLAB/tools/spm12/../standalone/run_spm12.sh /opt/matlab/R2016b run ${DATADIR}/GLM_TRwise/jobs/${subj}_GLM_TRwise_sdmodel.mat"         >> job
# #  echo "matlab"
# #  echo "spm_jobman('run', '${DATADIR}/GLM_TRwise/jobs/${subj}_GLM_TRwise_sdmodel.mat')"
#
# qsub job
# rm job

    #######
    # SLURM way
    echo '#!/bin/bash'                    > job.slurm
    echo "#SBATCH --job-name main_$subj"  >> job.slurm
    echo "#SBATCH --cpus-per-task 4"     >> job.slurm
    echo "#SBATCH --time 130:0:0"     >> job.slurm
    echo "#SBATCH --partition long" >> job.slurm
    echo "#SBATCH --mem 32GB"     >> job.slurm
    echo "#SBATCH --workdir /home/mpib/kloosterman/qsub" >> job.slurm
    echo "/home/mpib/kloosterman/MATLAB/tools/spm12/../standalone/run_spm12.sh /opt/matlab/R2016b run ${DATADIR}/5TRspertrial/jobs/${subj}_5TRspertrial_sdmodel.mat"  >> job.slurm
    sbatch job.slurm
    rm -f job.slurm

done

# /home/mpib/kloosterman/MATLAB/tools/spm12/../standalone/run_spm12.sh /opt/matlab/R2016b run /home/beegfs/kloosterman/projectdata/eyemem/5TRspertrial/jobs/sub-09_singletrial_5sres_sdmodel.mat