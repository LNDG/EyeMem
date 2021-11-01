#!/bin/bash

# Compile SPM12 with the variability toolbox added in

cd $TMPDIR

module load matlab/R2016b
matlab

restoredefaultpath
addpath('/home/mpib/kloosterman/MATLAB/tools/spm12')
addpath('/home/mpib/kloosterman/MATLAB/tools/spm12/matlabbatch')
addpath('/home/mpib/kloosterman/MATLAB/tools/spm12/config')
spm_make_standalone()