This folder contains scripts that were used to create 3D (x,y,z, without time dimension) nifti files containing BOLD SD data.

BOLD SD nifits were created for various conditions, e.g. remembered trials, ITIs, experimental conditions, ...
Each of these conditions has its own pair of scripts that produce 2 matlabbatches for input to the variability toolbox in SPM.
The pairs of scripts always are:

create_1stlevel_model_xxx.m
create_varTBX_model_xxx.m

The first script specifies a 1st level GLM in SPM for each participant, where infomation on events and their durations is specified.
The second script specifies parameters for the variability toolbox on how to proceed with output from the first level analysis. In this case, the lss (least-squares separate) method is chosen, so that each event is modelled separately and then SDs are taken over beta values for each event.