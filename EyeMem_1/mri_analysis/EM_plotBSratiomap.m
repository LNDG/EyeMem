
file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__86_80_earson_BfMRIbsr_lv1_Z>3.hdr';
mri_BSratio = ft_read_mri(file);
% mri_BSratio.mask = double(mri_BSratio.anatomy < 2.33 | mri_BSratio.anatomy > -2.33);
mri_BSratio.functional = mri_BSratio.anatomy;

% file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/MNI152_T1_3mm_brain.nii.gz';
file = '/Volumes/LNDG/Standards/MNI152_T1_0.5mm_brain.nii.gz';
mri_standard = ft_read_mri(file);


%% plot TODO also plot the negative cluster
cfg = [];
cfg.method        = 'slice';   % slice ortho
cfg.nslices = 15;
cfg.slicerange = [100 300];
cfg.funparameter  = 'functional';
cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [0 max(mri_BSratio.functional(:))];
cfg.funcolorlim   = [0 7];
cfg.opacitylim    = [0 2];
cfg.opacitymap    = 'rampup';
% cfg.atlas = '/Users/kloosterman/Library/CloudStorage/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/spm_anatomy/AllAreas_v18.mat'
ft_sourceplot(cfg, mri_BSratio, mri_standard)

saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.pdf')
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.png')
