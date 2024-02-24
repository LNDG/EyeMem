%% plot BS ratio maps for paper

% file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/MNI152_T1_3mm_brain.nii.gz';
file = '/Volumes/LNDG/Standards/MNI152_T1_0.5mm_brain.nii.gz';
% file = '/Volumes/LNDG/Standards/MNI152_T1_3mm_brain.nii.gz';
mri_standard = ft_read_mri(file);


%% TaskPLS SDBold.
% TODO atlas?

% file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__86_80_earson_BfMRIbsr_lv1_Z>3.hdr';
file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/taskPLS/gaze-specific/SDbold_vs_HMAX_youngold_44_41_BfMRIbsr_lv1.hdr';
mri_BSratio = ft_read_mri(file);
% ft_write_mri('/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/SDbold_vs_HMAX_youngold_44_41_BfMRIbsr_lv1.nii', mri_BSratio.anatomy,  'dataformat', 'nifti')
mri_BSratio.functional = mri_BSratio.anatomy;

cfg = [];
cfg.parameter = 'functional';
mri_BSratio = ft_sourceinterpolate(cfg, mri_BSratio, mri_standard);
% mri_BSratio.mask = mri_BSratio.functional; %double(mri_BSratio.functional > -3 | mri_BSratio.functional < 3);
mri_BSratio.mask = (mri_BSratio.functional > 3 | mri_BSratio.functional < -3);

cfg = [];
cfg.method        = 'slice';   % slice ortho
cfg.nslices = 10;
cfg.slicerange = [100 250];
% cfg.slicerange = [17 40];
cfg.funparameter  = 'functional';
cfg.maskparameter = 'mask';
cfg.opacitymap    = 'rampup';
cfg.funcolormap = colormaps(2);
cfg.title = 'TaskPLS SDbold 5 bins';
cfg.funcolorlim   = [3 5];
% cfg.atlas = '/Users/kloosterman/Library/CloudStorage/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/spm_anatomy/AllAreas_v18.mat'
ft_sourceplot(cfg, mri_BSratio, mri_standard); % hack in ft_sourceplot in 833: always 2 rows

% saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.pdf')
% saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.png')
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/taskPLS.pdf')
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/taskPLS.png')

%% BehavPLS SDBold.
close all

file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__86_80_earson_BfMRIbsr_lv1.hdr';
mri_BSratio = ft_read_mri(file);
mri_BSratio.functional = mri_BSratio.anatomy;
removefields(mri_BSratio, 'anatomy');

cfg = [];
cfg.parameter = 'functional';
mri_BSratio = ft_sourceinterpolate(cfg, mri_BSratio, mri_standard);
mri_BSratio.mask = (mri_BSratio.functional > 3 | mri_BSratio.functional < -3);

cfg = [];
cfg.method        = 'slice';   % slice ortho
cfg.nslices = 10;
cfg.slicerange = [100 250];
% cfg.slicerange = [17 40];
cfg.funparameter  = 'functional';
cfg.funcolorlim   = [-6 6];
cfg.maskparameter = 'mask';
% cfg.maskstyle = 'colormix';
cfg.maskstyle     = 'opacity';

cmap_cool = flipud(colormaps(3));
cmap_hot = colormaps(2);
cmap_black = zeros(round(length(cmap_hot)*2),3) + 1;
cmap = [cmap_cool; cmap_black; cmap_hot];
cfg.funcolormap = cmap;

cfg.opacitymap    = 'auto';
cfg.opacitylim    = 'auto';
% cfg.atlas = '/Users/kloosterman/Library/CloudStorage/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/spm_anatomy/AllAreas_v18.mat'
ft_sourceplot(cfg, mri_BSratio, mri_standard); %  hack in ft_sourceplot in 833: always 2 rows

saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.pdf')
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/behavpls_Z3.png')

%% TaskPLS MeanBold.
close all

file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/mean_5bins/fixednbins/taskPLS/gaze-specific/SDbold_vs_HMAX_youngold_45_42_BfMRIbsr_lv1.hdr';
mri_BSratio = ft_read_mri(file);
mri_BSratio.functional = mri_BSratio.anatomy;
removefields(mri_BSratio, 'anatomy');

cfg = [];
cfg.parameter = 'functional';
mri_BSratio = ft_sourceinterpolate(cfg, mri_BSratio, mri_standard);
mri_BSratio.mask = (mri_BSratio.functional > 3 | mri_BSratio.functional < -3);

ft_write_cifti('/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/mean_5bins', mri_BSratio, 'parameter', 'functional')

%%
cfg = [];
% cfg.method        = 'slice';   % slice ortho
cfg.method        = 'surface';   % slice ortho
cfg.surffile      = 'surface_inflated_right_caret.mat';

cfg.nslices = 10;
cfg.slicerange = [100 250];
% cfg.slicerange = [17 40];
cfg.funparameter  = 'functional';
cfg.funcolorlim   = [-8 8];
cfg.maskparameter = 'mask';
% cfg.maskstyle = 'colormix';
cfg.maskstyle     = 'opacity';

cmap_cool = flipud(colormaps(3));
cmap_hot = colormaps(2);
cmap_black = zeros(round(length(cmap_hot)*1.25),3) + 1;
cmap = [cmap_cool; cmap_black; cmap_hot];
cfg.funcolormap = cmap;

cfg.opacitymap    = 'auto';
cfg.opacitylim    = 'auto';
% cfg.atlas = '/Users/kloosterman/Library/CloudStorage/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/spm_anatomy/AllAreas_v18.mat'
ft_sourceplot(cfg, mri_BSratio, mri_standard); %  hack in ft_sourceplot in 833: always 2 rows

view ([90 0])
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/taskpls_meanBOLD.pdf')
saveas(gcf, '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/EyeMem/plots/taskpls_meanBOLD.png')

