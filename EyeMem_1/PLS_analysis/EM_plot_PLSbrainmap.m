function EM_plot_PLSbrainmap(pls, sesdat, source)
% plot brain maps for a behav PLS analysis

result = pls.result;

disp 'prepare for plotting'
SAV = 1;

LVno = 1;
BSRthresh = [-2 2];

disp 'reshape '
stat= [];
stat.stat   = squeeze(zeros(sesdat.st_dims));
stat.inside = squeeze(false(sesdat.st_dims));
stat.inside(sesdat.st_coords) = true;

stat.stat(stat.inside) = result.boot_result.compare_u(:,LVno);
stat.pos = source.pos; % TODO generate here so source not needed to load?
stat.dim = size(stat.stat);

disp 'add BSR mask'
stat.mask = zeros(size(stat.inside));
stat.mask(stat.inside) =  result.boot_result.compare_u(:,LVno) < BSRthresh(1) | result.boot_result.compare_u(:,LVno) > BSRthresh(2) ;
% stat.mask = stat.mask*0.5;  % TODO opacity not working yet
stat.maskdimord = 'pos';

disp 'get anatomy'
standardsfolder = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/A_standards';
anat = ft_read_mri(fullfile(standardsfolder, 'MNI152_T1_3mm_brain.nii.gz' ));

disp 'stick stats to anatomy'
cfg=[];
cfg.parameter = {'stat' 'mask'};
[interp] = ft_sourceinterpolate(cfg, stat, anat);

disp 'TODO opacity not working yet'
% interp.mask(isnan(interp.mask)) = 0; % 
% % % % TODO make mask seethrough
% % % % vol2plot.statmask(statmask == 1) = 0.1;
% % % % vol2plot.statmask(isnan(statmask(vol2plot.inside))) = 1;
% % % % vol2plot.statmaskdimord = 'pos';
% % % 

%%
disp 'plotting'
% close all
load colormap_jetlightgray.mat
cfg=[];
cfg.method = 'surface'; % slice ortho glassbrain vertex surface
cfg.funparameter = 'stat';
cfg.maskparameter = 'mask';
cfg.anaparameter = 'anatomy';
cfg.funcolormap = cmap;
cfg.surffile = 'surface_white_right.mat'; %surface_white_left
ft_sourceplot(cfg, interp)
title(sprintf('LV%d, p = %g', LVno, result.perm_result.sprob(LVno)))

SAV = 1;
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/plots';
if SAV
  outfile = sprintf('slicemap.png')
  saveas(gcf, fullfile(PREOUT, outfile ))
  cd(PREOUT)
end


%%
% title(sprintf('%s, LV%d, p = %g', whichgroup, LVno, result.perm_result.sprob(LVno)))
% hold on
% % view ([90 0]) 
% % figure %subplot(1,2,2)
% % ft_sourceplot(cfg, vol2plot )
% % view ([-90 0]) 
% if SAV
%   saveas(gcf, fullfile(PREOUT, sprintf('%s_%s_LV%d_%dbins.pdf', whichgroup, cfg.method, LVno, num_cond )))
%   saveas(gcf, fullfile(PREOUT, sprintf('%s_%s_LV%d_%dbins.png', whichgroup, cfg.method, LVno, num_cond )))
%   cd(PREOUT)
% end
