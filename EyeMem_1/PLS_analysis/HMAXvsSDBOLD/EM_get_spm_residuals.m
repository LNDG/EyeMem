function EM_get_spm_residuals(cfg)
% TODO run without spm_filter TODO check timeseries
nvoi = 5;

mripath = cfg.mripath;
disp(mripath)
cd(mripath)
trials = dir('*alltrials_*');
ntrials = length(trials);
subjstr = cfg.subjstr;
outfile = cfg.outfile;

mri = ft_read_mri(sprintf('../%s_sd_alltrials_sd.nii', subjstr));
mri.dim = [mri.dim nvoi ntrials ]; % reshape after collecting CHECK!
mri.anatomy = nan(mri.dim);

for itrial = 1:ntrials
  disp(itrial)
  cd(trials(itrial).name)
  load('SPM.mat')
  
  if ~ismac % only works on tardis
    allvols_spm = spm_write_residuals_edit(SPM, 0); % residual ts comes out, 2310 vols
  end
  
  %select VOIS of interest, load them;
  onset = round(SPM.Sess.U(1).ons); % interest
  onset=onset+5; % shift due to hrf
  volsoi = onset:onset+(nvoi-1);
  
  res_vols = dir('Res_*.nii');
  for ivoi = 1:nvoi
    volname = res_vols(volsoi(ivoi)).name;
    temp = ft_read_mri( volname );
    mri.anatomy(:,:,:,ivoi,itrial) = temp.anatomy;
  end
  
  % delete all residual vols
  if ~ismac % only  on tardis
    delete('Res_*.nii')
  end
  
  cd ..
  
end

mri.anatomy = mri.anatomy(:,:,:,:);
mri.dim = [mri.dim(1:3) prod(mri.dim(4:5))];

disp(outfile)
ft_write_mri(outfile, single(mri.anatomy),  'transform', mri.transform, 'dataformat', 'nifti'  )



% write nifti with 150*5 = 750 residuals, can go to EM_mri_to_ftsource