function EM_get_spm_residuals()
% TODO run without spm_filter
nvoi = 5;

for itrial = 1:150
  cd('sub-27_sd_alltrials_1-001')
  load SPM.mat
  
  spm_write_residuals_edit(SPM, 0); % residual ts comes out, 2310 vols
  
  %select VOIS of interest, load them;
  onset = round(SPM.Sess.U(1).ons); % interest
  volsoi = onset:onset+(nvoi-1);
  for ivoi = 1:nvoi
    ft_read_mri()
  end
  
  % delete all residual vols
  
  
  
  
  
  
end