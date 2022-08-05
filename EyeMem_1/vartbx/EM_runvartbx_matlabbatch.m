function EM_runvartbx_matlabbatch(cfg)

disp(cfg.batchfile)
load(cfg.batchfile)

spm_jobman('run', matlabbatch)