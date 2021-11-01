function EM_getBOLDvar(cfg)
% compute mean and sd over trials, and MSE

infile = cfg.infile;
outfile = cfg.outfile;
PREIN = cfg.PREIN;
subjno = cfg.subjno;
removeERP = cfg.removeERP;

disp(PREIN)
disp(subjno)
tic
load(infile)
toc

%% create hmax bins
% bin based on HMAX
% sort onsets based on hmax
[sortHMAX, sortinds] = sort(sourcetrl.trialinfo(:,10));  %hmax in 10, ascending, trial inds

disp 'make bins of trials based on hmax'
nbins = 3;
[~,~,cond_bins] = histcounts(sortinds, nbins);

for ibin = 1:nbins+1
  cfg=[];
  if ibin == nbins+1
    cfg.trials = 'all';
  else
    cfg.trials = cond_bins == ibin;
  end
  source_sel = ft_selectdata(cfg, sourcetrl);  
  
  switch removeERP
    case 'yes'
      [~, nrpt] = size(source_sel.pow);
      ERP = single(mean(source_sel.pow,2));
      source_sel.pow = source_sel.pow - repmat(ERP, 1, nrpt, 1 );
      if ismac
        %%
        %         figure; plot(squeeze(mean(sourcetrl.pow(sourcetrl.inside,:,:),2)))
        %         figure; plot(squeeze(mean(source_sel.pow(source_sel.inside,:,:))))
        tmp = source_sel;
%         tmp = sourcetrl;
        tmp.pow = squeeze(mean(tmp.pow,2));
%         tmp.pow = squeeze(std(tmp.pow,0,2));
        tmp.powdimord = 'pos_time';
        tmp.dim = [tmp.dim(1:3)];
        cfg=[];
        cfg.funparameter = 'pow';
        cfg.method = 'ortho'; % slice ortho glassbrain vertex
        load colormap_jetlightgray.mat
        %     cfg.funcolormap = cmap(129:end,:);
        cfg.funcolormap = cmap;
%         cfg.funcolorlim = 'zeromax';
            cfg.funcolorlim = 'maxabs';
%         cfg.funcolorlim = [0 2]% 'zeromax'
        cfg.location = [25 11 30];
        cfg.locationcoordinates = 'voxel';
        %     cfg.colorbar = 'yes'
        %     ft_sourceplot(cfg, tmp, anat)
        ft_sourceplot(cfg, tmp);
%%
      end
  end
  
  disp 'compute std and avg across trials'
  source(1,ibin) = source_sel;
%   source(1,ibin).dim = source(1,ibin).dim([1:3, 5]); % ONLY DIMS 1:3?
  source(1,ibin).dim = source(1,ibin).dim(1:3); % ONLY DIMS 1:3?
  source(1,ibin).powdimord = 'pos_time';
  source(1,ibin).pow = single(squeeze(mean(source_sel.pow,2)));
  source(2,ibin) = source(1,ibin);
  source(2,ibin).pow = single(squeeze(std(source_sel.pow,0,2)));

  %%
  disp 'compute MSE'  
  data = [];   % make data struct  % Required fields:  %   - time, trial, label
  for itrial = 1:size(source_sel.pow, 2)
    data.time{itrial} = source_sel.time;
    data.trial{itrial} = squeeze(source_sel.pow(source_sel.inside,itrial,:));
  end
  data.label = [];
  for ichan = 1:length(data.trial{1})
    data.label{ichan,1} = sprintf('%d', ichan);
  end
  
  cfg=[];
  cfg.toi = -3:14;
  cfg.timwin = 2;
  cfg.timescales = 1;
  cfg.filtmethod = 'no';
  cfg.recompute_r = 'perscale_toi_sp';
  cfg.coarsegrainmethod = 'pointavg';
  mse = ft_entropyanalysis(cfg, data);
  if ismac
    figure; plot(mse.time,squeeze(mean(mse.sampen(:,1,:))))
    figure; plot(mse.time,squeeze(mean(mse.r)))
  end
  source(3,ibin) = source(1,ibin);
  source(3,ibin).time = mse.time;
  source(3,ibin).pow = nan(length(source_sel.inside), length(mse.time));
  source(3,ibin).pow(source(3,ibin).inside,:) = mse.sampen;
  source(3,ibin).dim = source(3).dim(1:3);
  % add r parameter
  source(4,ibin) = source(3,ibin);
  source(4,ibin).pow(source(4,ibin).inside,:) = mse.r;
  %%

end

src = [];
src.source = source;
src.dimord = 'meas_bin';
src.meas = {'avg' 'std' 'MSE' 'rpar'};
src.bin = {'hmax1' 'hmax2' 'hmax3' 'non-binned'};
clear source

%% plot mse
if ismac
  plotdat = cat(3,source(3,:).pow);
  figure; plot(source(3,1).time, squeeze(mean(plotdat(source(1).inside,:,:))))
  legend({'1' '2' '3'})
  
  close all
  for ibin=1:nbins
    tmp = source(3,ibin);
    
    cfg=[];
    cfg.funparameter = 'pow';
    cfg.method = 'ortho'; % slice ortho glassbrain vertex
    load colormap_jetlightgray.mat
    %     cfg.funcolormap = cmap(129:end,:);
    cfg.funcolormap = cmap;
    cfg.funcolorlim = 'zeromax';
    %     cfg.funcolorlim = 'maxabs';
        cfg.funcolorlim = [0 2]% 'zeromax'
            cfg.location = [25 11 30];
    cfg.locationcoordinates = 'voxel';
    %     cfg.colorbar = 'yes'
    %     ft_sourceplot(cfg, tmp, anat)
    ft_sourceplot(cfg, tmp);
  end
end

disp(outfile);
save(outfile, 'src')

%%
if ismac
  tmp_std = source_std;
  cfg=[];
  cfg.funparameter = 'pow';
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  load colormap_jetlightgray.mat
  %     cfg.funcolormap = cmap(129:end,:);
  cfg.funcolormap = cmap;
  cfg.funcolorlim = 'zeromax';
  %     cfg.funcolorlim = 'maxabs';
  %   cfg.funcolorlim = [0 150]% 'zeromax'
  cfg.funcolorlim = 'zeromax' % 'zeromax'
  %     cfg.colorbar = 'yes'
  ft_sourceplot(cfg, tmp_std)
  
  tmp_avg = source_avg;
  cfg=[];
  cfg.funparameter = 'pow';
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  load colormap_jetlightgray.mat
  %     cfg.funcolormap = cmap(129:end,:);
  cfg.funcolormap = cmap;
  cfg.funcolorlim = 'maxabs';
  %   cfg.funcolorlim = [0 150]% 'zeromax'
  %     cfg.colorbar = 'yes'
  ft_sourceplot(cfg, tmp_avg)
  
  %   tmp.powdimord = 'pos';
  %   cfg=[];
  %   cfg.method = 'ortho'; % slice ortho glassbrain vertex
  %   cfg.funparameter = 'pow';
  %   %   cfg.funcolorlim = [-300 300];
  %   %   cfg.atlas = atlas;
  %   ft_sourceplot(cfg, tmp)
  
end
