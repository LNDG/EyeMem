function EM_getBOLDvar_LOTO(cfg)
% compute mean and sd over trials, and MSE using leave-one-trial-out
% approach

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

% %% create hmax bins
% % bin based on HMAX
% % sort onsets based on hmax
% [sortHMAX, sortinds] = sort(sourcetrl.trialinfo(:,10));  %hmax in 10, ascending, trial inds
%
% disp 'make bins of trials based on hmax'
% nbins = 3;
% [~,~,cond_bins] = histcounts(sortinds, nbins);
%
% for ibin = 1:nbins+1
% cfg=[];
% if ibin == nbins+1
%   cfg.trials = 'all';
% else
%   cfg.trials = cond_bins == ibin;
% end
% source_sel = ft_selectdata(cfg, sourcetrl);

switch removeERP
  case 'yes'
    ntrials = size(sourcetrl.trialinfo,1);
    ERP = single(mean(sourcetrl.pow,2));
    sourcetrl.pow = sourcetrl.pow - repmat(ERP, 1, ntrials, 1 );
    if ismac
      try
      %%
      %         figure; plot(squeeze(mean(sourcetrl.pow(sourcetrl.inside,:,:),2)))
      %         figure; plot(squeeze(mean(sourcetrl.pow(sourcetrl.inside,:,:))))
      tmp = sourcetrl;

      tmp = source(1);
      tmp.pow = squeeze(nanmean(tmp.pow));
      tmp.pow = tmp.pow(:,3);
      %         tmp.pow = squeeze(std(tmp.pow,0,2));
%       tmp.powdimord = 'pos_time';
      tmp.powdimord = 'pos';
      tmp.dim = [tmp.dim(1:3)];
      cfg=[];
      cfg.funparameter = 'pow';
      cfg.method = 'ortho'; % slice ortho glassbrain vertex
      load colormap_jetlightgray.mat
      %     cfg.funcolormap = cmap(129:end,:);
      cfg.funcolormap = cmap;
      %         cfg.funcolorlim = 'zeromax';
      %       cfg.funcolorlim = 'maxabs';
      cfg.funcolorlim = [-0.005 0.005];% 'zeromax'
      cfg.location = [25 11 30];
      cfg.locationcoordinates = 'voxel';
      %     cfg.colorbar = 'yes'
      %     ft_sourceplot(cfg, tmp, anat)
      ft_sourceplot(cfg, tmp);
      catch
      end
      %%
    end
end

%%
disp 'compute single trial MSE using LOTO'
data = [];   % make data struct  % Required fields:  %   - time, trial, label
ntrials = size(sourcetrl.pow, 2);
for itrial = 1:ntrials
  data.time{itrial} = sourcetrl.time;
  data.trial{itrial} = squeeze(sourcetrl.pow(sourcetrl.inside,itrial,:));
end
data.label = [];
for ichan = 1:length(data.trial{1})
  data.label{ichan,1} = sprintf('%d', ichan);
end

for itrial = 0:1%ntrials
  cfg=[];
  if itrial == 0
    cfg.trials = 'all';
  else
    trlinc = true(ntrials,1);
    trlinc(itrial) = false;
    cfg.trials = trlinc; % 0 1 1 1 etc
  end  
  cfg.toi = -1:6; %-3:14
  cfg.timwin = 2;
  cfg.timescales = 1;
  cfg.filtmethod = 'no';
  cfg.recompute_r = 'perscale_toi_sp';
  cfg.coarsegrainmethod = 'pointavg';
  mse = ft_entropyanalysis(cfg, data);
  if itrial == 0
    msealltrl = mse;
    mseLOTO = mse;
    mseLOTO.dimord = 'rpt_chan_timescales_time';
    mseLOTO.sampen = NaN([ntrials size(mse.sampen)]); %to keep single trials
    mseLOTO.r = NaN([ntrials size(mse.sampen)]); %to keep single trials
  else
    mseLOTO.sampen(itrial,:,:,:) = mse.sampen - msealltrl.sampen;
    mseLOTO.r(itrial,:,1,:) = mse.r - msealltrl.r; % TODO fix dimord rpara
  end
end

if ismac
  figure; plot(mse.time,squeeze(mean(mse.sampen(:,1,:))))
  figure; plot(mse.time,squeeze(mean(mse.r)))
end
source(1) = sourcetrl;
source(1).time = mse.time;
source(1).pow = nan(ntrials, length(sourcetrl.inside), length(mse.time));
source(1).pow(:,source(1).inside,:) = squeeze(mseLOTO.sampen);
source(1).dim = sourcetrl.dim(1:3);
source(1).powdimord = 'rpt_pos_time';
% add r parameter
source(2) = source(1);
source(2).pow(:,source(2).inside,:) = mseLOTO.r;
%%

% end

src = [];
src.source = source;
src.dimord = 'meas';
src.meas = {'MSE' 'rpar'};
disp(outfile);
save(outfile, 'src', '-v7.3')


%% 
% %% plot mse
% if ismac
%   plotdat = cat(3,source(3,:).pow);
%   figure; plot(source(3,1).time, squeeze(mean(plotdat(source(1).inside,:,:))))
%   legend({'1' '2' '3'})
%   
%   close all
%   for ibin=1:nbins
%     tmp = source(3,ibin);
%     
%     cfg=[];
%     cfg.funparameter = 'pow';
%     cfg.method = 'ortho'; % slice ortho glassbrain vertex
%     load colormap_jetlightgray.mat
%     %     cfg.funcolormap = cmap(129:end,:);
%     cfg.funcolormap = cmap;
%     cfg.funcolorlim = 'zeromax';
%     %     cfg.funcolorlim = 'maxabs';
%     cfg.funcolorlim = [0 2]% 'zeromax'
%     cfg.location = [25 11 30];
%     cfg.locationcoordinates = 'voxel';
%     %     cfg.colorbar = 'yes'
%     %     ft_sourceplot(cfg, tmp, anat)
%     ft_sourceplot(cfg, tmp);
%   end
% end
% 
% %%
% if ismac
%   tmp_std = source_std;
%   cfg=[];
%   cfg.funparameter = 'pow';
%   cfg.method = 'ortho'; % slice ortho glassbrain vertex
%   load colormap_jetlightgray.mat
%   %     cfg.funcolormap = cmap(129:end,:);
%   cfg.funcolormap = cmap;
%   cfg.funcolorlim = 'zeromax';
%   %     cfg.funcolorlim = 'maxabs';
%   %   cfg.funcolorlim = [0 150]% 'zeromax'
%   cfg.funcolorlim = 'zeromax' % 'zeromax'
%   %     cfg.colorbar = 'yes'
%   ft_sourceplot(cfg, tmp_std)
%   
%   tmp_avg = source_avg;
%   cfg=[];
%   cfg.funparameter = 'pow';
%   cfg.method = 'ortho'; % slice ortho glassbrain vertex
%   load colormap_jetlightgray.mat
%   %     cfg.funcolormap = cmap(129:end,:);
%   cfg.funcolormap = cmap;
%   cfg.funcolorlim = 'maxabs';
%   %   cfg.funcolorlim = [0 150]% 'zeromax'
%   %     cfg.colorbar = 'yes'
%   ft_sourceplot(cfg, tmp_avg)
%   
% end
