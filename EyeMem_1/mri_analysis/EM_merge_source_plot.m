function EM_merge_source_plot(source)
% plot avg and std time courses for voxels of interest
close all
%%
% baseline correction?
doBLC = 0;
regressSTD = 0;
for iage = 1:2
  for imeas = 3 [1,3,4] % 3:4%:3
    for ibin=4
      s = source.source(iage,imeas,ibin);
      s.dim = s.dim(1:3);
      if regressSTD
        sd = source.source(iage,2,4);
        sd.pow = sd.pow(:,:,6);
        cfg = [];
        cfg.confound = mean(sd.pow(:,:),2); % matrix, [Ntrials X Nconfounds], may not contain NaNs
        s = ft_regressconfound(cfg, s);
      end
      if doBLC
        cfg=[];
        cfg.latency = [-3 -2];
        cfg.avgovertime = 'yes';
        bl = ft_selectdata(cfg, s);
        s.pow(:,s.inside,:) = s.pow(:,s.inside,:) - repmat(bl.pow(:,bl.inside),1,1, length(s.time)) ;
%         cfg=[];
%         cfg.latency = [1 1];
%         s = ft_selectdata(cfg, s);
%         disp 'average over subjects'
%         s.pow = squeeze(mean(s.pow))';
%         s.dimord = 'pos';
        s.pow = squeeze(mean(s.pow));
        s.dimord = 'pos_time';
      else
        disp 'average over subjects'
        s.pow = squeeze(mean(s.pow));
        s.dimord = 'pos_time';
      end
      
      cfg=[];
      cfg.funparameter = 'pow';
      cfg.method = 'ortho'; % slice ortho glassbrain vertex
      load colormap_jetlightgray.mat
      %     cfg.funcolormap = cmap(129:end,:);
      cfg.funcolormap = cmap;
      if imeas == 1
        cfg.funcolorlim = [-1 1];
      elseif imeas == 2
        %       cfg.funcolorlim = 'zeromax';
        cfg.funcolorlim = [20 80];
      elseif imeas == 3
        cfg.funcolorlim = [0.7 1.2];
        %       cfg.funcolorlim = [0.8 1.4];
        %               cfg.funcolorlim = [-0.12 0.12];
        if doBLC
          cfg.funcolorlim = [-0.2 0.2];
        end
      elseif imeas == 4
        cfg.funcolorlim = [10 40]; % 'maxabs';
%         cfg.funcolorlim = [-10 10]; % 'maxabs'; % 'maxabs';
      end
      cfg.location = [25 11 31];
      cfg.locationcoordinates = 'voxel';
      cfg.colorbar = 'yes'; % does not work, line 1110
%       cfg.atlas = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/vtpm/vtpm.mat';
%       load(cfg.atlas)
%       cfg.roi = 'V1v'
%       s.coordsys = 'mni'; % nifti

      ft_sourceplot(cfg, s )
      title(sprintf('%s %s', source.age{iage}, source.meas{imeas} ))
      
      dat = squeeze(source.source(iage,imeas,1:3));
      dat = cat(4,dat.pow);
      dat = squeeze(mean(dat(:,sub2ind([60 72 60], 25, 11, 30),:,:)));
%       dat = squeeze(mean(dat(:,sub2ind([60 72 60], 31, 63, 30),:,:)));
      subplot(224); plot(source.source(iage,imeas,1).time, dat)
      legend({'1' '2' '3'}, 'Location', 'Southeast')
      ylim(cfg.funcolorlim)
      grid on
      
      %     subplot(2,2,4); ylim([0 100])
      %     subplot(2,2,3); hold on; colorbar
      
    end
  end
end