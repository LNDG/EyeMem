function EM_merge_source_plotRvsMSEcorr(source)
% plot avg and std time courses for voxels of interest
close all

for iage = 1:2
%   for imeas = 3 [1,3,4] % 3:4%:3
    for ibin=4
      disp 'correlate MSE and r para'
      sMSE = source.source(iage,3,ibin);
      sR = source.source(iage,4,ibin);
      s = sR;
      s.pow = corr_col(sMSE.pow(:,:), sR.pow(:,:));
      s.pow = reshape(s.pow', numel(s.inside), []);
      
%       [~,~,s.pow] = REGRESS(sMSE.pow(:,:), [sR.pow(:,:) zeros(size(sR.pow(:,:)));
      
      cfg=[];
      cfg.funparameter = 'pow';
      cfg.method = 'ortho'; % slice ortho glassbrain vertex
      load colormap_jetlightgray.mat
      %     cfg.funcolormap = cmap(129:end,:);
      cfg.funcolormap = cmap;
      cfg.location = [25 11 30]; % visual
%       cfg.location = [30 55 30]; % ACC
      cfg.locationcoordinates = 'voxel';
      cfg.colorbar = 'yes'; % does not work, line 1110
      %       cfg.atlas = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/vtpm/vtpm.mat';
      %       load(cfg.atlas)
      %       cfg.roi = 'V1v'
      %       s.coordsys = 'mni'; % nifti
      cfg.funcolorlim = [-1 0];

      ft_sourceplot(cfg, s )
      title(sprintf('%s %s', source.age{iage}, 'corr MSE vs rpara across subjects' ))
      
      tind = s.time == 1;
      dat1 = sMSE.pow(:, sub2ind([60 72 60], cfg.location(1), cfg.location(2), cfg.location(3) ),tind); % at t=0
      dat2 = sR.pow(:, sub2ind([60 72 60], cfg.location(1), cfg.location(2), cfg.location(3) ),tind);

      figure; scatter(dat1, dat2)
      xlabel('MSE'); ylabel('rpara'); lsline; axis square; box on
      title(sprintf('%s corr= %g at t=%d', source.age{iage}, corr(dat1, dat2), s.time(tind) ))

%       legend({'1' '2' '3'}, 'Location', 'Southeast')
%       ylim(cfg.funcolorlim)
%       grid on
      
      %     subplot(2,2,4); ylim([0 100])
      %     subplot(2,2,3); hold on; colorbar
      
    end
%   end
end