%% plot BS ratio maps for paper
% file = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/MNI152_T1_3mm_brain.nii.gz';
file = '/Volumes/LNDG/Standards/MNI152_T1_0.5mm_brain.nii.gz';
% file = '/Volumes/LNDG/Standards/MNI152_T1_3mm_brain.nii.gz';
mri_standard = ft_read_mri(file);


%% Set up
% TODO look into surface
% close all
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/';
files = {
  'std_1bins/fixednbins/taskPLS/non-gazespecific/SDbold_vs_HMAX_youngold_45_43_BfMRI'
  'mean_1bins/fixednbins/taskPLS/non-gazespecific/SDbold_vs_HMAX_youngold_45_43_BfMRI' % flip
  'std_5bins/fixednbins/taskPLS/non-gazespecific/SDbold_vs_HMAX_youngold_44_43_BfMRI'
  'std_5bins/fixednbins/taskPLS/gaze-specific/SDbold_vs_HMAX_youngold_44_41_BfMRI'
  'mean_5bins/fixednbins/taskPLS/gaze-specific/SDbold_vs_HMAX_youngold_45_42_BfMRI'     % flip
  'mean_5bins/fixednbins/taskPLS/non-gazespecific/SDbold_vs_HMAX_youngold_45_43_BfMRI'  % flip
  'std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__86_80_earson_BfMRI'
  };
flip_lv_sign = [1 -1 1 1 -1 -1 1]; % -1 means flip it

%% plot brain pics
method = {'slice'}; % slice  surface
BSthresh = [-2.33 2.33];
% BSthresh = [0 0];
minclustersize = 25; % nvoxels in func space

mkdir(fullfile(PREOUTplot, 'pdf'))
mkdir(fullfile(PREOUTplot, 'png'))
mkdir(fullfile(PREOUTplot, 'epsc2'))
for ifile = 1:length(files)
  load( fullfile(PREIN, [files{ifile} 'result.mat']), 'result'); % for pval
  mri_BSratio = ft_read_mri( fullfile(PREIN, [files{ifile} 'bsr_lv1.hdr']));
  mri_BSratio.functional = mri_BSratio.anatomy .* flip_lv_sign(ifile);   removefields(mri_BSratio, 'anatomy');
  % remove clusters smaller than minclustersize
  t=bwlabeln((mri_BSratio.functional > BSthresh(2) | mri_BSratio.functional < BSthresh(1))); % remove clusters < 25 voxels
  props=regionprops(t, 'PixelIdxList');
  clus_sizes = arrayfun(@(x) length(x.PixelIdxList), props);
  voxels_inds = vertcat(props(clus_sizes < minclustersize).PixelIdxList);
  mri_BSratio.functional(voxels_inds) = 0;
  
  for im = 1:length(method)
    % plotting
    hemi = {'left' 'right'};
    for ih= 1:2
      if strcmp(method, 'slice')
        cfg = [];
        cfg.parameter = 'functional';
        mri_BSratio = ft_sourceinterpolate(cfg, mri_BSratio, mri_standard);
        mri_BSratio.mask = (mri_BSratio.functional > BSthresh(2) | mri_BSratio.functional < BSthresh(1));
      elseif strcmp(method, 'surface')
        cfg = [];
        cfg.parameter = 'functional';
        %         cfg.interpmethod = 'smudge'; % sphere_avg takes long
        %         cfg.sphereradius = 1;
        %         load('surface_inflated_left.mat')
        mri_BSratio = ft_sourceinterpolate(cfg, mri_BSratio, mri_standard);
        mri_BSratio.mask = (mri_BSratio.functional > BSthresh(2) | mri_BSratio.functional < BSthresh(1));
      end
      
      cfg = [];
      cfg.method        = method{im};   % slice ortho
      cfg.surffile      = sprintf('surface_white_%s.mat', hemi{ih});
      %   cfg.surffile      = sprintf('surface_pial_%s.mat', hemi{ih});
      cfg.surfinflated      = sprintf('surface_inflated_%s.mat', hemi{ih});
      %   cfg.surfinflated      = sprintf('surface_inflated_%s_caret.mat', hemi{ih});
      
      cfg.nslices = 10;
      if max(mri_standard.dim) == 72;        cfg.slicerange = [17 40];      else; cfg.slicerange = [100 250];      end
      cfg.funparameter  = 'functional';
      cfg.funcolorlim   = [-8 8];
      cfg.maskparameter = 'mask';
      cfg.maskstyle     = 'opacity';  % colormix
      
      cmap_cool = flipud(colormaps(3));
      cmap_hot = colormaps(2);
      range_nonwhite = cfg.funcolorlim(2) - BSthresh(2);
      range_white = BSthresh(2);
      cmap_middle = ones(floor(length(cmap_hot)/range_nonwhite * range_white) * 2, 3); %*2 to make symmetrical towards blue
      cfg.funcolormap = [cmap_cool; cmap_middle; cmap_hot];
      
      cfg.opacitymap    = 'auto';
      cfg.opacitylim    = 'auto';
      %     cfg.sphereradius = 100;
      cfg.camlight = 'no';
      % cfg.atlas = '/Users/kloosterman/Library/CloudStorage/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/atlas/spm_anatomy/AllAreas_v18.mat'
      tok = tokenize(files{ifile}, '/');
      tit = sprintf('%s_',  tok{[3,1,4]}); % tok{1:end-1} 
      cfg.title = sprintf('%s LV1 p=%1.3f', tit, result.perm_result.sprob(1) );
      ft_sourceplot(cfg, mri_BSratio, mri_standard); %  hack in ft_sourceplot for slice in 833: always 2 rows
      % lt = light;
      name = sprintf('%s_', tok{1:end});
      if strcmp(cfg.method, 'surface')
        view ([270 0])
        saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('surf_%s%s_270.pdf', cfg.title, hemi{ih})))
        saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('surf_%s%s_270.png', cfg.title, hemi{ih})))
        view ([90 0])
        saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('surf_%s%s_90.pdf', cfg.title, hemi{ih})))
        saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('surf_%s%s_90.png', cfg.title, hemi{ih})))
      elseif strcmp(cfg.method, 'slice')
        title(cfg.title);
        f=gcf;
        f.Position = [   680   703   600   300 ];
        saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('%s_slice.pdf', name)))
        saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('%s_slice.png', name)))
        %         saveas(gcf, fullfile(PREOUTplot, 'epsc2', sprintf('slice_%s.eps', cfg.title)), 'epsc2')
      end
      if strcmp(method, 'slice'); break; end
    end
  end
  
end

%% plot Brain scores in separate figures
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/';
demean_subj = 0;
close all
for ifile = 1:6 %length(files)
  disp( fullfile(PREIN, [files{ifile} 'result.mat']));
  load( fullfile(PREIN, [files{ifile} 'result.mat']));
  nsub = sum(result.num_subj_lst);
  ncond = length(cond_name);
  
  idx = 1:result.num_subj_lst(1)*result.num_conditions;
  usc={};
  usc{1} = reshape(result.usc(idx,1), [], ncond ) / 10000 .* flip_lv_sign(ifile);
  usc{2} = reshape(result.usc(idx(end)+1:end,1), [], ncond )/ 10000 .* flip_lv_sign(ifile);
  if demean_subj
    usc{1} = usc{1} - mean(usc{1},2);
    usc{2} = usc{2} - mean(usc{2},2);
  end
  
  agecolors = [1 0.5 0.5; 0.5 0.5 1]; agenames = {'Younger' 'Older'};
  clear p ylims
  f = figure; f.Position =  [ 680   506   2*150   250 ];
  tiledlayout(1,2,'TileSpacing','compact');
  for iage = 1:2
    %     iplot=iplot+1;
    %     subplot(1, 6*2,iplot); hold on;
    nexttile
    plotSpread_incmarkersz( usc{iage}, 'distributionColors', repmat(agecolors(iage,:), ncond, 1) );
    p(iage) = plot(mean(usc{iage}), 'Color', agecolors(iage,:), 'Linewidth', 3);
    xlabel('HMAX bin');    ylabel('Brainscore (a.u.*10000)');    % legend(p, agenames); legend boxoff
    ax = gca;
    ylims(iage,:) = ax.YLim;
    if iage==2
      tok = tokenize(fileparts(files{ifile}), '/');
      title(sprintf('%s ', tok{[3,1,4]}), 'HorizontalAlignment', 'right'); 
      ax.YLim = [min(ylims(:)) max(ylims(:))];
      ax.YAxis.Visible = 0;
%       subplot(1, 6*2,iplot-1); 
%       ax=gca; ax.YLim = [min(ylims(:)) max(ylims(:))];
    end
  end
  name = sprintf('%s_', tok{1:end});
  saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('%s_Brainscores.pdf', name)))
  saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('%s_Brainscores.png', name)))

end


%% plot Brain scores all in one figure
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/';
demean_subj = 0;
% close all
f = figure; f.Position =  [ 680   506   6*300   500 ]; iplot=0;
tiledlayout(1,12,'TileSpacing','compact');
for ifile = 1:6 %length(files)
  disp( fullfile(PREIN, [files{ifile} 'result.mat']));
  load( fullfile(PREIN, [files{ifile} 'result.mat']));
  nsub = sum(result.num_subj_lst);
  ncond = length(cond_name);
  
  idx = 1:result.num_subj_lst(1)*result.num_conditions;
  usc={};
  usc{1} = reshape(result.usc(idx,1), [], ncond ) / 10000 .* flip_lv_sign(ifile);
  usc{2} = reshape(result.usc(idx(end)+1:end,1), [], ncond )/ 10000 .* flip_lv_sign(ifile);
  
  if demean_subj
    usc{1} = usc{1} - mean(usc{1},2);
    usc{2} = usc{2} - mean(usc{2},2);
  end
  
  agecolors = [1 0.5 0.5; 0.5 0.5 1]; agenames = {'Younger' 'Older'};
  clear p ylims
  for iage = 1:2
    %     iplot=iplot+1;
    %     subplot(1, 6*2,iplot); hold on;
    nexttile
    plotSpread_incmarkersz( usc{iage}, 'distributionColors', repmat(agecolors(iage,:), ncond, 1) );
    p(iage) = plot(mean(usc{iage}), 'Color', agecolors(iage,:), 'Linewidth', 3);
    xlabel('HMAX bin');    ylabel('Brainscore (a.u.*10000)');    % legend(p, agenames); legend boxoff
    ax = gca;
    ylims(iage,:) = ax.YLim;
    if iage==2
      tok = tokenize(fileparts(files{ifile}), '/');
      title(sprintf('%s ', tok{[3,1,4]}), 'HorizontalAlignment', 'right'); 
      ax.YLim = [min(ylims(:)) max(ylims(:))];
      ax.YAxis.Visible = 0;
%       subplot(1, 6*2,iplot-1); 
%       ax=gca; ax.YLim = [min(ylims(:)) max(ylims(:))];
    end
  end
end
%   name = sprintf('%s_', tok{1:end-1});
%   saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('Brainscores_%s.pdf', name)))
%   saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('Brainscores_%s.png', name)))
  saveas(gcf, fullfile(PREOUTplot, 'pdf', sprintf('Brainscores.pdf')))
  saveas(gcf, fullfile(PREOUTplot, 'png', sprintf('Brainscores.png')))





% %% load hmax per bin values for each subject
% PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/taskPLS/gaze-specific/';
% cd(PREIN)
% load SDbold_vs_HMAX_youngold_44_41_BfMRIresult.mat
% 
% ages = {'young', 'old'}; hmax_meanperbin=[];  subjlist={};
% for iage = 1:2
%   cd(ages{iage})
%   dirlist = dir('sub-*.mat');
%   for i = 1:length(dirlist)
%     disp(dirlist(i).name);
%     tmp = load(dirlist(i).name);
%     hmax_meanperbin{iage}(i,:) = tmp.hmax_meanperbin;
%     tok = tokenize(dirlist(i).name, '_');
%     subjlist{iage}{i} = tok{1};
%   end
%   cd ..
% end
% %%
% % plot hmax per bin
% figure; subplot(2,2,1); hold on
% plot(mean(hmax_meanperbin{1}))
% plot(mean(hmax_meanperbin{2}))
% legend({'YA' 'OA'})
% subplot(2,2,2); hold on
% plot(std(hmax_meanperbin{1}))
% plot(std(hmax_meanperbin{2}))
% legend({'YA' 'OA'})
% subplot(2,2,3); hold on
% plot(hmax_meanperbin{1}')
% plot(hmax_meanperbin{2}')
% legend({'YA' 'OA'})
% 
% 
% %%
% disp 'rmcorr hmax vs brain scores'
% usc_YA_demean = usc_YA - mean(usc_YA,2);
% usc_OA_demean = usc_OA - mean(usc_OA,2);
% hmax_YA_demean = hmax_meanperbin{1} - mean(hmax_meanperbin{1},2);
% hmax_OA_demean = hmax_meanperbin{2} - mean(hmax_meanperbin{2},2);
% [r_YA, p_YA] = corr(usc_YA_demean(:), hmax_YA_demean(:), 'type', 'Spearman')
% [r_OA, p_OA] = corr(usc_OA_demean(:), hmax_OA_demean(:), 'type', 'Spearman')
% f=figure;
% subplot(1,2,1); scatter(usc_YA_demean(:), hmax_YA_demean(:)); axis square; box on; lsline
% subplot(1,2,2); scatter(usc_OA_demean(:), hmax_OA_demean(:)); axis square; box on; lsline
% 
% 
