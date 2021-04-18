% function EM_pls_SDbold_vs_HMAX_plot(source)

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
  addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/plotting/plotSpread')
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; %yesno or 2afc
end

PREIN = fullfile(basepath, 'variability', 'ftsource', 'SDbold_vs_HMAX');
cd(PREIN)

whichgroup = 'young'
% whichgroup = 'old'
PREOUT = fullfile(PREIN, whichgroup, '0_plots');
mkdir(PREOUT)
switch whichgroup
  case 'both'
    load source_HMAXbinned.mat
  otherwise % YA or OA
    cd(fullfile(PREIN, whichgroup))
    subjlist = dir('pls_source_*.mat');
    allsource = {};
    for isub = 1:length(subjlist)
      disp(subjlist(isub).name)
      load(subjlist(isub).name)
      allsource{end+1} = source;
    end
    source = cell2mat(allsource);
    clear allsource
end

disp 'Match the inside bool aka common coordinates'
inside_allsubj = [source.inside];
common_voxels = all(inside_allsubj,2);
[source.inside] = deal(common_voxels);  % cast this to all subj source, note that pow field is not updated

disp 'append subjects'
cfg=[];
cfg.parameter = 'pow';
cfg.keepindividual = 'yes';
source = num2cell(source);
source = ft_sourcegrandaverage(cfg, source{:});  % num2cell(source)

if ismac
  tmp = source;
  tmp.powdimord = 'pos';
  %       vol=300;
  %   tmp.anatomy = tmp.anatomy(:,:,:,vol);
  tmp.pow = nanmean(tmp.pow,3);
  tmp.pow = transpose(nanmean(tmp.pow,1));
  cfg=[];
  cfg.method = 'slice'; % slice ortho glassbrain vertex
  cfg.funparameter = 'pow';
  cfg.funcolorlim = 'zeromax';% [-300 300];
  ft_sourceplot(cfg, tmp)
%   title(sprintf)
end

%%
disp 'run PLS analysis'

plsdat = source.pow(:,source.inside,:);
plsdat = permute(plsdat, [1 3 2] ); % dimord: subj_cond_pos
s = size(plsdat);
datamat_lst = {reshape(plsdat, prod(s(1:2)), s(3))}; % "subject in condition"

num_subj_lst = size(source.pow,1);
num_cond = length(source.freq);
option = [];
option.method = 1; % [1] | 2 | 3 | 4 | 5 | 6 TODO make strings
option.num_boot = 1000; %500 % ( single non-negative integer )
option.num_perm = 1000; %( single non-negative integer )

result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);
save result result

%%
disp 'plot SDbold data'
close all
plsdat = source.pow(:,source.inside,:);
for icond = 1:num_cond
  f = figure;
  f.Position = [         373         296        1905        1049];
  for isub = 1:size(plsdat,1)
    subplot(6,8,isub);
    dat = squeeze(plsdat(isub,:,icond));
    %     plot(dat)
    s = size(dat);
    dat2 = nan(ceil(sqrt(s(2))));
    dat2(1:s(2)) = dat;
    imagesc(dat2, [0 1000]);
    axis square
    set(gca,'visible','off')
  end
end
%TODO plot hist

%%
disp 'prepare for plotting'
SAV = 1;

LVno = 2;
BSRthresh = [-3 3];

disp 'reshape '
stat= [];
stat.stat = zeros(size(source.inside));
stat.stat(source.inside) = result.boot_result.compare_u(:,LVno);
stat.inside = source.inside;
stat.pos = source.pos;
stat.dim = source.dim;

disp 'add BSR mask'
stat.mask = zeros(size(stat.inside));
stat.mask(source.inside) =  result.boot_result.compare_u(:,LVno) < BSRthresh(1) | result.boot_result.compare_u(:,LVno) > BSRthresh(2) ;
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

%
disp 'plotting'
% close all
cfg=[];
cfg.method = 'slice'; % slice ortho glassbrain vertex surface
cfg.funparameter = 'stat';
cfg.maskparameter = 'mask';
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, interp)

title(sprintf('%s, LV%d, p = %g', whichgroup, LVno, result.perm_result.sprob(LVno)))
hold on
% view ([90 0]) 
% figure %subplot(1,2,2)
% ft_sourceplot(cfg, vol2plot )
% view ([-90 0]) 
if SAV
  saveas(gcf, fullfile(PREOUT, sprintf('%s_%s_LV%d_%dbins.pdf', whichgroup, cfg.method, LVno, num_cond )))
  saveas(gcf, fullfile(PREOUT, sprintf('%s_%s_LV%d_%dbins.png', whichgroup, cfg.method, LVno, num_cond )))
  cd(PREOUT)
end

% plot bar plot LV
brainscores = reshape(result.usc(:,LVno), num_subj_lst, num_cond);
behavscores = reshape(result.vsc(:,LVno), num_subj_lst, num_cond);

f = figure;
f.Position = [       1000         918        1035         420 ];
subplot(1,3,1);
handle = plotSpread( num2cell(brainscores, 1), 'distributionMarkers', 'o' );% b r circles: , 'distributionColors', [1 0.5 0.5; 0.5 0.5 1]
for icond = 1:num_cond
  line([icond-0.25 icond+0.25]', [mean(brainscores(:,icond)) mean(brainscores(:,icond))]',  'Color', 'k', 'Linewidth', 4)
end

title(sprintf('%s LV%d brainscores', whichgroup, LVno))
ylabel('Brain score')
subplot(1,3,2);
hold on; box on; axis square
% bar(mean(behavscores))
[r,p] = corr(source.freq', mean(brainscores)')
scatter(source.freq',  mean(brainscores)', 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 200);
% scatter(source.freq', mean(brainscores)');
title(sprintf('r = %1.2f, p = %g', r,p))
xlabel('HMAX value')
ylabel('Brain score')

% plot explained variance for LVs
subplot(1,3,3); hold on; 
title('Percentage of Explained variance')
lvvar = (result.s/sum(result.s))*100;
plot(lvvar); 
scatter(1:num_cond, lvvar,'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 200)
if SAV
  saveas(gcf, fullfile(PREOUT, sprintf('%s_brainscores_LV%d_%dbins.pdf', whichgroup, LVno, num_cond )))
  saveas(gcf, fullfile(PREOUT, sprintf('%s_brainscores_LV%d_%dbins.png', whichgroup, LVno, num_cond )))
  cd(PREOUT)
end
