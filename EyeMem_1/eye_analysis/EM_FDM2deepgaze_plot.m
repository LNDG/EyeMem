function EM_FDM2deepgaze_plot(maps)

addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/FACS')
addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/plotting')
% PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/plots';
PREOUT = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/plots';
SAV=1;
close all

%% plotting
disp 'time course of corr'
figure; hold on
agecolors = {'r' 'b'}; 
clear h
for iage= 1:2
  nsub = length(maps(iage).corr);
  h(iage) = shadedErrorBar(maps(iage).time, mean(maps(iage).corr), std(maps(iage).corr) / sqrt(nsub), agecolors{iage}, 1 );
end
permstats = 0
if permstats
  timelock=[];
  for iage = 1:2
    timelock(iage).time = maps(iage).time;
    timelock(iage).label = {'corrmap'};
    timelock(iage).fixdens = shiftdim(maps(iage).corr, -1);
    timelock(iage).dimord = 'chan_subj_time';
  end
  cfg = [];
  cfg.parameter = 'fixdens';
  cfg.method           = 'montecarlo';
  cfg.statistic        = 'indepsamplesT';  %depsamplesT ft_statfun_correlationT ft_statfun_correlationT_corrcol
  cfg.correctm         = 'cluster';  %'no'
  cfg.clusteralpha     = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.tail             = 0;
  cfg.clustertail      = 0;
  cfg.alpha            = 0.025;
  cfg.numrandomization = 10000;
  cfg.minnbchan        = 0;
  nsub(1) = size(timelock(1).fixdens,2);
  nsub(2) = size(timelock(2).fixdens,2);
  cfg.design = ones(2,sum(nsub));
  cfg.design(1,:) = [1:nsub(1) 1:nsub(2)];
  cfg.design(2,nsub(1)+1:end) = 2;
  cfg.ivar     = 2;
  
  stat = ft_timelockstatistics(cfg, timelock(1), timelock(2));
  if any(stat.mask)
    pl = plot_sig_bar(timelock(1).time, stat.mask, [], 6, [0 0 0]); % vpos, height, Color
  end
else
  [~,p]=ttest2(maps(1).corr,maps(2).corr);
end


legend([h.mainLine], {maps.agegroup}); legend boxoff
xlabel('Time from picture onset (s)')
ylabel( sprintf('%s correlation', maps(iage).corrtype))
title(sprintf('Corr subjectgaze vs. Deepgaze\n%dx%d AOIs, omit_centerAOI=%d',  maps(iage).nbins_x, maps(iage).nbins_y, maps(iage).omit_centerAOI))
if SAV
  outfile = sprintf('corrtc_%dx%d_AOIs_omit_centerAOI=%d.png', maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI);
  saveas(gcf, fullfile(PREOUT, outfile ))
end

%% time course of deepgaze correlation vs behavior correlation
f=figure; f.Position = [744 358 565 692];
for iage= 1:2
  subplot(2,1,iage)
  plot(maps(iage).time, maps(iage).corr2behav)
  legend({'Dprime all test trials' 'RT all test trials' 'Prop. remembered trials' 'RT remembered trials' 'RT all old trials'})
  legend boxoff
  xlabel('Time from picture onset (s)')
  ylabel(sprintf('%s correlation', maps(iage).corrtype))
  title(sprintf('%s (N=%d), correlation Deepgaze-to-subjectgaze match vs behavior', maps(iage).agegroup, size(maps(iage).behavior,1) ))
  ylim([-.4 .4])
end
if SAV
  outfile = sprintf('corrDG2SG_vs_behavior_%dx%d_AOIs_omit_centerAOI=%d.png', maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI);
  saveas(gcf, fullfile(PREOUT, outfile ))
end
cd(PREOUT)

disp 'TODO, update the below'
return
%%

disp 'plot scatters for each subject'
for iage = 1:2
  f=figure;
  f.Position = [ 2356         239        1243         818];
  for isub = 1:length(maps(iage).corr)
    subplot(6,8,isub);
%     scatter(maps(iage).fixation(isub,:)', maps(iage).deepgaze(isub,:)', '.');
%     heatscatter(maps(iage).fixation(isub,:)', maps(iage).deepgaze(isub,:)', PREOUT, 'test.png');
%     scatplot(maps(iage).fixation(isub,:)', maps(iage).deepgaze(isub,:)');
%     dscatter(maps(iage).fixation(isub,:)', maps(iage).deepgaze(isub,:)');
    dscatter(maps(iage).corrdat{isub}(:,1), maps(iage).corrdat{isub}(:,2));
    xlim([0 100])
    ylim([-0.05e-4 0.25e-4])
    title(sprintf('r=%1.2f', maps(iage).corr(isub,1)))
  end
  suptitstring = sprintf('%s, %dx%d AOIs, omit_centerAOI=%d %s',  maps(iage).agegroup, maps(iage).nbins_x, maps(iage).nbins_y, maps(iage).omit_centerAOI, maps(iage).corrtype);
  suptitle(suptitstring)  
  if SAV
    suptitstring(isspace(suptitstring)) = [];
    saveas(gcf, fullfile(PREOUT, [suptitstring '.png'] ))
%     saveas(gcf, fullfile(PREOUT, sprintf('%s_%s_corrpersub.png',  maps(iage).agegroup )))
  end
end

%%
disp 'plotspread bubble bar plot YA vs OA'
figure;
plotSpread_incmarkersz({maps.corr}, 'distributionColors', [1 0.5 0.5; 0.5 0.5 1]) % , 'categoryLabels', {maps.agegroup}
% mean lines
line([0.75 1.25]', [nanmean(maps(1).corr) nanmean(maps(1).corr)]',  'Color', 'r', 'Linewidth', 2)
line([1.75 2.25]', [nanmean(maps(2).corr) nanmean(maps(2).corr)]'',  'Color', 'b', 'Linewidth', 2)
% stats
[~,p] = ttest2(maps(1).corr, maps(2).corr);
ax=gca;
% text(1, ax.YLim(2), sprintf('p=%1.3f', p))
ax.XTickLabel = {maps.agegroup};
title(sprintf('%dx%d AOIs, omit_centerAOI=%d %s, p=%1.3f', maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI, maps(iage).corrtype, p))
xlim([0.5 2.5])
ylabel('Correlation per subj FDM vs Deepgaze')
if SAV
  outfile = sprintf('%dx%d_AOIs_omit_centerAOI=%d_spread.png', maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI)
  saveas(gcf, fullfile(PREOUT, outfile ))
end
cd(PREOUT)
% flipud?
% fix nans
