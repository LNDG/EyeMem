function [stat] = EM_FDM2deepgaze_plot(maps)

addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/FACS')
addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/plotting')
% PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/plots';
PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/plots';
SAV=1;
close all

%% plotting
disp 'time course of corr'
f = figure; 
f.Position = [ 1000         918         180         135 ];

hold on
agecolors = {'r' 'b'}; 
clear h
for iage= 1:2
  nsub = length(maps(iage).corr);
  h(iage) = shadedErrorBar(maps(iage).time, mean(maps(iage).corr), std(maps(iage).corr) / sqrt(nsub), agecolors{iage}, 0 );
end
ax=gca;
ax.XLim = [0 5];
ax.XTick = 0:5;
ax.XTickLabel = 0:5;
ax.FontSize = 8;

permstats = 1
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
%   cfg.alpha            = 0.99;
  cfg.numrandomization = 10000;
  cfg.minnbchan        = 0;
  nsub(1) = size(timelock(1).fixdens,2);
  nsub(2) = size(timelock(2).fixdens,2);
  cfg.design = ones(2,sum(nsub));
  cfg.design(1,:) = [1:nsub(1) 1:nsub(2)];
  cfg.design(2,nsub(1)+1:end) = 2;
  cfg.ivar     = 2;
  cfg.spmversion = 'spm12';
  
  stat = ft_timelockstatistics(cfg, timelock(1), timelock(2));
  if any(stat.mask)
    pl = plot_sig_bar(timelock(1).time, stat.mask, [], 3, [0 0 0]); % vpos, height, Color
  end
else
  [~,p]=ttest2(maps(1).corr,maps(2).corr);
end

leg = {'younger' 'older', 'p < 0.05, corrected'};
% legend([h.mainLine], {maps.agegroup}); legend boxoff
legend([h.mainLine], leg); legend boxoff
xlabel('Time from picture onset (s)')
ylabel( sprintf('Map similarity (%s’s r)', maps(iage).corrtype))

% title(sprintf('Corr subjectgaze vs. %s\n%dx%d AOIs, omit_centerAOI=%d', maps(iage).saliencymodel, maps(iage).nbins_x, maps(iage).nbins_y, maps(iage).omit_centerAOI))
if SAV
  outfile = sprintf('corrtc_%s_%dx%d_AOIs_omit_centerAOI=%d.png', maps(iage).saliencymodel, maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI);
  saveas(gcf, fullfile(PREOUT, outfile ))
  outfile = sprintf('corrtc_%s_%dx%d_AOIs_omit_centerAOI=%d.pdf', maps(iage).saliencymodel, maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI);
  saveas(gcf, fullfile(PREOUT, outfile ))
end

%% the below shows nothing interesting
return

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
  outfile = sprintf('corr_%s2SG_vs_behavior_%dx%d_AOIs_omit_centerAOI=%d.png',  maps(iage).saliencymodel, maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI);
  saveas(gcf, fullfile(PREOUT, outfile ))
end
cd(PREOUT)

ibehav=2; % 2 is RT study
itime=2;
iage = 2;
figure;
scatter(maps(iage).corrdatbehav{ibehav,itime}(:,1), maps(iage).corrdatbehav{ibehav,itime}(:,2))
[r,p]=corr(maps(iage).corrdatbehav{ibehav,itime}(:,1), maps(iage).corrdatbehav{ibehav,itime}(:,2), 'type', 'Spearman');
title(sprintf('r=%g p=%g', r, p))


alldat=[maps(iage).corrdatbehav{ibehav,itime}; maps(iage).corrdatbehav{ibehav,itime}];
corr(alldat(:,1), alldat(:,2), 'type', 'Spearman')


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
  outfile = sprintf('%dx%d_AOIs_omit_centerAOI=%d_spread.pdf', maps(1).nbins_x, maps(1).nbins_y, maps(1).omit_centerAOI)
  saveas(gcf, fullfile(PREOUT, outfile ))
end
cd(PREOUT)
% flipud?
% fix nans

