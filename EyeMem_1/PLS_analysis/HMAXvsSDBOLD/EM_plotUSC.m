function EM_plotUSC(result)

if nargin == 0
  load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_3bins/fixednbins/gaze-specific/SDbold_vs_HMAX_youngold_44_41_BfMRIresult.mat', ...
    'result')
end

PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/plots';

nsub = result.num_subj_lst;
ncond = result.num_conditions;
% inds = [1 nsub(1);
%   nsub(1)+1  nsub(1)*2;
%   nsub(1)+1  nsub(1)*2];
% 
% ind_start = linspace(1, nsub(1)*ncond, ncond+1);
% ind_end   = linspace(nsub(1), nsub(1)*ncond, ncond);

lv=1;

% age groups in different cells
uscdat{1} = result.usc(1:nsub(1)*ncond,lv);
uscdat{2} = result.usc(nsub(1)*ncond+1:end,lv);

plotdat = {};
% ageinds = 0;
for iage =1:2
  inds=0;
  for icond = 1:ncond
    curinds = inds+1:(inds+nsub(iage));
    plotdat{icond, iage} = uscdat{iage}(curinds,:);
    
    inds = curinds(end);
    
  end
end
flipit = 1;
if flipit
  plotdat = cellfun(@(x) x*-1, plotdat, 'UniformOutput', false);
end

%% 3 bar plots, separate per age group
% ylims = [14.2e4 15.4e4; 11.6e4 12.8e4];
ylims = [1e5 2e5; 1e5 2e5];
titleg = {'YA' 'OA'};
f = figure;
f.Position = [ 1000         918         400         135 ];
for iage=1:2
  subplot(1,2,iage); hold on
  %   bar(cellfun(@mean, plotdat(:,iage)))
  plotSpread(plotdat(:,iage))
  xdat=nan(2,3);
  ydat=nan(2,3);
  for icond=1:3
    xdat(:,icond) = [icond-0.25 icond+0.25];
    ydat(:,icond) = [mean(plotdat{icond,iage}); mean(plotdat{icond,iage})];
  end
%   if flipit
%     ydat = ydat*-1;
%   end
  line(xdat, ydat, 'linewidth', 4)
  
  % plot single subj lines
  plotlines=0;
  if plotlines
    ydat = [plotdat{:,iage}];
    xdat = repmat([1 2 3], size(ydat,1),1);
    line(xdat', ydat', 'linewidth', 1, 'Color', [0.5 0.5 0.5])
  end
  
  ylim(ylims(iage,:))
  title(titleg{iage})
  ylabel('Brain score')
end
saveas(gcf, 'brainscores.pdf')
saveas(gcf, 'brainscores.png')

%% YA vs OA, separate per 3 HMAX bins

ylims = [ 1.15e5 1.6e5];
titleg = {'Young' 'Older'};
f = figure;
f.Position = [ 1000         918         160         175 ];
%   subplot(1,2,iage); hold on
barmeans = cellfun(@mean, plotdat);
barsem = cellfun(@(x) std(x)/sqrt(length(x)), plotdat);
b = barweb(barmeans, barsem);
ylim(ylims)
% b = gca
b.bars(1).FaceColor = 'flat';
b.bars(2).FaceColor = 'flat';
b.bars(1).CData = [255 0 0; 255 0 0; 255 0 0];
b.bars(2).CData = [0 0 255; 0 0 255; 0 0 255];
b.bars(1).FaceAlpha = .4;
b.bars(2).FaceAlpha = .4;
b.bars(1).EdgeAlpha = 0;
b.bars(2).EdgeAlpha = 0;
ax = gca;
ax.XTickLabel = {'Low' 'Medium' 'High'};
legend(titleg, 'Location', 'northoutside'); legend boxoff
xlabel('Feature richness')
ylabel('Brain score')
ax.LineWidth = 0.5

  
%   plotSpread(plotdat')
% %   barweb(plotdat')
%   xdat=nan(2,3);
%   ydat=nan(2,3);
%   for icond=1:3
%     xdat(:,icond) = [icond-0.25 icond+0.25];
%     ydat(:,icond) = [mean(plotdat{icond,iage}); mean(plotdat{icond,iage})];
%   end
% %   if flipit
% %     ydat = ydat*-1;
% %   end
%   line(xdat, ydat, 'linewidth', 4)
%   
%   % plot single subj lines
%   plotlines=0;
%   if plotlines
%     ydat = [plotdat{:,iage}];
%     xdat = repmat([1 2 3], size(ydat,1),1);
%     line(xdat', ydat', 'linewidth', 1, 'Color', [0.5 0.5 0.5])
%   end
%   
%   
%   title(titleg{iage})
%   ylabel('Brain score')
% end
saveas(gcf, fullfile(PREOUT, 'brainscores2.pdf'))
saveas(gcf, fullfile(PREOUT, 'brainscores2.png'))
saveas(gcf, fullfile(PREOUT, 'brainscores2.eps'), 'epsc' )

