function EM_plotbehavioral(b)

% age_groups = {'young' 'old'};
% exp_phases = {'study'; 'test'};
% depvars = {'dprime' 'criterion' 'RT' 'propcorrect'  }; % 'RTsd' 'RTsd2'
%
%% plot modelfree
SAV=1;
behavnames = {
  {'dprime'};  {'criterion'};  { 'RT'} ;{ 'RTsd' };
  }; % all matrices
close all
nrow=8; ncol=5;

f = figure; iplot=0;
f.Position =[   680   467   85*ncol   100*nrow];
Fontsize = 6;
for im = 1:length(behavnames)
  for iphase = 2
    for iage = 1:2
      iplot=iplot+1;
      if isempty(behavnames{im})
        continue;
      end
      curb = getfield(b, behavnames{im}{:});
      ageinds = b.participants.group == b.age_groups{iage};
      data = squeeze(curb(ageinds,iphase, 6));  % dimord: 'subj_phase_cond'
      %   disp 'Drop NK1 high drift'
      %   data = data(2:end,:);
      subplot(nrow,ncol,iplot); hold on; % axis tight
      plotspreadfig(data, Fontsize, [], b.SUBJ);
      %       title(sprintf('%s %s\nruns %s', behavnames{im}{1}, diffleg{idiff}, avgtypestr{ia}), 'Fontsize', Fontsize-1)
    end
  end
end
if SAV
  %   saveas(gcf, fullfile(b.PREOUT, sprintf('behavior.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavior_%s_runs%s.png', diffleg{idiff}, avgtypestr{ia} )))
  cd(b.PREOUT)
end



minmax_vals = [0 3.5; -0.1 0.25; 0.5 1.4; 0.5 1]; % 0 0.4; 0 0.35;

% depvars = {'driftrate' 'boundsep' 'nondectime' };
% minmax_vals = [0 0.25; 0 0.2; -1 0.4; 0 0.4; 0 0.35; 0 1];

% depvars = {'RT_hits' 'RT_misses' 'RT_fas' 'RT_crs' };
% minmax_vals = [0 1.8; 0 1.8; 0 1.8; 0 1.8; 0 1.35; 0 1];

PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/plots/EyeMem/behavior';

close all
figure; set(gcf, 'Position', [0 -200 210*4 375*3])
iplot=0;

for idep = 1:length(depvars)
  for iage = 1:2
    %         subj_inds = behavior.agegroup == iage;
    subj_inds = behavior.participants.group == age_groups{iage}
    for iphase = 1:2 % study, test
      
      %             subj_dat{iphase, iage} = behavior.(exp_phases{iphase}).(depvars{idep})(subj_inds,6);
      subj_dat{iphase, iage} = cellfun(@mean, {behavior.(exp_phases{iphase})(subj_inds).(depvars{idep})})
      dat(iphase, iage) = nanmean(subj_dat{iphase, iage});
      %             dat_sem(iphase, iage) = nanstd(behavior.(exp_phases{iphase}).(depvars{idep})(subj_inds,6)) / sqrt(length(find(subj_inds)));
      dat_sem(iphase, iage) = nanstd(subj_dat{iphase, iage}) / sqrt(length(find(subj_inds)));
      length(find(subj_inds))
      %                     dat_sem(iphase, iage) = nanstd(behavior.(exp_phases{iphase}).(depvars{idep})(subj_inds,6));
      
    end
  end
  
  colors = cbrewer('qual', 'Set1', 10);
  % if we want each bar to have a different color, loop
  for iphase = 1:2
    iplot=iplot+1;
    subplot(5,4,iplot);  hold on;
    for iage = 1:2
      bar(iage, dat(iphase,iage), 'FaceColor',  colors(iage, :), 'EdgeColor', 'none', 'BarWidth', 0.6);
    end
    % show standard deviation on top
    h = ploterrbar(1:2, dat(iphase,:), [], dat_sem(iphase,:), 'k.', 'abshhxy', 0);
    set(h(1), 'marker', 'none'); % remove marker
    
    % label what we're seeing
    % if labels are too long to fit, use the xticklabelrotation with about -30
    % to rotate them so they're readable
    set(gca, 'xtick', [1 2], 'xticklabel', age_groups, ...
      'xlim', [0.5 2.5]);
    %         ylabel('Value'); % xlabel('Age group');
    ylim(minmax_vals(idep,:));
    
    dat1 = subj_dat{iphase,1};
    dat1 = dat1(~isnan(dat1));
    dat2 = subj_dat{iphase,2};
    dat2 = dat2(~isnan(dat2));
    datall = [mean(dat1) mean(dat2)];
    
    %         pval = randtest1d(dat1, dat2, 0, 1000);
    [~, pval] = ttest2(dat1, dat2);
    
    % if mysigstar gets 2 xpos inputs, it will draw a line between them and the sigstars on top
    mysigstar(gca, [1 2], max(datall) * 1.1, pval);
    
    title(sprintf('%s, %s N=%d, N=%d', depvars{idep}, exp_phases{iphase}, length(dat1), length(dat2) ));
    
    %         if strcmp(depvars(idep), 'dprime')
    %         correlate depvar for study and test
    if iphase == 2
      for iage = 1:2
        iplot=iplot+1;
        subplot(5,4,iplot);  hold on;
        scatter( subj_dat{1,iage}, subj_dat{2,iage} );
        
        [r, p] = corr(subj_dat{1,iage}, subj_dat{2,iage} );
        title(sprintf('%s r = %g, p = %g', age_groups{iage}, r, p ))
        if p < 0.05
          lsline
        end
        xlabel('study')
        ylabel('test')
      end
    end
    
  end
end
export_fig(fullfile(PREOUT, sprintf('behavior_%s.pdf', depvars{1})))
cd(PREOUT)
% barwebNK(dat, dat_sem, 0.75, exp_phases)%
% hold on; axis square; % after barweb call!
% legend(age_groups)

end

%% plot bars + lines
function ax = plotspreadfig(data, Fontsize, condlabels, SUBJ)
% lines between pts
line([ones(1,size(data,1)); ones(1,size(data,1))+1], [data(:,1) data(:,2)]', 'Color',[0.66 0.66 0.66 ], 'Linewidth', 0.5)

handle = plotSpread( {data(:,1) data(:,2)}, 'distributionMarkers', 'o', 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );% b r circles
set(handle{1}(1:2), 'MarkerSize', 3, 'Linewidth', 0.5)
% text(data(:,1)', data(:,2)', num2cell(SUBJ), 'Fontsize', Fontsize)

ax=gca;
ax.FontSize = Fontsize;
ax.XTickLabel = condlabels;

% mean lines
line([0.75 1.25]', [mean(data(:,1)) mean(data(:,1))]',  'Color', 'r', 'Linewidth', 2)
line([1.75 2.25]', [nanmean(data(:,2)) nanmean(data(:,2))]',  'Color', 'b', 'Linewidth', 2)

% stats
[~,p] = ttest(data(:,1), data(:,2));
text(1, ax.YLim(2), sprintf('p=%1.3f', p), 'Fontsize', Fontsize)
%   if p < 0.05 %   sigstar
xlim([0.5 2.5])
end
