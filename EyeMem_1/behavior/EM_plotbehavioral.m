function EM_plotbehavioral(b)

% exp_phases = {'study'; 'test'};
% depvars = {'dprime' 'criterion' 'RT' 'propcorrect'  }; % 'RTsd' 'RTsd2'
%
%% plot modelfree

% find out which SUBJ are still in the mix in fMRI
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/criterion/linearfit_fitcoeff1';;
cd(fullfile(PREIN))
SUBJ = [];
subjlist = dir('sub*_BfMRIsessiondata.mat');
for isub=1:length(subjlist)
  tmp = tokenize(subjlist(isub).name, '_');
  SUBJ(isub,1) = find(b.participants.participant_id == tmp{1});
end
SUBJ = sort(SUBJ);
SUBJage = b.participants.group(SUBJ,1);

SAV=1;
behavnames = {
  {'dprime'}; {'propcorrect'}; {'criterion'};  { 'RT'} ;{ 'RTsd' };{ 'omissions' };
  };
age_groups = {'young' 'old'};
age_colors = [1 0.5 0.5; 0.5 0.5 1];

close all
nrow=8; ncol=5;

f = figure; iplot=0;
f.Position =[   680   467   85*ncol   100*nrow];
Fontsize = 6;
for im = 1:length(behavnames)
  for iphase = 1:2
    for iage = 1:2
      if isempty(behavnames{im})
        continue;
      end
      curb = getfield(b, behavnames{im}{:});
      curb = curb(SUBJ,iphase, 6); % select SUBJ made it to fMRI
      
      
%       ageinds = b.participants.group == b.age_groups{iage};
      
      
      data{iage} = curb(SUBJage == age_groups{iage});

    end
    iplot=iplot+1;
    subplot(nrow,ncol,iplot); hold on; % axis tight
    handle = plotSpread_incmarkersz( data, 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );

    br = bar(transpose(cellfun(@mean, data)));
    br.FaceColor = 'flat';
    br.CData(1,:) = [255 0 0];
    br.CData(2,:) = [0 0 255];
    br.FaceAlpha = .4;
    br.EdgeAlpha = 0;
    shading flat
    
    [h,p] = ttest2(data{1}, data{2});
    sigstar({[1,2]},p);

    
%     barmeans = cellfun(@mean, data)';
%     barsem = cellfun(@(x) std(x)/sqrt(length(x)), data)';
%     b = barweb(barmeans, barsem, 0.5, {'Young' 'Older'});
% %     ylim(ylims)
%     % b = gca
%     b.bars(1).FaceColor = 'flat';
%     b.bars(2).FaceColor = 'flat';
%     b.bars(1).CData = [255 0 0; 255 0 0];
%     b.bars(2).CData = [0 0 255; 0 0 255];
%     b.bars(1).FaceAlpha = .4;
%     b.bars(2).FaceAlpha = .4;
%     b.bars(1).EdgeAlpha = 0;
%     b.bars(2).EdgeAlpha = 0;
%     
    title(sprintf('%s %s\nruns %s', behavnames{im}{1} ))
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior.png')))
  cd(b.PREOUT)
end

return

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
% line([ones(1,size(data,1)); ones(1,size(data,1))+1], [data(:,1) data(:,2)]', 'Color',[0.66 0.66 0.66 ], 'Linewidth', 0.5)

handle = plotSpread( {data}, 'distributionMarkers', 'o', 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );% b r circles
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
