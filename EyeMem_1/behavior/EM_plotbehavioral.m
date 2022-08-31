function EM_plotbehavioral(b)

% exp_phases = {'study'; 'test'};
% depvars = {'dprime' 'criterion' 'RT' 'propcorrect'  }; % 'RTsd' 'RTsd2'
%
% these subjects waited until after picture presentation until responding
% during test
% b.participants.participant_id(b.t>1)
% 
% ans = 
% 
%   9×1 string array
% 
%     "sub-15"
%     "sub-32"
%     "sub-58"
%     "sub-62"
%     "sub-65"
%     "sub-76"
%     "sub-84"
%     "sub-93"
%     "sub-95"

if nargin ==0
  load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat')
  b=behavior;
end

%% get data into structure with cells for age groups

% find out which SUBJ are still in the mix in fMRI
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/1TRspertrial/ftsource/std_3bins/fixednbins/behavPLSvsDDM/v/gaze-specific/linearfit_fitcoeff1';;
cd(fullfile(PREIN))
SUBJ = [];
subjlist = dir('sub*_BfMRIsessiondata.mat');
for isub=1:length(subjlist)
  tmp = tokenize(subjlist(isub).name, '_');
  SUBJ(isub,1) = find(b.participants.participant_id == tmp{1});
end
SUBJ = sort(SUBJ);
SUBJage = b.participants.group(SUBJ,1);

% behavnames = {  'dprime'; 'propcorrect'; 'criterion';   'RT' ; 'RTsd' ; 'omissions' ;  };
behavnames = { 'propcorrect';   'RT' ; 'propNo';  'criterion'; 'dprime';   'RTsd' ; % 'omissions' ; 
  'v'; 'a'; 't';   'z' ; 'dc' ; };
age_groups = {'Young' 'Older'};
phase_names = {'study' 'test'};
age_colors = [1 0.5 0.5; 0.5 0.5 1];

data = [];
data.behavnames = behavnames;
data.age_groups = age_groups;
data.phase_names = phase_names;
data.age_colors = age_colors;
SUBJagegroups = {'young' 'old'};
for iphase = 1:2
  for im = 1:length(behavnames)
    for iage = 1:2
      curb = getfield(b, behavnames{im});
      if size(curb,2) == 1
        curb = curb(SUBJ,1); % select SUBJ made it to fMRI
      else
        curb = curb(SUBJ,iphase, end); % select SUBJ made it to fMRI
      end
      %       data{im,iphase}{iage} = curb(SUBJage == age_groups{iage});
      data.(phase_names{iphase}).(behavnames{im}){iage} = curb(SUBJage == SUBJagegroups{iage});
    end
  end
end

%% plot modelfree behavior

behavnames = { 'propcorrect';   'RT' ; 'propNo';  };

SAV=1;
close all
nrow=2; ncol=3;

f = figure; iplot=0;
f.Position =[   680   467   85*ncol   150*nrow];
Fontsize = 6;

for iphase = 1:2
  disp(iphase)
  for im = 1:length(behavnames)
    disp(behavnames{im})
    iplot=iplot+1;
    subplot(nrow,ncol,iplot); hold on; % axis tight
    %     handle = plotSpread_incmarkersz( data{im,iphase}, 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );
    handle = plotSpread_incmarkersz( data.(phase_names{iphase}).(behavnames{im}), 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );
    
    bardat = cellfun(@mean, data.(phase_names{iphase}).(behavnames{im}) )
    barSD = cellfun(@std, data.(phase_names{iphase}).(behavnames{im}) )
    br = bar(transpose(bardat));
    br.FaceColor = 'flat';
    br.CData(1,:) = [255 0 0];
    br.CData(2,:) = [0 0 255];
    br.FaceAlpha = .4;
    br.EdgeAlpha = 0;
    shading flat
    ax=gca; ax.XTickLabel = age_groups; ax.XTickLabelRotation = 45;
    ax.FontSize = 8;
    
    [~,p] = ttest2(data.(phase_names{iphase}).(behavnames{im}){1}, data.(phase_names{iphase}).(behavnames{im}){2})
    if p < 0.05
      sigstar({[1,2]},p);
    end
%     title(sprintf('%s\n%s\n', phase_names{iphase}, behavnames{im} ))
    if strcmp( behavnames{im}, 'propNo') 
      br.BaseLine.BaseValue = 0.5;
      ylabel('Proportion No')
      title(sprintf('\n\nBias\n'))
      ylim([0 1])
    end
    if strcmp( behavnames{im}, 'propcorrect')
      ylim([0.5 1])
      ylabel('Proportion correct')
      title(sprintf('\n\nAccuracy\n' ))
    end
    if strcmp( behavnames{im}, 'RT')
      ax.YTick = [1 1.5 2];
      ylim([0.5 2])
      ylabel('Reaction time (s)')
%       title(sprintf('%s\n', phase_names{iphase} ))
      title(sprintf('%s\nSpeed\n', phase_names{iphase} ))
    end
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_modelfree.pdf' )))
%   saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_modelfree.eps' )), 'epsc')
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_modelfree.png')))
  cd(b.PREOUT)
end

%% plot modelfree behavior, tasks grouped

behavnames = { 'propcorrect';   'RT' ; 'propNo';  };

SAV=1;
close all
nrow=1; ncol=3;

f = figure; iplot=0;
f.Position =[   680   467   120*ncol   150*nrow];
Fontsize = 6;

% for iphase = 1:2
  for im = 1:length(behavnames)
    iplot=iplot+1;
    subplot(nrow,ncol,iplot); hold on; % axis tight
    
    plotdat = [data.study.(behavnames{im});  data.test.(behavnames{im})];
    
    br = bar(transpose(cellfun(@mean, plotdat' )));
    br(1).FaceColor = 'flat';
    br(1).CData=[1 0 0];
    br(1).FaceAlpha = .4;
    br(1).EdgeAlpha = 0;
    br(2).FaceColor = 'flat';
    br(2).CData = [0 0 1];
    br(2).FaceAlpha = .4;
    br(2).EdgeAlpha = 0;
    shading flat
    ax=gca; ax.XTick = [1 2]; ax.XTickLabel = phase_names; ax.XTickLabelRotation = 45;
    ax.FontSize = 8;

    hold on
    Xlocs = [br(1).XEndPoints br(2).XEndPoints];
%     handle = plotSpread_incmarkersz(plotdat(:) , 'distributionColors', [1 0.5 0.5; 1 0.5 0.5;  0.5 0.5 1; 0.5 0.5 1], ...
%       'xValues', Xlocs, 'xMode', 'auto');
%     handle = plotSpread_incmarkersz(plotdat(:) , 'distributionColors', [1 0.25 0.25; 1 0.25 0.25;  0.25 0.25 1; 0.25 0.25 1], ...
%       'xValues', Xlocs, 'xMode', 'auto');
    handle = plotSpread_incmarkersz(plotdat(:) , 'distributionColors', [1 0 0; 1 0 0; 0 0 1; 0 0 1], ...
      'xValues', Xlocs, 'xMode', 'auto');
    
    if strcmp( behavnames{im}, 'propNo')
      br(1).BaseLine.BaseValue = 0.5;
      br(2).BaseLine.BaseValue = 0.5;
      ylabel('Prop. No responses')
      title(sprintf('\n\nBias\n'))
      ylim([0 1])
    end
    if strcmp( behavnames{im}, 'propcorrect')
      ax=gca; ax.YTick = [0.6 0.8 1];
      ylim([0.5 1])
      ylabel('Prop. correct')
      title(sprintf('\n\nAccuracy\n' ))
    end
    if strcmp( behavnames{im}, 'RT')
      ax=gca; ax.YTick = [1 1.5 2];
      ylim([0.5 2])
      ylabel('Reaction time (s)')
%       title(sprintf('%s\n', phase_names{iphase} ))
%       title(sprintf('%s\nSpeed\n', phase_names{iphase} ))
      title(sprintf('\n\nSpeed\n' ))
      lh = legend(handle{1}(1:2), age_groups)
    end
    p=[];
    [~,p(1)] = ttest2(data.study.(behavnames{im}){1}, data.study.(behavnames{im}){2});
    [~,p(2)] = ttest2(data.test.(behavnames{im}){1}, data.test.(behavnames{im}){2});
    p(p>0.05) = NaN;
    hs = sigstar({ Xlocs([1 3]) Xlocs([2 4]) },p);
    if strcmp( behavnames{im}, 'RT')
      lh = legend(br, age_groups); legend boxoff; lh.Position = [0.4472 0.6436 0.1722 0.1500];
    end
  end
% end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_modelfree_grouped.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_modelfree_grouped.png')))
  cd(b.PREOUT)
end

%% plot DDM behavior

behavnames = { 'v';   'a' ; 't';  'z' ; 'dc';  };
behavtits = {'Drift rate' 'Boundary\nseparation' 'Non-decision\ntime' 'Starting point\nbias' 'Drift bias'};
SAV=1;
close all
nrow=2; ncol=5;

f = figure; iplot=0;
f.Position =[   680   467   75*ncol   150*nrow];
Fontsize = 6;

for iphase = 2
  disp(phase_names{iphase})
  for im = 1:length(behavnames)
    disp(behavnames{im})
    iplot=iplot+1;
    subplot(nrow,ncol,iplot); hold on; % axis tight
    %     handle = plotSpread_incmarkersz( data{im,iphase}, 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );
    handle = plotSpread_incmarkersz( data.(phase_names{iphase}).(behavnames{im}), 'distributionColors', [1 0 0; 0 0 1] );
    
    bardat = cellfun(@mean, data.(phase_names{iphase}).(behavnames{im}) )
    br = bar(transpose(bardat));
    br.FaceColor = 'flat';
    br.CData(1,:) = [255 0 0];
    br.CData(2,:) = [0 0 255];
    br.FaceAlpha = .4;
    br.EdgeAlpha = 0;
    shading flat
    ax=gca; ax.XTickLabel = age_groups; ax.XTickLabelRotation = 45;
    ax.FontSize = 8;
    
    title(sprintf([behavtits{im} '\n']))
    if strcmp( behavnames{im}, 'z')
      ylim([0.4 0.65])
      br.BaseLine.BaseValue = 0.5;
    end
    if strcmp( behavnames{im}, 'v')
      ylabel('Parameter estimate')
    end
    if strcmp( behavnames{im}, 'z')
      ax.YTick = [0.4 0.5 0.6];
    end
    [~,p] = ttest2(data.(phase_names{iphase}).(behavnames{im}){1}, data.(phase_names{iphase}).(behavnames{im}){2})
    if p>0.05; p=NaN; end
    sigstar({[1,2]},p);
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_DDM.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_DDM.png')))
  cd(b.PREOUT)
end


%% correlate behavior with each other: dprime vs RT

corrpairs = {'dprime' 'RT';  'dprime' 'criterion';  'RT' 'criterion'; };
corrtype = 'Spearman';

f=figure;iplot=0;
f.Position =        [1000         498         500            350];
for iage = 1:2
  for im = 1:size(corrpairs,1)
    iplot=iplot+1;
    subplot(2,3,iplot)
    corrdat = [data.(phase_names{iphase}).(corrpairs{im,1}){iage} data.(phase_names{iphase}).(corrpairs{im,2}){iage}];
    [r,p]=corr( corrdat(:,1), corrdat(:,2), 'type', corrtype );
    scatter( corrdat(:,1), corrdat(:,2), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 50, 'LineWidth', 1)
    title(sprintf('%s %s\nr=%1.2f, p=%1.3f', age_groups{iage}, corrtype, r,p))
    if p<0.05
      lsline
    end
    xlabel(corrpairs(im,1))
    ylabel(corrpairs(im,2))
    box on; axis square
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_corrSDT.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_corrSDT.png')))
  cd(b.PREOUT)
end

%% correlate behavior with each other study vs test
%TODO FIX, run this in spearman
SAV=1;
corrmeas = {'propcorrect', 'RT', 'propNo'};

f=figure;iplot=0;
f.Position =        [1000         498         961         840];
for iage = 1:2
  for im = 1:length(corrmeas)
    iplot=iplot+1;
    subplot(4,3,iplot)
    
    corrdat = [data.(phase_names{1}).(corrmeas{im}){iage} data.(phase_names{iphase}).(corrmeas{im}){iage}];
    
    [r,p]=corr( corrdat(:,1), corrdat(:,2), 'type', corrtype );
    scatter( corrdat(:,1), corrdat(:,2), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 50, 'LineWidth', 1)
    title(sprintf('%s %s %s\nr=%1.2f, p=%1.3f', corrmeas{im}, age_groups{iage}, corrtype, r,p))
    if p<0.05
      lsline
    end
    ylabel('test')
    xlabel('study')
    box on; axis square
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_corr_studyvstest.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior_corr_studyvstest.png')))
  cd(b.PREOUT)
end


%% correlate behavior with each other: DDM params

corrpairs = {'v' 't';  'v' 'a';  'v' 'dc';  'v' 't'; ...
  'a' 't';  'dc' 'a';  'v' 'dc';  'v' 't';};
corrtype = 'Spearman';
% corrtype = 'Pearson';
close all
% corrdat={};
corrTable = {};
for iage = 1:2
%   corrdat{iage} = [data.test.v{iage} data.test.a{iage} data.test.t{iage} data.test.z{iage} data.test.dc{iage}];  
  corrTable{iage} = table(data.test.v{iage}, data.test.a{iage}, data.test.t{iage}, data.test.z{iage}, data.test.dc{iage}, 'VariableNames', {'v','a','t','z', 'dc'});
end
corrTable{3} = [corrTable{1}; corrTable{2}];
age_groups{3} = 'All ages';

f=figure;iplot=0;
f.Position =        [779         870        1275         350];
for iage = 1:3
  subplot(1,3,iage)
  corrplot( corrTable{iage}, 'type', corrtype )
  title(sprintf('%s, %s correlation matrix', age_groups{iage}, corrtype))
end

if SAV
  orient landscape
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_corrDDMbehav_%s.pdf', corrtype )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_corrDDMbehav_%s.png', corrtype )))
  cd(b.PREOUT)
end

%% plot RT distributions

edges = [0:0.1:2.5];

histdat = {};
for iage = 1:2
  for iphase = 1:2
    curb = b.RTsingletrial{iphase}(SUBJ,:);  % select SUBJ made it to fMRI
    for iage = 1:2
      data = curb(SUBJage == SUBJagegroups{iage},:);
      for isub = 1:size(data,1)
        histdat{iage}(isub,:,iphase) = histcounts(data(isub,:), edges);
      end
    end
  end
end

%
age_colors = {'r' 'b'};
phaseleg = {'study' 'test'};
f=figure;iplot=0;
f.Position =        [1000         498         400         150];
for iphase = 1:2
  iplot=iplot+1;
  subplot(1,2,iplot); hold on
  histmean=[];
  histsem=[];
  clear s
  for iage = 1:2
    histmean(iage,:) = mean(histdat{iage}(:,:,iphase));
%     histsem(iage,:) = std(histdat{iage}(:,:,iphase))/sqrt(size(histdat{iage},1));
    histsem(iage,:) = std(histdat{iage}(:,:,iphase));
    s(iage) = shadedErrorBar(edges(2:end)',histmean(iage,:)',histsem(iage,:)', {age_colors{iage}, 'Linewidth', 1.5}, 1);
  end
  ax=gca;
  ax.YLim = [0 ax.YLim(2)];
  ax.LineWidth = 1;
%   pl = plot(edges(2:end), histmean);
%   pl(1).Color = 'r';
%   pl(2).Color = 'b';
  legend([s.mainLine], age_groups)
  legend boxoff
  title(phaseleg{iphase})
  xlabel('Time (s)')
  ylabel('Frequency of occurrence')
end
  if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_RTdistributions.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_RTdistributions.png')))
  cd(b.PREOUT)
end

