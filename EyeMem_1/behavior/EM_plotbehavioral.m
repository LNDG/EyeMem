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
%     "su8b-58"
%     "sub-62"
%     "sub-65"
%     "sub-76"
%     "sub-84"
%     "sub-93"
%     "sub-95"

%% get data into structure with cells for age groups

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

% behavnames = {  'dprime'; 'propcorrect'; 'criterion';   'RT' ; 'RTsd' ; 'omissions' ;  };
behavnames = {  'dprime'; 'propcorrect'; 'criterion';   'RT' ; 'RTsd' ; 'omissions' ; 
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

%% plot SDT and DDM behavior
SAV=1;
close all
nrow=4; ncol=6;

f = figure; iplot=0;
f.Position =[   680   467   85*ncol   150*nrow];
Fontsize = 6;

for iphase = 1:2
  for im = 1:length(behavnames)
    iplot=iplot+1;
    subplot(nrow,ncol,iplot); hold on; % axis tight
    %     handle = plotSpread_incmarkersz( data{im,iphase}, 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );
    handle = plotSpread_incmarkersz( data.(phase_names{iphase}).(behavnames{im}), 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );
    
    br = bar(transpose(cellfun(@mean, data.(phase_names{iphase}).(behavnames{im}) )));
    br.FaceColor = 'flat';
    br.CData(1,:) = [255 0 0];
    br.CData(2,:) = [0 0 255];
    br.FaceAlpha = .4;
    br.EdgeAlpha = 0;
    shading flat
    ax=gca; ax.XTickLabel = age_groups; ax.XTickLabelRotation = 45;
    ax.FontSize = 8;
    
    [h,p] = ttest2(data.(phase_names{iphase}).(behavnames{im}){1}, data.(phase_names{iphase}).(behavnames{im}){2});
    if p < 0.05
      sigstar({[1,2]},p);
    end
    title(sprintf('%s\n%s\n', phase_names{iphase}, behavnames{im} ))
    if strcmp( behavnames{im}, 'z')
      ylim([0.4 0.65])
      plot(ax.XLim, [0.5 0.5], 'k', 'Linewidth', 0.5)
      br.BaseLine.BaseValue = 0.5;
    end
  end
  iplot=iplot+1;
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('EM_behavior.png')))
  cd(b.PREOUT)
end



%% correlate behavior with each other dprime vs RT

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

%% correlate behavior with each other study vs test
%TODO FIX, run this in spearman

f=figure;iplot=0;
f.Position =        [1000         498         961         840];
for iage = 1:2
  for im = 1:size(data,1)
    iplot=iplot+1;
    subplot(4,4,iplot)
    [r,p]=corr( data{im,1}{iage},  data{im,2}{iage} );
    scatter( data{im,1}{iage},  data{im,2}{iage}, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 50, 'LineWidth', 1)
    title(sprintf('%s\n%s r=%1.2f', behavnames{im}{:}, age_groups{iage}, r))
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

