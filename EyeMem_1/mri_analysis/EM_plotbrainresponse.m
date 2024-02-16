function source = EM_plotbrainresponse()
 
PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource';
cd(PREIN)
 
files = dir('source_sub*.mat');
load participantinfo.mat % TODO make this reliable
% agefolder = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
 
% load(files(1).name, 'pow')
load common_coords.mat
source_cell_mean =  {{} {}};   source_cell_sdbetween = {{} {}};   source_cell_sdwithin =  {{} {}};   source_cell_sdall =  {{} {}};   
indYA = 0; indOA = 0;
for ifile = 1:length(files) % TODO split YA and OA
  tok = tokenize(files(ifile).name, '_');
  [~, subj] = fileparts(tok{2});
  ind = Participants.participant_id == subj;
  agegroup = Participants(ind,:).group;
  if agegroup == 'young'
    iage = 1; 
  else
    iage = 2;
  end
  
  disp(files(ifile).name);
  source = load(files(ifile).name);
  source.inside = common_coords;
  source.time = 1:5;
  source.powdimord = 'pos';

  source_SDbetweentrials = source; % "evoked": trl to trl variability in evoked responses
  source_SDwithintrials = source; % "induced": variability within trials
  source_SDbetweentrials.pow = mean(source_SDbetweentrials.pow,3); % average within each trial to isolate across trl var
  source_SDbetweentrials.pow = repmat(source_SDbetweentrials.pow, [1 1 5]); % tile it, not sure if necessary
  source_SDwithintrials.pow = source_SDwithintrials.pow - source_SDbetweentrials.pow; % remove across trial variability to isolate within

  source_SDbetweentrials.pow = std(source_SDbetweentrials.pow(:,:),0,2); % now take SDs
  source_SDwithintrials.pow = std(source_SDwithintrials.pow(:,:),0,2);

  source_SDall = source;
  source_SDall.pow = std(source_SDall.pow(:,:),0,2); % SD across everything: 
  
  source_Meanall = source;
  source_Meanall.pow = mean(source_Meanall.pow(:,:),2); % mean across everything: 
  
  source_cell_mean{iage}{end+1} = source_Meanall;
  source_cell_sdbetween{iage}{end+1} = source_SDbetweentrials;
  source_cell_sdwithin{iage}{end+1} = source_SDwithintrials;
  source_cell_sdall{iage}{end+1} = source_SDall;
    
  clear source
end
 
for iage = 1:2
  %   source{1} = ft_sourcegrandaverage(cfg, source_cell{:,1});
  %   source{2} = ft_sourcegrandaverage(cfg, source_cell{:,2});
  %   source{3} = ft_sourcegrandaverage(cfg, source_cell{:,3});
  %   source{4} = ft_sourcegrandaverage(cfg, source_cell{:,4});
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  source_cell_mean{iage} = ft_sourcegrandaverage(cfg, source_cell_mean{iage}{:});
  source_cell_sdbetween{iage} = ft_sourcegrandaverage(cfg, source_cell_sdbetween{iage}{:});
  source_cell_sdwithin{iage} = ft_sourcegrandaverage(cfg, source_cell_sdwithin{iage}{:});
  source_cell_sdall{iage} = ft_sourcegrandaverage(cfg, source_cell_sdall{iage}{:});
end
%avg over subjects
% for iage = 1:2
%   source_cell_mean{iage}.pow = mean(source_cell_mean{iage}.pow)';
%   source_cell_sdbetween{iage}.pow = mean(source_cell_sdbetween{iage}.pow)';
%   source_cell_sdwithin{iage}.pow = mean(source_cell_sdwithin{iage}.pow)';
%   source_cell_sdall{iage}.pow = mean(source_cell_sdall{iage}.pow)';
%   
%   source_cell_mean{iage}.dimord = 'pos';
%   source_cell_sdbetween{iage}.dimord = 'pos';
%   source_cell_sdwithin{iage}.dimord = 'pos';
%   source_cell_sdall{iage}.dimord = 'pos';
% end

% average over age groups 
cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
source_cell_mean{3} = ft_sourcegrandaverage(cfg, source_cell_mean{:});
source_cell_sdbetween{3} = ft_sourcegrandaverage(cfg, source_cell_sdbetween{:});
source_cell_sdwithin{3} = ft_sourcegrandaverage(cfg, source_cell_sdwithin{:});
source_cell_sdall{3} = ft_sourcegrandaverage(cfg, source_cell_sdall{:});

cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'pow';
% source{5} = ft_math(cfg, source{2}, source{3});
source_cell_mean{4} = ft_math(cfg, source_cell_mean{1:2});
source_cell_sdbetween{4} = ft_math(cfg, source_cell_sdbetween{1:2});
source_cell_sdwithin{4} = ft_math(cfg, source_cell_sdwithin{1:2});
source_cell_sdall{4} = ft_math(cfg, source_cell_sdall{1:2});



save source_average source_cell_mean source_cell_sdbetween source_cell_sdwithin source_cell_sdall


%%
load colormap_jetlightgray.mat

close all
cfg=[];
cfg.funparameter = 'pow';
cfg.method = 'slice'; % slice ortho glassbrain vertex
cfg.funcolormap = cmap;
cfg.funcolorlim = 'maxabs';
% [~,i]= max(source{1}.pow(:));
% cfg.location = [-30 -90 18]; % set at max mean response voxel, 130361
% cfg.location = [-6 -78 51]; % set at max SDall response voxel, 178113
cfg.location = [0 -75 45]; % set at max SDall response voxel, 178113
ft_sourceplot(cfg, source_cell_sdwithin{1});   title('SD within YA')
ft_sourceplot(cfg, source_cell_sdwithin{2});   title('SD within OA')
ft_sourceplot(cfg, source_cell_sdwithin{4});   title('SD within YA-OA')
ft_sourceplot(cfg, source_cell_sdbetween{1});   title('SD between YA')
ft_sourceplot(cfg, source_cell_sdbetween{2});   title('SD between OA')
ft_sourceplot(cfg, source_cell_sdbetween{4});   title('SD between YA-OA')
ft_sourceplot(cfg, source_cell_mean{1});   title('mean YA')
ft_sourceplot(cfg, source_cell_mean{2});   title('mean OA')
ft_sourceplot(cfg, source_cell_mean{4});   title('mean YA-OA')

cfg.funcolorlim = 'zeromax';
cfg.funcolormap = cmap(129:end,:);
ft_sourceplot(cfg, source{2});   title('SDbetweentrials')

cfg.funcolorlim = 'zeromax';
ft_sourceplot(cfg, source{3});   title('SDwithintrials')

cfg.funcolorlim = 'zeromax';
ft_sourceplot(cfg, source{4});   title('SDfull')

%% plot time course example subject
% voxelind = 130361;
% voxelind = 178113;
voxelind = 169531;
voxelind = 170624;
voxelind = 75764;
source_subj = load('source_sub-35.mat');
source_subj.time = 0:4;

f = figure; subplot(2,2,1); plot(source_subj.pow(voxelind,:)); hold on
title(sprintf('Time series 150 trials voxel %d', voxelind))

subplot(2,2,2); plot(source_subj.time, squeeze(source_subj.pow(voxelind,:,1:5))); hold on
plot(source_subj.time, mean(squeeze(source_subj.pow(voxelind,:,1:5))), 'Linewidth', 2, 'Color', [0 0 0]);
title(sprintf('Single trial time series'))

subplot(2,2,3); plot(source_subj.time, squeeze(source_subj.pow(voxelind,:,1:5)) - ...
  mean(squeeze(source_subj.pow(voxelind,:,1:5)),2)); hold on
title(sprintf('within-trial variability'))

subplot(2,2,4); plot(source_subj.time, ...
  repmat(mean(squeeze(source_subj.pow(voxelind,:,1:5)),2), [1 5] )); hold on
title(sprintf('between-trial variability'))

saveas(f, fullfile(PREOUT, sprintf('voxeltimeseries_%d.pdf', voxelind)))




