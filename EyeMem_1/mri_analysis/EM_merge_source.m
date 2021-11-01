function source = EM_merge_source()
% merge source made by EM_selectVOIs_fMRI to compute avg, plot, etc

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%   backend = 'slurm';
      backend = 'torque';
  %   backend = 'local';
  compile = 'no';
end

PREIN = fullfile(basepath, 'variability', 'VOIsel', 'ERPremoved');
cd(PREIN)

source = [];
agedirs = {'YA' 'OA' 'allsubj'};
sourceall = cell(2,1);
for iage = 1:2
  cd(fullfile(PREIN, agedirs{iage}))
  
  % SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?
  SUBJ = dir('source_sub-*');

  for isub = 1:length(SUBJ)
    disp(SUBJ(isub).name)
    [~,a]= fileparts(SUBJ(isub).name);
    a = tokenize(a, '_');
    source.SUBJ{iage}(isub,1) = a(2);
    
    load(SUBJ(isub).name)
    sourceall{iage}(isub,:,:) = num2cell(src.source);
  end
end

%%
disp 'collect subjects'
% clear temp
cfg=[];
cfg.keepindividual = 'yes';
for iage = 1:2
  for imeas = 1:size(src.source,1)
    for ibin = 1:size(src.source,2)
      temp = ft_sourcegrandaverage(cfg, sourceall{iage}{:, imeas, ibin} );
      source.source(iage,imeas,ibin) = rmfield(temp, 'unit');
    end
  end
end
clear sourceall
%%
%  TODO average over subj first, then contrast (indep samples...)
% disp 'YA-OA contrast'
% cfg=[];
% cfg.parameter = 'pow';
% cfg.operation = 'subtract';
% temp(3,:) = arrayfun(@(x,y) ft_math(cfg, x,y), temp(1,:), temp(2,:));


% disp 'bin3-bin1 contrast'
% cfg=[];
% cfg.parameter = 'pow';
% cfg.operation = 'subtract';
% % temp(:,:,5) = arrayfun(@(x,y) ft_math(cfg, x,y), temp(:,:,3), temp(:,:,3));
% source.source(:,:,5) = arrayfun(@(x,y) ft_math(cfg, x,y), source.source(:,:,3), source.source(:,:,1));

% deal with output
source.SUBJ{3} = vertcat(source.SUBJ{:});
source.dimord = ['age_' src.dimord];
source.age = agedirs;
source.meas = src.meas; %{'avg' 'std' 'mse'};
source.bin = src.bin; %{'avg' 'std' 'mse'};
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
source.behavior = behav;

% outfile = fullfile(PREIN, 'source_all.mat');
% disp(outfile)
% save(outfile)

