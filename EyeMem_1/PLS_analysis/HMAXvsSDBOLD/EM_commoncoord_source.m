function EM_commoncoord_source(cfg)
% run from runMIBmeg_analysis

subjlist = cfg.subjlist;
PREIN = cfg.PREIN;

disp 'Match the inside bool aka common coordinates'

cd(subjlist(1).folder)

allsource = {};
for isub = 1:length(subjlist)
  disp(subjlist(isub).name)
  source = load(subjlist(isub).name, 'inside');
  allsource{end+1} = source;
end
source = cell2mat(allsource);
clear allsource

disp 'Match the inside bool aka common coordinates'
inside_allsubj = [source.inside];
common_coords = all(inside_allsubj,2);
[source.inside] = deal(common_coords);  % cast this to all subj source, note that pow field is not updated

disp 'Apply Gray matter mask'
standardsfolder = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/A_standards/tissuepriors/';
GMmask = ft_read_mri(fullfile(standardsfolder, 'avg152T1_gray_MNI_3mm.nii' ));

source2=load(subjlist(isub).name);
source2.pow = mean(source2.pow(:,:),2);
source2.powdimord = 'pos';

disp 'stick mask to anatomy'
cfg=[];
cfg.parameter = {'anatomy'};
[interp] = ft_sourceinterpolate(cfg, GMmask, source2);

disp 'binarize again'
% figure;histogram(interp.anatomy(:))
interp.anatomy(interp.anatomy <= 0.5) = 0;
interp.anatomy(interp.anatomy >= 0.5) = 1;
interp.anatomy = interp.anatomy == 1; % make bool

disp 'convert mask volume to source' %this puts the transform field in pos, so mm coord stay intact
source_mask = ft_checkdata(interp, 'datatype', {'source'}, 'feedback', 'yes', 'hasunit', 'yes');

sum(common_coords)
common_coords = common_coords & source_mask.anatomy;
sum(common_coords)
source2.inside = common_coords;

if ismac
  disp 'plotting'
  % close all
  cfg=[];
  cfg.method = 'slice'; % slice ortho glassbrain vertex surface
  cfg.funparameter = 'pow';
  cfg.maskparameter = 'inside';
  cfg.anaparameter = 'anatomy';
  ft_sourceplot(cfg, source2, interp)
end


disp 'save common_coords'
save(fullfile(PREIN, 'common_coords.mat'), 'common_coords')
% [source.inside] = deal(common_coords);  % cast this to all subj source, note that pow field is not updated
