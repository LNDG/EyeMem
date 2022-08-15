function EM_mri_to_ftsource(cfg)
% put the Eyemem mri data into fieldtrip source. This will (1) make creating the common
% coordinates easier (using pos and inside), and (2) allows rpt and time dims, so trials and time
% points can be easily selected for subsequent PLS analysis.
% inside should come from a Gray Matter mask, so it is the same for all
% subjects.

% example in ft_datatype_source:
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [6732x1 logical]      boolean vector that indicates at which positions the source activity was estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          time: 0.100                 the latency at which the activity is estimated (in seconds)
%          freq: 30                    the frequency at which the activity is estimated (in Hz)
%           pow: [6732x120 double]     the estimated power at each source position
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 repetitions (i.e. trials)
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure

% 
% selecting inside voxels and put them back after PLS:
% IND = find(inside) gives the IND of pos that are included. pow(IND) =
% BSR from compare_u field then puts it back.

close all

mrifile = cfg.mrifile;
outfile = cfg.outfile;
subjno = cfg.subjno;
eyefile = cfg.eyefile;
nTRpertrial = cfg.nTRpertrial;

disp 'read mri from file'
mri = ft_read_mri(mrifile);

% 1st the data. TR's of the last trial might be missing: add nan trial up to 150
disp 'Fill up last trial with nans'
nTRreq = nTRpertrial*150;
nTR = mri.dim(4);
if nTR ~= nTRreq % TR's missing
  nTRmissing = nTRreq - nTR;
  nanTR = nan( [ mri.dim(1:3) nTRmissing] );
  mri.anatomy = cat(4, mri.anatomy, nanTR);
  mri.dim = [mri.dim(1:3) nTRreq];
end

if ismac
  vol=300;
  tmp = mri;
  tmp.dim = tmp.dim(1:3);
%   tmp.anatomy = tmp.anatomy(:,:,:,vol); 
  tmp.anatomy = nanmean(tmp.anatomy,4); 
  cfg=[];
  cfg.method = 'ortho'; % slice
  ft_sourceplot(cfg, tmp)
end

disp 'convert volume to source' %this puts the transform field in pos, so mm coord stay intact
source = ft_checkdata(mri, 'datatype', {'source'}, 'feedback', 'yes', 'hasunit', 'yes');

disp 'rename anatomy to pow'
source = renameStructField(source , 'anatomy', 'pow');

disp 'fix inside field: make false outside the brain'
source.inside = true(size(source.inside,1),1);
source.inside(mean(source.pow,2) == 0) = false;

disp 'reshape to 150 trials x nTRpertrial TRs'
source.pow = reshape(source.pow, [], nTRpertrial, 150); 
source.pow =  permute(source.pow, [1 3 2]); % flip last 2 dims
source.powdimord = 'pos_rpt_time';  
source.dim = source.dim(1:3);

if ismac
  vol=300;
  tmp = source;
  tmp.powdimord = 'pos';
%   tmp.anatomy = tmp.anatomy(:,:,:,vol); 
  tmp.pow = nanmean(tmp.pow(:,:),2); 
  cfg=[];
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  cfg.funparameter = 'pow';
  cfg.funcolorlim = [-1000 1000]; % [-300 300]
  ft_sourceplot(cfg, tmp)

end

load(eyefile);
source.trialinfo = data.trialinfo;
% source.cfg = mri.hdr;

disp(outfile); disp(source)
% save(outfile, 'source') 
save(outfile, '-struct', 'source') % save as separate variables

% if ismac
%   disp 'Kolmogorov-Smirnov test for normality'
%   [h,p] = kstest(seldat(seldat~=0));
%   figure; histogram(seldat(seldat~=0))
%   title('Example subject, all brain voxels')
% end
