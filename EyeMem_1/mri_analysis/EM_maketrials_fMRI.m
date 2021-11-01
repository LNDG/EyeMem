function EM_maketrials_fMRI(cfg)
% select vois after creating trials
%
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
latency = cfg.latency;
% eventshift = cfg.eventshift;

disp(mrifile)
disp 'read mri from file'
tic
mri = ft_read_mri(mrifile);
toc

disp 'convert volume to source' %this puts the transform field in pos, so mm coord stay intact
source = ft_checkdata(mri, 'datatype', {'source'}, 'feedback', 'yes', 'hasunit', 'yes');

disp 'rename anatomy to pow'
source = renameStructField(source , 'anatomy', 'pow');

% disp 'fix inside field: make false outside the brain'
% source.inside = true(size(source.inside,1),1);
% source.inside(mean(source.pow,2) == 0) = false;

disp 'only include voxels common to all subjects'
if ismac
  load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/common_coords.mat
else
  load /home/mpib/kloosterman/projectdata/eyemem/variability/ftsource/common_coords.mat
end  
% source.inside = repmat(common_coords, 1, size(source.inside,2)); % not
% needed for every time point
source.inside = common_coords; % pow not yet updated

% disp 'reshape to 150 trials x 5 TRs'
% source.pow = reshape(source.pow, [], 5, 150);
% source.pow =  permute(source.pow, [1 3 2]); % flip last 2 dims
% source.powdimord = 'pos_rpt_time';
% source.dim = source.dim(1:3);
%%
if ismac
  standardsfolder = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/A_standards';
  anat = ft_read_mri(fullfile(standardsfolder, 'MNI152_T1_3mm_brain.nii.gz' ));
  
%   disp 'stick stats to anatomy'
%   cfg=[];
%   cfg.parameter = {'stat' 'mask'};
%   [interp] = ft_sourceinterpolate(cfg, stat, anat);

  
  %   atlas = ft_read_atlas('TTatlas+tlrc.HEAD'); % TTatlas+tlrc.HEAD
  tmp = source;
  tmp.powdimord = 'pos';
  %   tmp.anatomy = tmp.anatomy(:,:,:,vol);
    tmp.pow = std(tmp.pow(:,:),0,2);
%   tmp.pow = mean(tmp.pow(:,:),2);
  tmp.dim = tmp.dim(1:3);
  cfg=[];
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  cfg.funparameter = 'pow';
  %   cfg.funcolorlim = [-300 300];
  %   cfg.atlas = atlas;
  ft_sourceplot(cfg, tmp)
end

%%
load(eyefile); % data, has events

disp 'update triggers so they match the concatenated data'
nvolsperrun = 462;  % 2310/462 = 5;
onsetshift = (data.trialinfo(:,12)-1) .* nvolsperrun;  % runno-1 * 462 = shifts due to concatenation
piconsets = round(data.trialinfo(:,9)) + onsetshift; % see EMcheckfMRItimings for event timing checks

if isnumeric(latency)
  % 1. Make stim-locked timecourses ("ERPs")
  begtim = -4;
  endtim = 15;
  ntps = diff([begtim endtim])+1;
  trl = [piconsets+begtim piconsets+endtim];
  
  sourcetrl = source;
  sourcetrl.pow = NaN([size(sourcetrl.pow,1), length(trl), ntps]);
  sourcetrl.powdimord = 'pos_rpt_time';
%   sourcetrl.dim = [60 72 60 150 ntps]; % only dims 1:3???
  sourcetrl.dim = [60 72 60]; % only dims 1:3???
  sourcetrl.time = begtim:endtim;
  for itrial = 1:length(trl)
    startind = trl(itrial,1);
    endind = trl(itrial,end);
    if endind > size(source.pow,2)
      endind = size(source.pow,2);
    end
    sourcetrl.pow(:, itrial, 1:length(startind:endind)) = source.pow(:, startind:endind);
  end
  sourcetrl.trialinfo = data.trialinfo;
  
  %% plotting trial time courses
  if ismac
    avg = mean(sourcetrl.pow(sourcetrl.inside,:,:),2);
    sd = squeeze(std(avg));
    avg = squeeze(mean(avg));
    figure; shadedErrorBar(begtim:endtim, avg, sd/sqrt(149))
    
    % dat = source.pow(source.inside,:);
    % figure; histogram(dat(1,:));
    voxoi = 108754;
    vox = nanmean(sourcetrl.pow(voxoi,:,:),2);
    voxsd = nanstd(sourcetrl.pow(voxoi,:,:),0,2);
    figure; shadedErrorBar(begtim:endtim, squeeze(vox), squeeze(voxsd/sqrt(149)));
    grid on; set(gca, 'XTick', -20:20); title(voxoi)
    
    vox = squeeze(sourcetrl.pow(voxoi,:,:));
    figure; plot(begtim:endtim, vox(1:1:end,:));
    grid on; set(gca, 'XTick', -20:20)
    title(subjno)
    %% plot of avg or std trial
    tmp = sourcetrl;
    tmp.pow = squeeze(mean(tmp.pow, 2));
%     tmp.pow = squeeze(std(tmp.pow, 0, 2));
    tmp.powdimord = 'pos_time';
    tmp.dim = tmp.dim(1:3);
    %%
    cfg=[];
    cfg.funparameter = 'pow';
    cfg.method = 'ortho'; % slice ortho glassbrain vertex
    load colormap_jetlightgray.mat
%     cfg.funcolormap = cmap(129:end,:);
    cfg.funcolormap = cmap;
    cfg.funcolorlim = 'zeromax';
    cfg.funcolorlim = 'maxabs';
%     cfg.funcolorlim = [0 150]% 'zeromax'
%     cfg.colorbar = 'yes'
%     ft_sourceplot(cfg, tmp, anat)
    ft_sourceplot(cfg, tmp)
    %%
  end
  
else
  sourcetrl = source;
end

if ~ismac; clear source; end

sourcetrl.pow = single(sourcetrl.pow);
disp(outfile);
disp(sourcetrl)
save(outfile, 'sourcetrl', '-v7.3')
% save(outfile, '-struct', 'source') % save as separate variables

%%
% % deconvolution tryout: just fancy averaging only helps with several, jittered events
% y = source.pow(87158,:)'; % visual cortex voxel with high SD
% ntp = 8;
% X = zeros(length(y), ntp+1);
% X(:,end) = 1; % constant term
% for i = 1:ntp
%   X(piconsets+i-1, i) = 1;
% end
% figure; imagesc(X); ylim([1 50])
% b = regress(y, X)
% figure; plot(b(1:ntp))

