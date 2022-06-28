function EM_hmax_fmrislope
% correlate hmax index with SDBOLD to see if SDBOLD tracks stimulus
% complexity

basepath = '/Users/kloosterman/beegfs/projectdata/eyemem';
cd(basepath)

PREIN_mri = fullfile(basepath, 'variability', 'all_stimuli_5sres');
PREIN_eye = fullfile(basepath, 'preproc', 'eye');

% young:
subjlist_YA = {
  
'sub-07'    'sub-09'    'sub-10'    'sub-11'    'sub-15'    'sub-17'    'sub-18'    'sub-19'    'sub-20'    'sub-22'    'sub-31'    'sub-33' ...
'sub-39'    'sub-44'    'sub-47'    'sub-49'    'sub-50'    'sub-51'    'sub-52'    'sub-53'    'sub-55'    'sub-56'    'sub-57'    'sub-58' ...
'sub-59'    'sub-62'    'sub-64'    'sub-65'    'sub-67'    'sub-68'    'sub-69'    'sub-70'    'sub-71'    'sub-72'    'sub-73'    'sub-74' ...
'sub-76'    'sub-77'    'sub-78'    'sub-80'    'sub-83'    'sub-90'    'sub-92'    'sub-93'    'sub-95'    'sub-97'    'sub-98'};

subjlist = [];
subjlist.YA = dir(fullfile(PREIN_mri, 'YA', 'sub-*'));
subjlist.OA = dir(fullfile(PREIN_mri, 'OA', 'sub-*'));

niicorr = cell(1,2);
niicorr{1} = nan(length(subjlist.YA),60,72,60);
niicorr{2} = nan(length(subjlist.OA),60,72,60);
% niicorr = nan(nsub,2,60,72,60); % dimord: sub age x y z 
% dprime.test =

ageleg = {'YA' 'OA'};
corr_hmaxfmrislope = [];
for iage=1:2
  cursubjlist = subjlist.(ageleg{iage});
  nsub = length(cursubjlist);
  for isub = 1:nsub
    
    disp(cursubjlist(isub).name)
    %   cd(subjlist(isub).name)
    
    niilist = dir( fullfile(cursubjlist(isub).folder, cursubjlist(isub).name, '*.nii'));
    
    clear niimat
    for in = 1:length(niilist) % load mri
      nii = load_untouch_nii( fullfile(niilist(in).folder, niilist(in).name));
      niimat(in,:,:,:) = nii.img;
    end
    % load behav
    load( fullfile(PREIN_eye, ageleg{iage}, sprintf('eye_%s', cursubjlist(isub).name)) )
    runinfo = ft_findcfg(data(1).cfg, 'runinfo');
    behavior = runinfo.behavior;
    dprime.test{iage}(isub,1) = nanmean(behavior.test.dprime,2);
    
    corrtype = 'Pearson'; % corr_col only does pearson for now
    niimat2d = reshape(niimat, in, []);
    % slow way using corr
    %   niicorr2d = nan(1,size(niimat2d,2));
    %   for icol = 1:size(niimat2d,2)
    %     if all(~isnan(niimat2d(:,icol)))
    %       %     disp(icol)
    %       niicorr2d(1,icol) = corr(niimat2d(:,icol), transpose(1:length(niilist)), 'type', corrtype);
    %     end
    %   end
    
    % use corr_col, faster
    niicorr2d = corr_col(niimat2d, repmat(transpose(1:length(niilist)), 1, size(niimat2d,2) ));
    
    sz = size(niimat);
    %     if ismember(subjlist(isub).name, subjlist_YA)
    %       %     niicorr(isub,1,:,:,:) = reshape(niicorr2d, sz(2:4));
    %       %     niicorr{1} = [niicorr{1}   reshape(niicorr2d, sz(2:4))];
    %     else
    %       niicorr(isub,2,:,:,:) = reshape(niicorr2d, sz(2:4));
    %     end
    niicorr{iage}(isub,:,:,:) = reshape(niicorr2d, sz(2:4)); %% nan(length(subjlist.YA), 60,72,60); % niicorr{1} = nan(length(subjlist.YA),60,72,60);

    cd ..
  end
  % correlate hmax vs sdbold corr with dprime across subj: best adapters
  % perform best?
  niicorr2d = reshape(niicorr{iage}, nsub, []);
  temp = corr_col(niicorr2d, repmat(dprime.test{iage}, 1, size(niicorr2d,2)));
  sz = size(niicorr{iage});
  corr_hmaxfmrislope(iage,:,:,:) = reshape(temp, sz(2:4));
end

% niicorrmean = squeeze(nanmean(niicorr, 1));
% niicorrmean(3,:,:,:) = mean(niicorrmean);
% niicorrmean(4,:,:,:) = niicorrmean(1,:,:,:) - niicorrmean(2,:,:,:);
niicorrmean = cellfun(@(x) squeeze(nanmean(x, 1)), niicorr, 'uni', 0 );
% niicorrmean{3} = % TODO avg young and old
niicorrmean{4} = niicorrmean{1} - niicorrmean{2};

%% plot hmax vs sdBOLD corr
close all
load /Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/plotting/colormap_jetlightgray.mat
for iage = [1,2,4] % 1:4
  f = figure; f.Position = [ 680         562        1000         500];
  %   colormap(cmap)
  iplot=0;
  for islice = 10:5:50
    iplot=iplot+1;
    subplot(3,3,iplot)
%     if iage == 3
%       imagesc(squeeze(niicorrmean(iage,:,:,islice)), [-0.025 0.025]);
%     else
%       imagesc(squeeze(niicorrmean(iage,:,:,islice)), [-0.05 0.05]);
%     end
    if iage == 3
      imagesc(squeeze(niicorrmean{iage}(:,:,islice)), [-0.025 0.025]);
    else
      imagesc(squeeze(niicorrmean{iage}(:,:,islice)), [-0.05 0.05]);      
    end
    colorbar
  end
end
title(corrtype)

%% plot hmax vs sdBOLD corr vs dprime
close all
load /Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/plotting/colormap_jetlightgray.mat
for iage = 1:2
  f = figure; f.Position = [ 680         562        1000         500];
  %   colormap(cmap)
  iplot=0;
  for islice = 10:5:50
    iplot=iplot+1;
    subplot(3,3,iplot)

    imagesc(squeeze(corr_hmaxfmrislope(iage,:,:,islice)), [0.25 0.3]);  % [-0.25 0.25]
    
    colorbar
  end
end
title(corrtype)
%% save out nifti
outnii = nii;
for iage = 1:4
  nii.img = squeeze(niicorrmean{iage});
  save_untouch_nii(nii, fullfile(PREIN_mri, sprintf('hmax_fmrislope%d',iage)))% sprintf('hmax_fmrislope%d',iage))
end
  cd(PREIN_mri)

%% save out nifti fmrislope corr
outnii = nii;
for iage = 1:2
  nii.img = squeeze(corr_hmaxfmrislope(iage,:,:,:));
  save_untouch_nii(nii, fullfile(PREIN_mri, sprintf('hmax_fmrislope_vs_dprime%d',iage)))
end
  cd(PREIN_mri)

%regress
% corr(niimat())
