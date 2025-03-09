% function EM_GLMsingle()

path = '/Users/kloosterman/Documents/GitHub/GLMsingle';
cd(path); addpath(path); setup;

eyefile = '/Users/kloosterman/projectdata/eyemem/preproc/eye/YA/eye_sub-09.mat';
eye_data = load(eyefile);
stimtimes = round(reshape(eye_data.trialinfo.fMRI_trigger, 30, 5));
RT = round(eye_data.trialinfo.RT_quizpic/1000);

data = {}; design = {}; design2 = {}; extraregressors = {};
for irun = 1:5
  niftifile = sprintf('/Users/kloosterman/projectdata/eyemem/data/B_data/sub-09/func/sub-09_task-eyemem_run-0%d_bold.nii.gz', irun);
  mri = ft_read_mri(niftifile);
  mri.anatomy = single(mri.anatomy);
  mri_data = mri.anatomy(:,:,:,13:end);
  mri_mirror = mri.anatomy(:,:,:,end-11:end);
  mri_mirror = flip(mri_mirror,4);
  data{end+1} = cat(4, mri_data, mri_mirror);

  designmat = false(size(data{1},4), size(stimtimes,1)); 
  designmat2 = false(size(data{1},4), size(stimtimes,1)*5); 
  extra_regr = zeros(size(data{1},4), 2); 
  ltr = 0;
  for ipic = 1:30
    % designmat(stimtimes(ipic,irun):stimtimes(ipic,irun)+4, ipic) = true; % model 5 TRs as trials
    designmat(stimtimes(ipic,irun):stimtimes(ipic,irun)+0, ipic) = true; % model 1st TR as trial
    % designmat(stimtimes(ipic,irun)+4, ipic) = true; % model 5th TR as trial
    % designmat(stim_rnd(ipic):stim_rnd(ipic)+4,ipic) = true; % model 5 TRs per trial 1s shift
    % designmat(stim_rnd(ipic),ipic) = true; % single trial beta
    % if ipic ~= 30
      extra_regr(stimtimes(ipic,irun)+4+2,1) = true; % Quiz pic
      extra_regr(stimtimes(ipic,irun)+4+2+RT,2) = 1; % Quiz pic
    % end
    for itr = 0:4
      designmat2(stimtimes(ipic,irun)+itr, ltr+1) = true; % model 5 TRs as trials
      ltr = ltr+1;
    end

  end
  design{end+1} = designmat;
  design2{end+1} = designmat2;
  extraregressors{end+1} = extra_regr;
end

figure; tiledlayout(1,5)
for d = 1:5  
  nexttile
  imagesc(design{d}); colormap gray; drawnow
  xlabel('Conditions')
  ylabel('TRs')
  title(sprintf('Design matrix for run%i',d))
  axis image
end

% figure(2);clf
% imagesc(data{1}(:,:,1,1));
% colormap(gray);
% axis equal tight;
% c=colorbar;
% title('fMRI data (first volume)');
% set(gcf,'Position',[418   412   782   605])
% axis off
% c.Label.String = 'T2*w intensity';
% set(gca,'FontSize',15)


%%
opt = struct('wantmemoryoutputs',[0 0 0 0], 'wantfracridge', 0, 'wantglmdenoise', 0, 'wantlibrary', 1, 'chunknum', prod(mri.dim));
% opt.extraregressors = extraregressors';
stimdur = 5;
tr = 1;
outputdir = '/Users/kloosterman/projectdata/eyemem/preproc';
tic 
[results] = GLMestimatesingletrial(design,data,stimdur,tr,[outputdir '/GLMsingle5TRs5TRspertrial_separateconds'],opt);
toc

%% TODO compare GLMsingle with LSS

% lss_file = '/Users/kloosterman/projectdata/eyemem/variability2/5TRspertrial/ftsource/source_sub-09.mat';
% 
% source = load(lss_file)