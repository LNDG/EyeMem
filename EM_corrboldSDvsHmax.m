function EM_corrboldSDvsHmax()
% load beta niftis per trial, bin and correlate with bin number
% TODO run on cluster per subj?
PREIN = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/SDBOLDpertrial/singletrial_5sres';
addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/NIFTI_toolbox')

PREOUT = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/PLS/HMAXvsSDBOLD';
cd(PREIN)
subjlist = dir('sub*');
plotit = 0;

sd_corr = [];
for isub = 1:length(subjlist)
  cd(subjlist(isub).name)
  cd('beta_series')
  niilist = dir('*_beta_series.nii');
  ntrials = length(niilist);
  im_append = [];
  ibin = 0;
  for itrial = 1:ntrials
    %     load nii    
    nii = dir(sprintf('*trial%d_beta_series.nii', itrial));
    disp(nii.name)
    niidat = load_nii(niilist(itrial).name);
    im_append = cat(4, im_append, niidat.img);
    % append nii
    if mod(itrial,5) == 0 || itrial == ntrials % 5 trials at a time or last trial
      disp('get sd')
      ibin = ibin+1;
      % save SDnii per cond (bin)
      outimg = std(im_append,1,4); % get SD over 25 samples, SD weighted by N observations 
      niiout = make_nii(outimg, [3 3 3]);
      [~, niiname] = fileparts(nii.name);
      save_nii(niiout, fullfile(PREOUT, [niiname '_SD.nii'] ))
      im_append = []; % empty for the next 5 trials

      if plotit
        figure; ipl=0;
        for i=10:3:50
          ipl=ipl+1;
          subplot(4,4,ipl)
          imagesc(squeeze(niiout.img(:,:,i))); % , [0 1500]
          colorbar
        end
      end
      
    end
  end
  
  % TODO save nifti with SDbold per HMAX bin -> put in PLS!!!
  
%   % correlate for each voxel across bins the bin no with SD
%   siz= size(im_sd);
%   bin_ind_mat = repmat(transpose(1:30), 1, siz(2:4));  size(bin_ind_mat)
%   bin_ind_mat = permute(bin_ind_mat, [1 3 4 2]); % size(bin_ind_mat)
%   sd_corr(isub,:,:,:) = reshape(corr_col(im_sd(:,:), bin_ind_mat(:,:)), siz(2:4));
%   if plotit    
%     figure; ipl=0;
%     for i=10:3:50
%       ipl=ipl+1;
%       subplot(4,4,ipl)
%       imagesc(squeeze(sd_corr(isub,:,:,i)), [-1 1]);
%     end
%   end
  
  cd(PREIN)
  
end




