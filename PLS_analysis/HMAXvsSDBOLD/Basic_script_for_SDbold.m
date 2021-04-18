addpath /software/jimmy/plsgui;
addpath('/home/randy/matlab/old_pls/');
addpath /arrakis/mcintosh_lab/natasa/parse_melodic;
addpath /arrakis/mcintosh_lab/natasa/lego_matlab;
workdir = '/iris/grady_lab/charisa/doug/';
pls_dir = '/iris/grady_lab/charisa/doug/PLS/';
atlasdir = '/flora/mcintosh_lab/natasa/MultiTask/CommonAtlas/nl/';


group_names = {'Young','Old'};%set group structure
patterns = {'Y','O'};%set key character for determining groups
cd(workdir);
for g=1:numel(group_names)
  junk = dir([patterns{g} '*']);%not sure what this is for...identifying when things are wrong? 'dir' returns multiple bits of info about a file labelled as junk (e.g., filename, data, directory, etc).
  groups{g} = {junk.name};
end
cd(pls_dir);%watch as change from workdir to pls_dir....cross reference to current folder structure on RRI servers.
subjects = [groups{1},groups{2}];

conditions = {'fix','cond1','cond2','cond3','cond4','cond5','cond6','cond7'};%set all relevant condition names

%maskfile = [atlasdir 'avg4_4mm_brainmask'];
maskfile = [atlasdir 'avg4_4mm_GM_mask'];%determine the location and name of the mask file

%%% get st_coords from the mask
mask = load_img(maskfile,'l');%look for anything in the mask that is a 1...denotes a logical index for valid values. Must get this load_img command from Natasa though.
st_coords = reshape(find(mask),1,[]);%this command reshapes into a single row vector, covering all elements (denoted by '[]'), and the 'find' returns only the cell/position number of all valid values (i.e., =1 in the mask).
%for example...if cell positions (i.e., linear index positions..which
%always go down the columns) 1,2,4, and 6 had "1" values in the mask, then
%the reshape would return a single row vector with nothing but values
%1,2,4,6 in it!
mask =  load_img([atlasdir 'wm_roi_4mm'],'l');
wm_coords = reshape(find(mask),1,[]);
mask =  load_img([atlasdir 'csf_roi_4mm'],'l');
csf_coords = reshape(find(mask),1,[]);
clear mask;


% $$$ datatypes = {'raw_smooth','deno_smooth'}; 
% $$$ vartypes = {'mean1_allcorrected','stdev_allcorrected'};

datatypes = {'deno_smooth'}; %this datatype/vartype structure is nice if we
%are running multiple analysis pipeline types, etc. If not, could abandon
%these loops.
vartypes = {'mean1_allcorrected'};

for dt = 1:numel(datatypes)%these successive loops work on datatypes, 
    %vartypes, and subjects. If have only one type, just loop over subjs.
  datatype = datatypes{dt};%datatypes was defined as a cell array above.
  %Now, 'datatypes{dt}' returns the value from that cell array...e.g., 
  %deno_smooth is in position 1, thus is also first in the loop, and then 
  %we will use it below for naming.
  for vt = 1:numel(vartypes)
    vartype = vartypes{vt};
    pattern = [vartype '_' datatype];%this establishes same naming for all
    %subjects...see below.
    for s=1:numel(subjects)
      subj = subjects{s};
      %if exist([subj '_' pattern '_BfMRIdatamat.mat'],'file'), continue;end
      
      a = load([workdir '/charisaPLS/' subj '_BfMRIsession.mat']);%this
      %loads a particular session file.
      a.session_info.pls_data_path = pls_dir;%this will add info to 'a' for
      %tracking purposes, and for commands below.
      a.session_info.datamat_prefix = [subj '_' pattern];%stores common
      %datamat prefix
      a.session_info.file_pattern = '*nii';%
      
      switch datatype
      
       case 'raw_smooth'
        for run=1:a.session_info.num_runs%everything following a. is part 
            %of the pls session file structure, and can thus be called
            %here.
          a.session_info.run(run).data_path = [workdir subj];
          a.session_info.run(run).data_files = {['smooth_CommonAtlas_mc_tshift_cond' num2str(run) '.nii']};
          a.session_info.run(run).file_pattern = ['smooth_CommonAtlas_mc_tshift_cond' num2str(run) '.nii*'];
        end
      
       case 'deno_smooth'
        for run=1:a.session_info.num_runs
          a.session_info.run(run).data_path = [workdir subj '/melodic_CommonAtlas_mc_tshift_cond' num2str(run) '.ica'];
          a.session_info.run(run).data_files = {'smooth_pipeline_denoised_data.nii'};
          a.session_info.run(run).file_pattern = 'smooth_pipeline_denoised_data.nii*';
        end
      end
      save([pls_dir subj '_' pattern '_BfMRIsession.mat'],'-struct','a','-mat');%essentially, we have 
      %just added some infromation to the session files that were already
      %created, and saved them again here with new naming convention 
      %accoridng the the loops above. Note that original session files must
      %follow appropriate naming structure as evident in line 61 above. Now
      %onto to datamat prep and creation.
      
      %%%% first set up common fields for datamats - same acrosss subjects.
      %these come from batch PLS I believe.
      clear tmp;
      tmp.create_ver = '5.0807231';%not sure where value comes from; check.
      tmp.st_win_size = 1;%PLS command...must check.
      tmp.st_coords = st_coords;%these are defined above as a mask of voxels
      %to use. All others are excluded (e.g., in a GM mask).
      tmp.st_origin = [0 0 0];%sets the origin...PLS default?
      tmp.st_voxel_size = [4 4 4];%vox size defined
      tmp.st_dims = [40 48 1 34];%dimensions of matrix defined for PLS
      tmp.st_evt_list = [1:numel(conditions)];%how many conditions?
      tmp.SingleSubject = 0;%not a single subj design
      tmp.behavname = {};%no behav names or data just yet, but set up the
      %field anyway...
      tmp.behavdata = [];
      tmp.normalize_volume_mean = 0;%we do our own normalizations below
      tmp.create_datamat_info.brain_mask_file = [maskfile '.img'];%brain
      %mask identified (defined above)
      tmp.create_datamat_info.brain_coord_thresh = 0;%no ROIs here...check
      tmp.create_datamat_info.consider_all_voxels_as_brain = 0;%we don't...
      tmp.create_datamat_info.num_skipped_scans = 0;%no skipped scans
      tmp.create_datamat_info.run_idx = [1:7];%seven runs here
      tmp.create_datamat_info.ignore_slices = [];%no slices to ignore, so 
      %keep all
      tmp.create_datamat_info.normalize_volume_mean = 0;%no normalization 
      %here...
      tmp.create_datamat_info.normalize_with_baseline = 1;%yes, normalize
      %to baseline (but changes below, so must override?)
      tmp.create_datamat_info.merge_across_runs = 1;%says to merge across 
      %runs, but we don't do this below anyway...is overridden...??
      tmp.create_datamat_info.single_subject_analysis = 0;%not single subj
      
      
      %%% get subject data based on variability estimate
      %%% and set up appropriate subject specific fields
      tmp.st_sessionFile = [pls_dir subj '_' pattern '_BfMRIsession.mat'];
      %this grabs all relevant session files the fit our naming criteria
      %as defined above.
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % create this subject's datamat
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      tmp.st_datamat = zeros(numel(conditions),numel(st_coords)); %(cond voxel)
      
      % intialize cond specific scan count for populating cond_data
      clear count cond_data block_scan;
      for cond = 1:numel(conditions)
          count{cond} = 0;
      end
      
      switch vartype
        
       case {'stdev_allcorrected'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % within each block express each scan as deviation from block's 
        % temporal mean.Concatenate all these deviation values into one 
        % long condition specific set of scans that were normalized for 
        % block-to-block differences in signal strength. In the end calculate
        % stdev across all normalized cond scans
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % for each condition identify its scans within  runs 
        % and prepare where to put cond specific normalized data
        for cond = 1:numel(conditions)
          tot_num_scans = 0;
          for run = 1:a.session_info.num_runs
            onsets = a.session_info.run(run).blk_onsets{cond}+1;% +1 is because we need matlab indexing convention
            %onsets = a.session_info.run(run).blk_onsets{cond};
            lengths = a.session_info.run(run).blk_length{cond};
            for block = 1:numel(onsets)
              block_scans{cond}{run}{block} = onsets(block)-1+[1:lengths(block)];
              this_length = lengths(block);
              if max(block_scans{cond}{run}{block}>a.session_info.run(run).num_scans)
                disp(['bljak ' subj ' something wrong in block onset lengths']);
                block_scans{cond}{run}{block} = intersect(block_scans{cond}{run}{block},[1:a.session_info.run(run).num_scans]);
                this_length = numel(block_scans{cond}{run}{block});
              end
              tot_num_scans = tot_num_scans + this_length;
            end
          end
          cond_data{cond} = zeros(numel(st_coords),tot_num_scans);
        end
        
        for run = 1:a.session_info.num_runs          
          % load nifti file for this run
          fname = [a.session_info.run(run).data_path '/' a.session_info.run(run).data_files{1}];
          nii = load_nii(fname); %(x y z time)
          img = double(reshape(nii.img,[],size(nii.img,4))); % (voxel time)
          clear nii;
          switch vartype
           case 'stdev_allcorrected'
            % load unsmoothed version of the data and get wm and csf time series
            loc = findstr('smooth_',fname);
            if ~isempty(loc)
              fname1 = [fname(1:(loc-1)) fname((loc+length('smooth_')):end)]
            else fname1=fname; end
            nii = load_nii(fname1);
            img1 = double(reshape(nii.img,[],size(nii.img,4))); 
            wm_ts = reshape(mean(img1(wm_coords,:),1),[],1); 
            csf_ts = reshape(mean(img1(csf_coords,:),1),[],1);
            clear nii img1;
            
            mc_ts = load([workdir subj '/cond' num2str(run) '_mcparams.txt']);
            img = img(st_coords,:);
            temporal_mean = mean(img,2);
            img = residualize([mc_ts wm_ts csf_ts],img')' + repmat(temporal_mean,[1 size(img,2)]);           
          end
          
          for cond = 1:numel(conditions)
            for block = 1:numel(block_scans{cond}{run})
              block_data = img(:,block_scans{cond}{run}{block});% (vox time)
              % normalize block_data to global block mean = 100. 
              block_data = 100*block_data/mean(mean(block_data));
              % temporal mean of this block
              block_mean = mean(block_data,2); % (vox) - this should be 100
              % express scans in this block as  deviations from block_mean
              % and append to cond_data
              good_vox = find(block_mean);              
              for t = 1:size(block_data,2)
                count{cond} = count{cond} + 1;
                cond_data{cond}(good_vox,count{cond}) = block_data(good_vox,t) - block_mean(good_vox);
              end
            end
          end
        end
          
        % now calc stdev across all cond scans
        for cond = 1:numel(conditions)
          tmp.st_datamat(cond,:) = squeeze(std(cond_data{cond},0,2))';
        end
        %all values get saved in approp datamat below, outside of the
        %loops.
        
       case {'mean1_allcorrected'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % within each block express each scan as perc deviation 
        % from block's initial scan. Calculate mean of each block.
        % In the end calculate mean across all blocks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % for each condition identify its scans within  runs 
        % and prepare where to put cond specific normalized data
        for cond = 1:numel(conditions)
          num_blocks  = 0;
          for run = 1:a.session_info.num_runs
            onsets = a.session_info.run(run).blk_onsets{cond}+1; % +1 is because we need matlab indexing convention
            %onsets = a.session_info.run(run).blk_onsets{cond}; 
            lengths = a.session_info.run(run).blk_length{cond};
            for block = 1:numel(onsets)
              block_scans{cond}{run}{block} = onsets(block)-1+[1:lengths(block)];
              this_length = 1;
              if max(block_scans{cond}{run}{block})>a.session_info.run(run).num_scans | min(block_scans{cond}{run}{block})<0
                disp(['bljak ' subj ' something wrong in block onset lengths']);
                block_scans{cond}{run}{block} = intersect(block_scans{cond}{run}{block},[1:a.session_info.run(run).num_scans]);
              end
              if numel(block_scans{cond}{run}{block})<=1, this_length = 0; end
              num_blocks = num_blocks + this_length;
            end
          end
          cond_data{cond} = zeros(numel(st_coords),num_blocks);
        end

        for run = 1:a.session_info.num_runs          
          % load nifti file for this run
          fname = [a.session_info.run(run).data_path '/' a.session_info.run(run).data_files{1}];
          nii = load_nii(fname); %(x y z time)
          img = double(reshape(nii.img,[],size(nii.img,4))); % (voxel time)
          clear nii;
          switch vartype         
           case {'mean1_allcorrected'}
            % load unsmoothed version of the data and get wm and csf time series
            loc = findstr('smooth_',fname);
            if ~isempty(loc)
              fname1 = [fname(1:(loc-1)) fname((loc+length('smooth_')):end)]
            else fname1 = fname;end
            nii = load_nii(fname1);
            img1 = double(reshape(nii.img,[],size(nii.img,4))); 
            wm_ts = reshape(mean(img1(wm_coords,:),1),[],1);
            csf_ts = reshape(mean(img1(csf_coords,:),1),[],1);
            clear nii img1;            
            mc_ts = load([workdir subj '/cond' num2str(run) '_mcparams.txt']);
            img = img(st_coords,:);
            temporal_mean = mean(img,2);
            img = residualize([mc_ts wm_ts csf_ts],img')' + repmat(temporal_mean,[1 size(img,2)]);
          end %switch
          
          for cond = 1:numel(conditions)
            for block = 1:numel(block_scans{cond}{run})
              block_data = img(:,block_scans{cond}{run}{block});% (vox time)
              if size(block_data,2)<=1, continue; end;
              % normalize block_data to global block mean = 100;
              block_data = 100*block_data/mean(mean(block_data));
              % first scan of this block
              block_init = block_data(:,1);
              % in case block_init is 0 for some voxel, then we'll leave that cond_data as 0
              % otherwise express scans in this block as deviation from block_init
              % calculate block mean across scans 2:end and append to cond_data
              count{cond} = count{cond} + 1;
              good_vox = find(block_init);
              scans = 2:size(block_data,2);
              switch vartype
               case {'mean1_allcorrected'}
                cond_data{cond}(good_vox,count{cond}) = mean(100 * (block_data(good_vox,scans) - repmat(block_init(good_vox),[1 numel(scans)])) ./ repmat(block_init(good_vox),[1 numel(scans)]), 2); % here calculate mean as perc change from baseline then I'll end up flattening variance across voxels - not sure if this is desirable or not
              end
            end
          end % for cond
        end % for run
        
        % now calc mean across all block means
        for cond = 1:numel(conditions)
          tmp.st_datamat(cond,:) = squeeze(mean(cond_data{cond},2))';
        end 
          
       otherwise       
        error('batch: wrong varaibility estimate');
      end 
      clear data;
      save([subj '_' pattern '_BfMRIdatamat.mat'],'-struct','tmp','-mat');
 
    end
  end
end

