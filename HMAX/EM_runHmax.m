function [c2,c1,bestBands,bestLocations,s2,s1] = EM_runHmax(cfg)
%    run Hmax on eyemem pics
%
% [c2,c1,bestBands,bestLocations,s2,s1] = example(exampleImages,saveFolder)
%
% An example code showing how HMAX is initalized and used. This function
% will call all the relevant functions to generate C2 activations for the
% patches and images provided.
%
% args:
%
%     exampleImages: a cell array. Each cell should contain the full
%                    path to an image.
%
%     saveFolder: a string. Directory to save the output in.
%
% returns:
%
%     c2,s2,c1,s1: see C1.m and C2.m
%
%     filters: the gabor filters used for computing S1 responses

pictures_path = cfg.pictures_path;
saveFolder = cfg.saveFolder;
piclist = cfg.piclist;
condition = cfg.condition;

%% If no arguments are provided, these are the default variables.
if isempty(pictures_path) && isempty(saveFolder)
%   saveFolder = './output/';
  %   load('exampleImages.mat');
  pictures_path = '/Volumes/LNDG/Projects/EyeMem/eyemem/study_information/D_paradigm/fMRI_exp_final/stimuli_640x480/fractals';
%   saveFolder = fullfile(pictures_path, 'hmax')
  saveFolder = fullfile(pictures_path);
end

%% Preprocess the images.
% Creates a cell array with each cell containing a grayscaled
% representation of one image. Data type should be double, not uint8.
cd(pictures_path)

if isempty(piclist)
  disp('Running for all pics in folder')
  delete('.*')
  piclist = dir('*.bmp');
  if isempty(piclist);
    piclist = dir('*.png');
  end
  npics = length(piclist);
  
  for ipic = 1:npics
    picdat{ipic} = double(rgb2gray(imread(piclist(ipic).name)));
  end
else
  disp('using provided piclist')
end

npics = length(piclist);
picno = nan(npics,1);
for ipic = 1:npics
  picdat{ipic} = double(rgb2gray(imread(piclist{ipic})));
  if strfind(piclist{ipic}, 'image')
    [~,a] = fileparts(piclist{ipic}(6:end));
%     piclist{ipic} = piclist{ipic}(6:end);
  else
    [~,a] = fileparts(piclist{ipic});
  end  
  picno(ipic,1) = str2num(a);
end

%% Initialize S1 gabor filters and C1 parameters
fprintf('initializing S1 gabor filters\n');
orientations = [90 -45 0 45]; % 4 orientations for gabor filters
RFsizes      = 7:2:37;        % receptive field sizes
div          = 4:-.05:3.25;    % tuning parameters for the filters' "tightness"
[filterSizes,filters,c1OL,~] = initGabor(orientations,RFsizes,div);

fprintf('initializing C1 parameters\n')
c1Scale = 1:2:18; % defining 8 scale bands
c1Space = 8:2:22; % defining spatial pooling range for each scale band

%% Load the universal patch set.
fprintf('Loading the universal patch set\n')
load('universal_patch_set.mat','patches','patchSizes');

nPatchSizes = size(patchSizes,2);


%% For each patch calculate responses
fprintf('calculating unit responses\n');

hmaxout =[];
hmaxout.c1median = nan(npics, 8);
% hmaxout.c1median = nan(npics, 8);
for ipic = 1:npics % cf https://doi.org/10.1101/249029

%   [c2,c1,bestBands,bestLocations,s2,s1] = extractC2forCell...
%     (filters,filterSizes,c1Space,c1Scale,c1OL,patches,   picdat(ipic),  nPatchSizes,patchSizes(1:3,:));
  [c2,c1] = extractC2forCell...
    (filters,filterSizes,c1Space,c1Scale,c1OL,patches,   picdat(ipic),  nPatchSizes,patchSizes(1:3,:));

  c1_ori4 = cellfun(@(x) x(:,:,4), c1{1}, 'uni', 0); % take ori 4
  hmaxout.c1median(ipic,:) = cell2mat(cellfun(@(x) median(x(:)), c1_ori4, 'uni', 0)); % take median
  hmaxout.c2median(ipic,:)  = cellfun(@median, c2, 'uni', 0); % take median
%   hmaxout.picname = piclist(ipic).name;
  clear c2 c1 
end
hmaxout.picname = piclist;
hmaxout.picno = picno;


%% Save the output
% Note that c1, s2, and s1 do not get saved.
save( fullfile(saveFolder, sprintf('hmax_%s.mat', condition)), 'hmaxout' );

