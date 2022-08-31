function EM_rmdirs()

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; %yesno or 2afc
  backend = 'slurm';
%       backend = 'torque';
%     backend = 'local';
  compile = 'no';
end

dir2rem = fullfile(basepath, 'variability2/1TRspertrial')
cd(dir2rem)
dirlist= dir('sub-*');

for idir = 1:length(dirlist)
  if contains(dirlist(idir).name, 'temp')
    disp(dirlist(idir).name)
%     rmdir(dirlist(idir).name)
    disp('deleted')
  else
    disp(dirlist(idir).name)
    cd(dirlist(idir).name)
%     rmdir('spm', 's')
    disp('deleted')
    cd ..
  end
end
