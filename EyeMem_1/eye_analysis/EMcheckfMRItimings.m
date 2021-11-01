cd /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye/OA
close all
list = dir('eye*.mat')
piconsets =[];
for i=1:length(list)
  disp(list(i).name)
  load(list(i).name)
  if length(data.trialinfo) ~= 150
    disp('not 150 trials')
    continue
  end
%   piconsets{i} = round(data.trialinfo([1,31,61,91,121],9)) + 5; %12 subtracted, bu trigger only sent at 4 ?!
  piconsets{i} = data.trialinfo([1,31,61,91,121],9); %12 subtracted, bu trigger only sent at 4 ?!
  
  
end
piconsets = [piconsets{:}]

unique(piconsets(:))
figure; histogram(piconsets(:))
title('range should be between 6 and 11')
