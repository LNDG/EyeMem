function EM_decodefMRI_plot()

agedirs = {'YA' 'OA'};

d=cell(2,1);
for iage = 1:2
  cd(sprintf('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/LOTO/%s/1_decodingMSE', agedirs{iage}))
  subjlist = dir('*.mat');
  
  for isub = 1:length(subjlist)
    disp(subjlist(isub).name)
    load(subjlist(isub).name)
    if sum(accuracy) > 0
      d{iage}(end+1,:) = accuracy;
    end
  end
end
%%
figure;
for iage = 1:2
  % plot([-1:6], d)
  hold on
  plot([-1:6], mean(d{iage}), 'Linewidth', 4)
end
legend(agedirs)