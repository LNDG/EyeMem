% EM_plot_gazeHMAXvsSDbold
addpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/winsor')

load('SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIresult.mat')

nsub = result.num_subj_lst();
ncond = result.num_conditions;

bsYA = reshape(result.usc(1:nsub(1)*ncond,1), [nsub(1) 3]);
bsOA = reshape(result.usc(nsub(1)*ncond+1:sum(nsub)*ncond,1), [nsub(2) 3]);

do_winsor=0;
if do_winsor
  bsYA = LNDG_winsorize(bsYA);
  bsOA = LNDG_winsorize(bsOA);
end

figure; barweb(mean(bsYA), std(bsYA)/sqrt(length(bsYA)))
figure; barweb(mean(bsOA), std(bsOA)/sqrt(length(bsOA)))
title(['do_winsor' do_winsor])

bsALL = [mean(bsYA); mean(bsOA)];
bsALLsem = [std(bsYA)/sqrt(length(bsYA)); std(bsOA)/sqrt(length(bsOA))];
X = categorical({'YA' 'OA'});
X = reordercats(X, {'YA' 'OA'});
figure; b = barweb(bsALL, bsALLsem);
ylabel('Brain score')
title(sprintf('do_winsor %d', do_winsor))

% ylim([1.5e5 2.5e5])
% legend({char(1:5)})

% RMANOVA with factors Agegroup and binno
agecol = repmat({'YA'}, [length(bsYA) 1])
% t1 = table(agecol, bsYA(:,1),bsYA(:,2),bsYA(:,3),bsYA(:,4),bsYA(:,5), ...
%   'VariableNames',{'agegroup','bin1','bin2','bin3','bin4','bin5'})
% agecol = repmat({'OA'}, [length(bsOA) 1]);
% t2 = table(agecol, bsOA(:,1),bsOA(:,2),bsOA(:,3),bsOA(:,4),bsOA(:,5), ...
%   'VariableNames',{'agegroup','bin1','bin2','bin3','bin4','bin5'});
t1 = table(agecol, bsYA(:,1),bsYA(:,2),bsYA(:,3), ...
  'VariableNames',{'agegroup','bin1','bin2','bin3'})
agecol = repmat({'OA'}, [length(bsOA) 1]);
t2 = table(agecol, bsOA(:,1),bsOA(:,2),bsOA(:,3), ...
  'VariableNames',{'agegroup','bin1','bin2','bin3'});
t=[t1;t2];
% Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bin1-bin3~agegroup ','WithinDesign',Meas)
ranova(rm)

% PREOUT = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/PLS_analysis/HMAXvsSDBOLD'
% writetable(t, fullfile(PREOUT, sprintf('gazeHMAXvsIQRbold_wins%d.csv', do_winsor)))



