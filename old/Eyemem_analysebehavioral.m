basepath = '/Users/kloosterman/Dropbox/PROJECTS/Memory_Eyetracking/Eye_only_pilot/data/preproc';
cd(basepath)

youngfemale = [6 2 9 12 24 12];
youngmale = [10 8 14 23 25];

oldfemale = [1 4 20 13 17 18 16];
oldmale = [5 3 7 15 11 21 19 22];

young = [2     6     8     9    10    12    12    14    23    24    25];
old = [ 1     3     4     5     7    11    13    15    16    17    18    19    20    21    22];

SUBJ = {  'S13'	'S14'	'S15'	'S16'	'S20'	'S21'	'S222'	'S23'	'S24'	'S25'	};
SUBJno = [  13	14	15	16	20	21	22	23	24	25	];

stimcategories = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'};

behav=[];
for isub = 1:length(SUBJ)
    cd(SUBJ{isub})
    load([SUBJ{isub} '_test_stim_data.mat'] )
    for icat=1:5
        
        trialinfo = data.trialinfo(data.trialinfo(:,2) == icat, :);
        
        behav.accuracy(isub,icat) = length(find(trialinfo(:,5) == trialinfo(:,6))) / length(trialinfo);  % NEW=1, OLD=2
        
        
        behav.Hitrate(isub,icat) = length(find( trialinfo(:,5) == 2 & trialinfo(:,6) == 2 )) ... 
        /  length(find( trialinfo(:,5) == 2 ));
 
        behav.FArate(isub,icat) = length(find( trialinfo(:,5) == 1 & trialinfo(:,6) == 2 )) ... 
        /  length(find( trialinfo(:,5) == 2 )); 
        
        behav.RT(isub,icat) = mean(trialinfo(:,7)) / 1000;
    end
    cd ..
end

behav.Hitrate(behav.Hitrate == 1) = 0.95;
behav.FArate(behav.FArate == 0) = 0.05;

behav.dprime = norminv(behav.Hitrate) - norminv(behav.FArate);
behav.criterion = -0.5 * (norminv(behav.Hitrate) + norminv(behav.FArate));

%% plot
close all
accuracy = [behav.accuracy];
% accuracy = [behav.criterion];
acc_old = mean(accuracy(ismember(SUBJno, old),:));
acc_young = mean(accuracy(ismember(SUBJno, young),:));

% RT = [behav.dprime];
RT = [behav.RT];
RT_old = mean(RT(ismember(SUBJno, old),:));
RT_young = mean(RT(ismember(SUBJno, young),:));

set(0, 'defaultaxesfontsize', 14)
figure; subplot(2,1,1);
bar([acc_young' acc_old' ] )
set(gca, 'XTicklabel', stimcategories, 'Tickdir', 'out')
legend({'Young', 'Old'})
title('% correct')



subplot(2,1,2);
bar([RT_young' RT_old' ] )
set(gca, 'XTicklabel', stimcategories, 'Tickdir', 'out')
legend({'Young', 'Old'})
title('dprime')



