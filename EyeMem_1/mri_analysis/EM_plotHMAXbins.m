%% EM_plotHMAXbins TODO do in taskPLS folder!

PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/taskPLS/gaze-specific/';
cd(PREIN)
% load SDbold_vs_HMAX_youngold_44_41_BfMRIresult_10000perms.mat
load SDbold_vs_HMAX_youngold_44_41_BfMRIresult.mat
% load id_list.mat
% disp 'TODO make this robust'
% id_list = id_list{:};
% idx = strcmp(id_list,'sub-68'); % drop
% id_list(idx) = [];

% dirlistYA = dir('young/sub-*.mat');
% dirlistOA = dir('old/sub-*.mat');
ages = {'young', 'old'};

hmax_meanperbin=[];  subjlist={};
for iage = 1:2
  cd(ages{iage})
  dirlist = dir('sub-*.mat');
  for i = 1:length(dirlist)
    disp(dirlist(i).name);
    tmp = load(dirlist(i).name);
    hmax_meanperbin{iage}(i,:) = tmp.hmax_meanperbin;
    tok = tokenize(dirlist(i).name, '_');
    subjlist{iage}{i} = tok{1};
  end
  cd ..
end

% dirlist = dir('sub-*.mat');
% hmax_meanperbin=[];  subjlist={};
% for i = 1:length(dirlist)
%   disp(dirlist(i).name);
%   tmp = load(dirlist(i).name);
%   hmax_meanperbin(i,:) = tmp.hmax_meanperbin;
%   tok = tokenize(dirlist(i).name, '_');
%   subjlist{i,1} = tok{1};
% end
% 

nsub = sum(result.num_subj_lst);
vsc = reshape(result.vsc(:,1), nsub, [] );
nbins = size(vsc,2);
vsc_YA = vsc(1:result.num_subj_lst(1),:);
vsc_OA = vsc(result.num_subj_lst(1)+1:end,:);
figure; subplot(2,1,1)
plot(mean(vsc_YA)); hold on; plot(mean(vsc_OA));  title('vsc (Designscores)')

nsub = sum(result.num_subj_lst);
% usc = reshape(result.usc(:,1), nsub, [] );
% usc = reshape(result.boot_result.usc2(:,1), nsub, [] );
% usc = reshape(result.boot_result.usc2(:,1), 5, [] )';
idx = 1:result.num_subj_lst(1)*result.num_conditions;
usc_YA = reshape(result.usc(idx,1), [], 5 );
usc_OA = reshape(result.usc(idx(end)+1:end,1), [], 5 );

figure; subplot(2,1,1)
% plot(mean(vsc_YA)); hold on; plot(mean(vsc_OA));  title('vsc (Designscores)')
plot(result.v(1:5,1)); hold on;  plot(result.v(6:10,1)); title('v (design scores)')
subplot(2,2,3);
plot(mean(usc_YA)); hold on; 
subplot(2,2,4);
plot(mean(usc_OA));  
title('usc (Brainscores)')

%%
disp 'rmcorr hmax vs brain scores'
usc_YA_demean = usc_YA - mean(usc_YA,2);
usc_OA_demean = usc_OA - mean(usc_OA,2);
hmax_YA_demean = hmax_meanperbin{1} - mean(hmax_meanperbin{1},2);
hmax_OA_demean = hmax_meanperbin{2} - mean(hmax_meanperbin{2},2);
[r_YA, p_YA] = corr(usc_YA_demean(:), hmax_YA_demean(:), 'type', 'Spearman')
[r_OA, p_OA] = corr(usc_OA_demean(:), hmax_OA_demean(:), 'type', 'Spearman')
f=figure; 
subplot(1,2,1); scatter(usc_YA_demean(:), hmax_YA_demean(:)); axis square; box on; lsline
subplot(1,2,2); scatter(usc_OA_demean(:), hmax_OA_demean(:)); axis square; box on; lsline

%%
mri = ft_read_mri('SDbold_vs_HMAX_youngold_44_41_BfMRIbsr_lv1.hdr');
% ft_write_
