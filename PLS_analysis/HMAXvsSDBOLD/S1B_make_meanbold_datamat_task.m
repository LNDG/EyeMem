function S2_make_meanbold_datamat_task()

% Create the mean datamats based on the info mats.

% 180228 | adapted from STSW rest to task

pn.root         = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/mri/task/';
pn.plstoolbox   = [pn.root, 'analyses/B4_PLS_preproc2/T_tools/'];  addpath(genpath(pn.plstoolbox));
pn.PLSfiles     = [pn.root, 'analyses/B4_PLS_preproc2/B_data/SD_STSWD_task_v2/'];

cd(pn.PLSfiles); % files will be output in the CD directory

% N = 43 YA + 53 OA;
IDs = {'1117';'1118';'1120';'1124';'1125';'1126';'1131';'1132';'1135';'1136';...
    '1151';'1160';'1164';'1167';'1169';'1172';'1173';'1178';'1182';'1214';'1215';...
    '1216';'1219';'1223';'1227';'1228';'1233';'1234';'1237';'1239';'1240';'1243';...
    '1245';'1247';'1250';'1252';'1257';'1261';'1265';'1266';'1268';'1270';'1276';'1281';...
    '2104';'2107';'2108';'2112';'2118';'2120';'2121';'2123';'2125';'2129';'2130';...
    '2131';'2132';'2133';'2134';'2135';'2139';'2140';'2145';'2147';'2149';'2157';...
    '2160';'2201';'2202';'2203';'2205';'2206';'2209';'2210';'2211';'2213';'2214';...
    '2215';'2216';'2217';'2219';'2222';'2224';'2226';'2227';'2236';'2237';'2238';...
    '2241';'2244';'2246';'2248';'2250';'2251';'2252';'2258';'2261'};

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);
    BATCHPATH = [pn.PLSfiles, IDs{indID}, '_task_PLS_info.mat'];
    fprintf('Creating PLS_data_mat for subject %s \n', IDs{indID});
    batch_plsgui(BATCHPATH);
    fprintf('Finished creating PLS_data_mat for subject %s \n', IDs{indID});
end

end