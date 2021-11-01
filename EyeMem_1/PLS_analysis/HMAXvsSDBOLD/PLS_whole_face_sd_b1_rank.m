clear all
clc
%% Initialize things for model runs...
addpath /Volumes/LNDG/Projects/SocAnx/data/training_testing/DG_checks/TrainingX/Functions_and_scripts/; %must make sure updated PLS software is in your path already. 

%set data paths. In general, should always have result files placed in same
%folder as subject files. Unify????
%TrainTestModelPath = '/Volumes/FB-LIP/LNDG/Projects/SocAnx/data/training_testing/DG_checks/TrainingX/TrainTest_models/';
SubjFilePath = '/Volumes/LNDG/Projects/SocAnx/data/preproc/conditions_n104/condition_53_leosplayground_wholebrain/';%'/Volumes/FB-LIP/LNDG/Projects/SocAnx/data/training_testing/DG_checks/TrainingX/subj_files/';

% we'll also need the PLS toolbox
addpath(genpath('/Volumes/LNDG/Projects/SocAnx/PLSplayground/PLS/'));
% also need Crossval_error
addpath('/Volumes/LNDG/Projects/SocAnx/scripts/CrossVal_Error/')
%set working directory to subj folder and have all models run there...
cd (SubjFilePath)

%% list all IDs, remove outliers if needed, and then randomize.
id_list_full = {'1013tcjp','1016zhqs','1019nngp','1023dwsg','1030wljg','1037twcm','1042kpnq',...
    '1043bxws','1044pgfc','1049dmtf','1052lpcj','1058nrpj','1061nzmt','1065dznl','1066lgcd',...
    '1068gxts','1076mdvr','1089prsj','1091gfnz','1107ltqq','1136xdkb','1141cbxz','1143xjkp',...
    '1149scss','1151thfc','1173pjzz','1181wndz','1184qmqr','1190nnzq','1198ngcd','1201flck',...
    '1207zpjb','1214dvml','1216rnmd','1220krwj','1222nhdb','1226rsmf','1231pwjg','1234mzgm',...
    '1240hxdr','1244zbkq','1254pjtc','1275fppr','1293gdpq','1295vhps','1299kxzw'};
%id_list_full = {'SD_1013tcjp_seven','SD_1016zhqs_seven','SD_1019nngp_seven','SD_1023dwsg_seven','SD_1030wljg_seven','SD_1037twcm_seven','SD_1042kpnq_seven','SD_1043bxws_seven','SD_1044pgfc_seven','SD_1049dmtf_seven','SD_1052lpcj_seven','SD_1058nrpj_seven','SD_1061nzmt_seven','SD_1065dznl_seven','SD_1066lgcd_seven','SD_1068gxts_seven','SD_1076mdvr_seven','SD_1089prsj_seven','SD_1091gfnz_seven','SD_1107ltqq_seven','SD_1136xdkb_seven','SD_1141cbxz_seven','SD_1143xjkp_seven','SD_1149scss_seven','SD_1151thfc_seven','SD_1173pjzz_seven','SD_1181wndz_seven','SD_1184qmqr_seven','SD_1190nnzq_seven','SD_1198ngcd_seven','SD_1201flck_seven','SD_1207zpjb_seven','SD_1214dvml_seven','SD_1216rnmd_seven','SD_1220krwj_seven','SD_1222nhdb_seven','SD_1226rsmf_seven','SD_1231pwjg_seven','SD_1234mzgm_seven','SD_1240hxdr_seven','SD_1244zbkq_seven','SD_1254pjtc_seven','SD_1275fppr_seven','SD_1293gdpq_seven','SD_1295vhps_seven','SD_1299kxzw_seven'};%full ID list
%id_list_full = {'1013tcjp_seven','1016zhqs_seven','1019nngp_seven','1023dwsg_seven','1030wljg_seven','1037twcm_seven','1042kpnq_seven','1043bxws_seven','1044pgfc_seven','1049dmtf_seven','1052lpcj_seven','1058nrpj_seven','1061nzmt_seven','1065dznl_seven','1066lgcd_seven','1068gxts_seven','1076mdvr_seven','1089prsj_seven','1091gfnz_seven','1107ltqq_seven','1136xdkb_seven','1141cbxz_seven','1143xjkp_seven','1149scss_seven','1151thfc_seven','1173pjzz_seven','1181wndz_seven','1184qmqr_seven','1190nnzq_seven','1198ngcd_seven','1201flck_seven','1207zpjb_seven','1214dvml_seven','1216rnmd_seven','1220krwj_seven','1222nhdb_seven','1226rsmf_seven','1231pwjg_seven','1234mzgm_seven','1240hxdr_seven','1244zbkq_seven','1254pjtc_seven','1275fppr_seven','1293gdpq_seven','1295vhps_seven','1299kxzw_seven'};%full ID list

id_list_NoOutliers_idx = find(~ismember(id_list_full, {'NaN'}));%if not outliers to remove, insert 'NaN' in place of outlier IDs.  

rand_seed = rng('shuffle');%logs the random seed, that can later be reused to replicate a given randperm call (and thus and exact id permutation order) in; the line of code above. If want to restore, run: rng(rand_seed) to restore.
id_list_rand = randperm(numel(id_list_NoOutliers_idx))';%randomized index of full list of IDs to then select folds from


% the next step is redundant since id_list_NoOutliers_idx only contains the
% numbers 1:Numberofsubjectswithoutoutlier
% no need to apply a random since random vector and its application are the
% same.
id_list_randidx = id_list_NoOutliers_idx(id_list_rand);
%% optional: load the full model results and see where we're at in comparison
fullmod_res = load('/Users/waschke/Documents/Projects/SocAnx/Full_mod_redvox_n46.mat');

wholebrain_mod = load('/Volumes/LNDG/Projects/SocAnx/data/preproc/conditions_n104/condition_53_leosplayground/condition_53_behav_lsassr_delta_post_b1_n46_wholebrain_BfMRIresult.mat');

%% Paste a column vector of behavioural data of interest from full sample, and convert to cell arrays of strings for use.
delta_lsassr_mat{1} = [-17;-52;-67;4;-33;-13;-14;-43;-36;-61;-16;-16;-25;-66;3;-53;-51;-35;-62;-17;-32;-42;...
    -19;-4;-26;-6;-33;-16;-37;-21;-34;3;-12;-21;-40;-48;-74;-37;-45;-28;-70;-47;-30;-41;-5;-13];%lsassr_delta_post_b1
%to get all values greater than 1 to make all calculations easy. No need for this if data are all positive.
 delta_lsassr_mat{1} = (delta_lsassr_mat{1}*-1)+max(delta_lsassr_mat{1}+1);


delta_lsassr_mat{2} = [-32;-62;-85;-7;-12;-7;2;-53;-21;-54;-17;-6;-18;-49;16;-37;-48;-31;-63;-13;-12;-30;...
    -19;-16;-21;-12;1;-14;-42;-21;-33;15;-25;-18;-35;-24;-91;-41;-38;-29;-75;-49;-24;-48;-4;-25];%lsassr_delta_post_b2
%to get all values greater than 1 to make all calculations easy. No need for this if data are all positive.
 delta_lsassr_mat{2} = (delta_lsassr_mat{2}*-1)+max(delta_lsassr_mat{2}+1);


 baseline_lsassr_mat{1} = [86;78;88;52;80;53;111;88;101;81;57;67;118;117;88;62;78;88;67;49;52;99;73;46;61;93;62;...
     64;82;90;79;93;44;79;57;79;85;78;97;53;104;55;58;72;66;66]; %lsassr_b1



%baseline_lsassr_mat_z{1} = normalize(baseline_lsassr_mat{1});
 baseline_lsassr_mat{2} = [101;88;106;63;59;47;95;98;86;74;58;57;111;100;75;46;75;84;68;45;32;87;73;58;56;99;28;...
    62;87;90;78;81;57;76;52;55;102;82;90;54;109;57;52;79;65;78]; %lsassr_b2

%baseline_lsassr_mat_z{2} = normalize(baseline_lsassr_mat{2});

b = 1; %behavioural variable chosen for indexing below. This will also match delta to baseline scores when needed for regressions below.

behavior_variables = {'lsassr_delta_post_b1','lsassr_delta_post_b2','lsassr_b1','lsassr_b2'}; %for naming behav vars in PLS
behavior_chosen = behavior_variables(b); 



%% Set non-dynamic fields for PLSmodeltxtfilegenerator, etc.

%For relatively fixed options, could instead use inputParser function in Matlab in future instead of
%forcing setting of all manually here...

pls_option = '3';%behavioral PLS
pls_model_types = {'taskPLS','behavPLS'};
if pls_option =='1' %task PLS
    pls_model_chosen = pls_model_types(1);
else if pls_option == '3'
        pls_model_chosen = pls_model_types(2);
    end
end

mean_type = '0';%default for task PLS
cormode = '8';%default
num_perm = '1000';%min 1000!
num_split = '0';%default
num_boot = '1000';%min 1000!
boot_type = 'strat';%default
clim = '95';%default
save_data = '0';%to save stacked datamat. Usually off.
condition_names = {'mean_faces_total_ses-01','mean_1st_face_total_ses-01','mean_1st_face_1sthalf_ses-01',...
    'mean_1st_face_2ndhalf_ses-01','mean_2nd_face_total_ses-01','mean_2nd_face_1sthalf_ses-01',...
    'mean_2nd_face_2ndhalf_ses-01','mean_restingstate_total_ses-01','mean_restingstate_1st_80_ses-01',...
    'mean_restingstate_2nd_80_ses-01','mean_restingstate_1st_40_ses-01','mean_restingstate_2nd_40_ses-01',...
    'mean_restingstate_3rd_40_ses-01','mean_restingstate_4th_40_ses-01','mean_restingstate_1st_20_ses-01',...
    'mean_restingstate_2nd_20_ses-01','mean_restingstate_3rd_20_ses-01','mean_restingstate_4th_20_ses-01',...
    'mean_restingstate_5th_20_ses-01','mean_restingstate_6th_20_ses-01','mean_restingstate_7th_20_ses-01',...
    'mean_restingstate_8th_20_ses-01','mean_fixation_total_ses-01','mean_fixation_1st_ses-01','mean_fixation_2nd_ses-01',...
    'mean_fixation_3rd_ses-01','mean_faces_total_ses-02','mean_1st_face_total_ses-02','mean_1st_face_1sthalf_ses-02',...
    'mean_1st_face_2ndhalf_ses-02','mean_2nd_face_total_ses-02','mean_2nd_face_1sthalf_ses-02','mean_2nd_face_2ndhalf_ses-02',...
    'mean_restingstate_total_ses-02','mean_restingstate_1st_80_ses-02','mean_restingstate_2nd_80_ses-02',...
    'mean_restingstate_1st_40_ses-02','mean_restingstate_2nd_40_ses-02','mean_restingstate_3rd_40_ses-02',...
    'mean_restingstate_4th_40_ses-02','mean_restingstate_1st_20_ses-02','mean_restingstate_2nd_20_ses-02',...
    'mean_restingstate_3rd_20_ses-02','mean_restingstate_4th_20_ses-02','mean_restingstate_5th_20_ses-02',...
    'mean_restingstate_6th_20_ses-02','mean_restingstate_7th_20_ses-02','mean_restingstate_8th_20_ses-02',...
    'mean_fixation_total_ses-02','mean_fixation_1st_ses-02','mean_fixation_2nd_ses-02','mean_fixation_3rd_ses-02',...
    'sd_faces_total_ses-01','sd_1st_face_total_ses-01','sd_1st_face_1sthalf_ses-01','sd_1st_face_2ndhalf_ses-01',...
    'sd_2nd_face_total_ses-01','sd_2nd_face_1sthalf_ses-01','sd_2nd_face_2ndhalf_ses-01','sd_restingstate_total_ses-01',...
    'sd_restingstate_1st_80_ses-01','sd_restingstate_2nd_80_ses-01','sd_restingstate_1st_40_ses-01',...
    'sd_restingstate_2nd_40_ses-01','sd_restingstate_3rd_40_ses-01','sd_restingstate_4th_40_ses-01',...
    'sd_restingstate_1st_20_ses-01','sd_restingstate_2nd_20_ses-01','sd_restingstate_3rd_20_ses-01',...
    'sd_restingstate_4th_20_ses-01','sd_restingstate_5th_20_ses-01','sd_restingstate_6th_20_ses-01',...
    'sd_restingstate_7th_20_ses-01','sd_restingstate_8th_20_ses-01','sd_fixation_total_ses-01',...
    'sd_fixation_1st_ses-01','sd_fixation_2nd_ses-01','sd_fixation_3rd_ses-01','sd_faces_total_ses-02',...
    'sd_1st_face_total_ses-02','sd_1st_face_1sthalf_ses-02','sd_1st_face_2ndhalf_ses-02','sd_2nd_face_total_ses-02',...
    'sd_2nd_face_1sthalf_ses-02','sd_2nd_face_2ndhalf_ses-02','sd_restingstate_total_ses-02',...
    'sd_restingstate_1st_80_ses-02','sd_restingstate_2nd_80_ses-02','sd_restingstate_1st_40_ses-02',...
    'sd_restingstate_2nd_40_ses-02','sd_restingstate_3rd_40_ses-02','sd_restingstate_4th_40_ses-02',...
    'sd_restingstate_1st_20_ses-02','sd_restingstate_2nd_20_ses-02','sd_restingstate_3rd_20_ses-02',...
    'sd_restingstate_4th_20_ses-02','sd_restingstate_5th_20_ses-02','sd_restingstate_6th_20_ses-02',...
    'sd_restingstate_7th_20_ses-02','sd_restingstate_8th_20_ses-02','sd_fixation_total_ses-02',...
    'sd_fixation_1st_ses-02','sd_fixation_2nd_ses-02','sd_fixation_3rd_ses-02'}; %condition names as you have them in sessiondatamat files.

%% #condition_number condition_name
%	1	'mean_faces_total_ses-01'           %	53	'sd_faces_total_ses-01'
%	2	'mean_1st_face_total_ses-01'        %	54	'sd_1st_face_total_ses-01'
%	3	'mean_1st_face_1sthalf_ses-01'      %	55	'sd_1st_face_1sthalf_ses-01'
%	4	'mean_1st_face_2ndhalf_ses-01'      %	56	'sd_1st_face_2ndhalf_ses-01'
%	5	'mean_2nd_face_total_ses-01'        %	57	'sd_2nd_face_total_ses-01'
%	6	'mean_2nd_face_1sthalf_ses-01'      %	58	'sd_2nd_face_1sthalf_ses-01'
%	7	'mean_2nd_face_2ndhalf_ses-01'      %	59	'sd_2nd_face_2ndhalf_ses-01'
%	8	'mean_restingstate_total_ses-01'	%	60	'sd_restingstate_total_ses-01'
%	9	'mean_restingstate_1st_80_ses-01'	%	61	'sd_restingstate_1st_80_ses-01'
%	10	'mean_restingstate_2nd_80_ses-01'	%	62	'sd_restingstate_2nd_80_ses-01'
%	11	'mean_restingstate_1st_40_ses-01'	%	63	'sd_restingstate_1st_40_ses-01'
%	12	'mean_restingstate_2nd_40_ses-01'	%	64	'sd_restingstate_2nd_40_ses-01'
%	13	'mean_restingstate_3rd_40_ses-01'	%	65	'sd_restingstate_3rd_40_ses-01'
%	14	'mean_restingstate_4th_40_ses-01'	%	66	'sd_restingstate_4th_40_ses-01'
%	15	'mean_restingstate_1st_20_ses-01'	%	67	'sd_restingstate_1st_20_ses-01'
%	16	'mean_restingstate_2nd_20_ses-01'	%	68	'sd_restingstate_2nd_20_ses-01'
%	17	'mean_restingstate_3rd_20_ses-01'	%	69	'sd_restingstate_3rd_20_ses-01'
%	18	'mean_restingstate_4th_20_ses-01'	%	70	'sd_restingstate_4th_20_ses-01'
%	19	'mean_restingstate_5th_20_ses-01'	%	71	'sd_restingstate_5th_20_ses-01'
%	20	'mean_restingstate_6th_20_ses-01'	%	72	'sd_restingstate_6th_20_ses-01'
%	21	'mean_restingstate_7th_20_ses-01'	%	73	'sd_restingstate_7th_20_ses-01'
%	22	'mean_restingstate_8th_20_ses-01'	%	74	'sd_restingstate_8th_20_ses-01'
%	23	'mean_fixation_total_ses-01'        %	75	'sd_fixation_total_ses-01'
%	24	'mean_fixation_1st_ses-01'          %	76	'sd_fixation_1st_ses-01'
%	25	'mean_fixation_2nd_ses-01'          %	77	'sd_fixation_2nd_ses-01'
%	26	'mean_fixation_3rd_ses-01'          %	78	'sd_fixation_3rd_ses-01'
%	27	'mean_faces_total_ses-02'           %	79	'sd_faces_total_ses-02'
%	28	'mean_1st_face_total_ses-02'        %	80	'sd_1st_face_total_ses-02'
%	29	'mean_1st_face_1sthalf_ses-02'      %	81	'sd_1st_face_1sthalf_ses-02'
%	30	'mean_1st_face_2ndhalf_ses-02'      %	82	'sd_1st_face_2ndhalf_ses-02'
%	31	'mean_2nd_face_total_ses-02'        %	83	'sd_2nd_face_total_ses-02'
%	32	'mean_2nd_face_1sthalf_ses-02'      %	84	'sd_2nd_face_1sthalf_ses-02'
%	33	'mean_2nd_face_2ndhalf_ses-02'      %	85	'sd_2nd_face_2ndhalf_ses-02'
%	34	'mean_restingstate_total_ses-02'	%	86	'sd_restingstate_total_ses-02'
%	35	'mean_restingstate_1st_80_ses-02'	%	87	'sd_restingstate_1st_80_ses-02'
%	36	'mean_restingstate_2nd_80_ses-02'	%	88	'sd_restingstate_2nd_80_ses-02'
%	37	'mean_restingstate_1st_40_ses-02'	%	89	'sd_restingstate_1st_40_ses-02'
%	38	'mean_restingstate_2nd_40_ses-02'	%	90	'sd_restingstate_2nd_40_ses-02'
%	39	'mean_restingstate_3rd_40_ses-02'	%	91	'sd_restingstate_3rd_40_ses-02'
%	40	'mean_restingstate_4th_40_ses-02'	%	92	'sd_restingstate_4th_40_ses-02'
%	41	'mean_restingstate_1st_20_ses-02'	%	93	'sd_restingstate_1st_20_ses-02'
%	42	'mean_restingstate_2nd_20_ses-02'	%	94	'sd_restingstate_2nd_20_ses-02'
%	43	'mean_restingstate_3rd_20_ses-02'	%	95	'sd_restingstate_3rd_20_ses-02'
%	44	'mean_restingstate_4th_20_ses-02'	%	96	'sd_restingstate_4th_20_ses-02'
%	45	'mean_restingstate_5th_20_ses-02'	%	97	'sd_restingstate_5th_20_ses-02'
%	46	'mean_restingstate_6th_20_ses-02'	%	98	'sd_restingstate_6th_20_ses-02'
%	47	'mean_restingstate_7th_20_ses-02'	%	99	'sd_restingstate_7th_20_ses-02'
%	48	'mean_restingstate_8th_20_ses-02'	%	100	'sd_restingstate_8th_20_ses-02'
%	49	'mean_fixation_total_ses-02'        %	101	'sd_fixation_total_ses-02'
%	50	'mean_fixation_1st_ses-02'          %	102	'sd_fixation_1st_ses-02'
%	51	'mean_fixation_2nd_ses-02'          %	103	'sd_fixation_2nd_ses-02'
%	52	'mean_fixation_3rd_ses-02'          %	104	'sd_fixation_3rd_ses-02'

%%
%                '1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104
 selected_cond = '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0'; %example to select sd_faces_total_ses-01 at this point.

% unfortunately, the PLS code really needs a character array like the one
% above, hm.
% maybe one can change a couple of things in the PLS parser / txtfile
% generator to make this easier. for now let's stick to it

 % come up with a more convenient way of specifying this list of 0s and 1s
% (or only one 1 for now)
% smth like
% conofint = 83;
% conds = zeros(1,numofcondstotal); conds(conofint) = 1;
% condchar=num2str(conds);

  test_cond = ['0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 '...
     '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 '...
     '0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 '...
     '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 '...
     '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0']; %example to select sd_faces_total_ses-01 at this poin
 
 %selected_cond = '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0'; %example to select sd_faces_total_ses-02 at this point.
%selected_cond = '1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0'; %example to select mean_faces_total_ses-01 at this point.
%selected_cond = '0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0'; %example to select mean_faces_total_ses-02 at this point.


%%
%set some fields as a whole that can be called when required.
brain_conds = str2num(selected_cond);%convert selected_cond cell to numbers
% this isn't used again, might be a nice check for the number of variables
% though
sum_brain_conds = sum(brain_conds);

% obv, all of this can be done in less expressions. 
% redundance is key for now in order to be sure whats going on
cond_idx = find(brain_conds);%parse conditions of interest for use in PLS file naming below
conds_chosen = condition_names(cond_idx);

%% create save path for model files/results
%SavePath = ([SubjFilePath,'seed',num2str(rand_seed.Seed),'_',conds_chosen{:},'_',behavior_chosen{:},'/']);

%SavePath = ([TrainTestModelPath,'seed',num2str(rand_seed.Seed),'_',conds_chosen{:},'_',behavior_chosen{:},'/']);
%%%THIS DOESN"T WORK AT ALL. Matlab complains it is all read-only,
%%%directories don't exist, etc. Matlab admin prvileges? How could anything
%%%be written if that's the case?

%% Now loop over creation of training set PLS models, and then compute brain scores in test data.

for i = 1

   
    

    %get fold-based training set indices and sort
    train_fold_idx{i} = 1:46;
    %generate training set ID list for each fold from full id list
    id_train_fold{i} = id_list_full(train_fold_idx{i});
    
    %now generate some variable names that PLS needs to see to generate
    %training data based models...
    txtfilename =    [conds_chosen{:},'_',pls_model_chosen{:},'_',behavior_chosen{:},'_trainall',num2str(i),'_BfMRIanalysis_rank.txt'];
    resultfilename = [conds_chosen{:},'_',pls_model_chosen{:},'_',behavior_chosen{:},'_trainall',num2str(i),'_BfMRIresult_rank.mat'];
    id_list = id_train_fold{i};
    
    %filter behavior variable(s) of interest
    % maybe change naming schemes a bit? make sure no confusion takes place
    behavior_data_train{i} = delta_lsassr_mat{b}(train_fold_idx{i});
    
    % this converts column vector of numbers to cell array of strings to make
    % the PLSmodeltxtfilegenerator work a little more easily :-)
    behavior_data = compose('%d',behavior_data_train{i}); 
    behavior_name = behavior_chosen;%to write to PLS model textfile
    
    %now create PLS model files dynamically
    % maybe tweek the PLStxtfilegenerator a bit in order to make things a
    % bit more intuitive? A structure as input instaed of all the Cells?
    % txt generation in matlab will always stay messy. Let's see if its's
    % worth the trouble.
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,...
        mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,...
        selected_cond,behavior_data,behavior_name)
    
    %run PLS model on training data
    batch_plsgui (txtfilename); %Could also parallelize, as in "Subjectlist_and_parfor_batchpls_subjectfilecreation.m".
    
    %load PLS training set model result you just generated and get brain
    %scores
    load([SubjFilePath,resultfilename],'result');
    
    % example extracting training-set brain scores from 1st LV only 
    % (and adding constant to put all in positive range). 
    % Adapt if more than LV of interest/signif in your model...
%     train_BS{i} = result.usc(:,1) + 100; 
    train_BS{i} = result.usc(:,1) *-1;     


    % thresholding goes here. 
    % not crucial right now, might become useful, though
%     thresh{i} = find(abs(result.boot_result.compare_u)>3.0); %find BS threshold to use as mask
    
%     %compute new training set brain scores for vals BSR +/- 3.0
%     for j=1:numel(id_train_fold{i})
%         load([SubjFilePath,id_train_fold{1,i}{1,j},'_BfMRIsessiondata.mat'],'st_datamat');%could save massive overhead by pre-loading all relevant st_datamat data into matrix beforehand...
%         train_BS_thresh{i}(j,1) = (st_datamat(cond_idx,thresh{i}) * result.u(thresh{i})) + 100;%Adding 100 to make all scores in same positive range as above for train_BS. BE CAREFUL HERE IF MORE THAN ONE CONDITION ARE CHOSEN! Have to adjust this code for that scenario...
%     end
     
    
    %Compute brain scores for test set from training data PLS model results
    % apply weights from trining data to test data.
    % Did it by hand and everything seems to be flawless
   
end


