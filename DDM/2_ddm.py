import os, itertools
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.anova import AnovaRM
from joblib import Parallel, delayed
from IPython import embed as shell

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set(style='ticks', font='Arial', font_scale=1, rc={
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.linewidth': 0.25,
    'xtick.major.width': 0.25,
    'ytick.major.width': 0.25,
    'ytick.major.width': 0.25,
    'ytick.major.width': 0.25,
    'ytick.major.pad' : 2.0,
    'ytick.minor.pad' : 2.0,
    'xtick.major.pad' : 2.0,
    'xtick.minor.pad' : 2.0,
    'axes.labelpad' : 4.0,
    'axes.titlepad' : 6.0,
    } )
sns.plotting_context()

def params_melt(params, model_settings, flip_b=False):
    try:
        params = params.loc[:,params.columns!='Unnamed: 0']
    except:
        pass
    
    # flip b & z:
    if flip_b:
        params_overall = pd.DataFrame({'subj_idx': np.array(df_emp.groupby(['subj_idx']).first().index.get_level_values('subj_idx')),
                        'c': np.array(df_emp.groupby(['subj_idx']).apply(analyses_tools.behavior, 'c'))})
        b_columns = params.columns[[p[0]=='b' for p in params.columns]]
        for b_column in b_columns:
            params['{}'.format(b_column)] = params[b_column]
            for subj in params['subj_idx'].unique():
                if params_overall.loc[params_overall['subj_idx'] == subj, 'c'].values < 0:
                    params.loc[params['subj_idx'] == subj, '{}'.format(b_column)] = params.loc[params['subj_idx'] == subj, '{}'.format(b_column)] * -1
        z_columns = params.columns[[p[0]=='z' for p in params.columns]]
        for z_column in z_columns:
            params['{}'.format(z_column)] = params[z_column]
            for subj in params['subj_idx'].unique():
                if params_overall.loc[params_overall['subj_idx'] == subj, 'c'].values < 0:
                    params.loc[params['subj_idx'] == subj, '{}'.format(z_column)] = params.loc[params['subj_idx'] == subj, '{}'.format(z_column)] * -1
    
    # melt:
    params = params.melt(id_vars=['subj_idx'])
    for i in range(params.shape[0]):
        variable = "".join(itertools.takewhile(str.isalpha, params.loc[i,'variable']))
        if variable in model_settings['depends_on']: 
            conditions = model_settings['depends_on'][variable]
            if conditions is not None:
                if len(conditions) == 2:
                    params.loc[i,conditions[0]] = int(params.loc[i,'variable'][-3])
                    params.loc[i,conditions[1]] = int(params.loc[i,'variable'][-1])
                elif len(conditions) == 1:
                    params.loc[i,conditions[0]] = int(params.loc[i,'variable'][-1])
        params.loc[i,'variable'] = variable    
    
    return params

# all model elements:
# -------------------

project_dir = '/home/jwdegee/atx/'
exp_name = 'atx'
n_jobs = 24
analysis_step = 0

# dirs:
os.makedirs(os.path.join(project_dir, 'fits'), exist_ok=True)
os.makedirs(os.path.join(project_dir, 'figs'), exist_ok=True) 

# load data:
df = pd.read_csv(os.path.join(project_dir, 'data', '{}.csv'.format(exp_name)))

# compute T_dur:
T_dur = df['rt'].max()+1

# set options:
model_settings = [
    # pyddm:
    {'fit': 'pyddm', 'depends_on': {'a':None, 'u':None, 'v':None, 't':None, 'z':None, 'b':None, 'k':None}, 'start_bias': True, 'drift_bias': True, 'leak': False, 'urgency':False, 'T_dur':T_dur}, #0
    {'fit': 'pyddm', 'depends_on': {'a':['drug'], 'u':None, 'v':['drug', 'difficulty'], 't':['drug'], 'z':['drug'], 'b':['drug'], 'k':None}, 'start_bias': True, 'drift_bias': True, 'leak': False, 'urgency':False, 'T_dur':T_dur}, #1
    {'fit': 'pyddm', 'depends_on': {'a':['drug'], 'u':None, 'v':['drug', 'difficulty'], 't':['drug'], 'z':['drug'], 'b':['drug'], 'k':None}, 'start_bias': True, 'drift_bias': True, 'leak': False, 'urgency':True, 'T_dur':T_dur}, #1
    {'fit': 'pyddm', 'depends_on': {'a':['drug'], 'u':['drug'], 'v':['drug', 'difficulty'], 't':['drug'], 'z':['drug'], 'b':['drug'], 'k':None}, 'start_bias': True, 'drift_bias': True, 'leak': False, 'urgency':True, 'T_dur':T_dur}, #1
    ]

versions = [1,2,3]
for version in versions:

    df_emp = df.copy()
    df_emp['response2'] = np.array(df_emp['response'])
    df_emp.loc[(df_emp['simon']==1)&(df_emp['response2']==0), 'response2'] = 3
    df_emp.loc[(df_emp['simon']==1)&(df_emp['response2']==1), 'response2'] = 0
    df_emp.loc[(df_emp['simon']==1)&(df_emp['response2']==3), 'response2'] = 1
    df_emp['response'] = np.array(df_emp['response2'])
    
    # fit model:
    if analysis_step == 0:
        
        if model_settings[version]['fit'] == 'pyddm':
            from accumodels import pyddm_tools
            df_emp.loc[df_emp['stimulus']==0, 'stimulus'] = -1 # pyddm expects -1 and 1 as stimuli identifiers
            res = Parallel(n_jobs=n_jobs, verbose=1, backend='loky')(delayed(pyddm_tools.fit_model)
                            (data, model_settings[version], subj_idx) for subj_idx, data in df_emp.groupby(['subj_idx']))
            
            
            
            params = pd.concat(res).reset_index(drop=True)
        params.to_csv(os.path.join(project_dir, 'fits', '{}_{}.csv'.format(exp_name, version)))

    elif analysis_step == 1:
        
        from accumodels import pyddm_tools, plot_tools

        # simulate data:
        params = pd.read_csv(os.path.join(project_dir, 'fits', '{}_{}.csv'.format(exp_name, version)))
        
        df_emp.loc[df_emp['stimulus']==0, 'stimulus'] = -1 # pyddm expects -1 and 1 as stimuli identifiers
        df_sim = pd.concat(Parallel(n_jobs=n_jobs, verbose=1, backend='loky')(delayed(pyddm_tools.simulate_data)(data, params, model_settings[version], subj_idx, 100000)
                            for subj_idx, data in df_emp.groupby(['subj_idx']))).reset_index()
        df_sim.to_csv(os.path.join(project_dir, 'fits', '{}_{}_df_sim.csv'.format(exp_name, version)))
        
        # model fit:
        if not 'drug' in df_sim.columns:
            df_emp['drug'] = 1
            df_sim['drug'] = 1
        for b, d in df_emp.groupby(['drug']):
            fig = plot_tools.summary_plot_group(df_emp.loc[df_emp['drug']==b,:], df_sim.loc[df_sim['drug']==b,:])
            fig.savefig(os.path.join(project_dir, 'figs', '{}_{}_model_fit_{}.pdf'.format(exp_name, version, b)))

# analyses:
if analysis_step == 2:
    
    bics = []
    resids = []
    sdt_emps = []
    sdt_sims = []
    df_sims = []
    for version in [1,2,3]:
        
        # load:
        params = pd.read_csv(os.path.join(project_dir, 'fits', '{}_{}.csv'.format(exp_name, version)))
        print('BIC version {} = {}'.format(version, params['bic'].mean()))
        params.loc[:, [p[0]=='a' for p in params.columns]] = params.loc[:, [p[0]=='a' for p in params.columns]] * 2 
        params = params_melt(params, model_settings[version])
        # df_sim = pd.read_csv(os.path.join(project_dir, 'fits', '{}_{}_df_sim.csv'.format(exp_name, version)))
        # df_sim.loc[df_sim['stimulus']==-1, 'stimulus'] = 0 # we now expects 0 and 1 as stimuli identifiers
        # df_sims.append(df_sim)
        
        from statsmodels.stats.anova import AnovaRM
        aovrm = AnovaRM(params.loc[(params['variable']=='v'),:], 'value', 'subj_idx', within=['drug', 'difficulty'])
        res = aovrm.fit()
        print(res)

        params.loc[(params['variable']=='v')&(params['difficulty']==0), 'variable'] = 'v_easy'
        params.loc[(params['variable']=='v')&(params['difficulty']==0), 'variable'] = 'v_easy'
        params.loc[(params['variable']=='v')&(params['difficulty']==1), 'variable'] = 'v_hard'
        params.loc[(params['variable']=='v')&(params['difficulty']==1), 'variable'] = 'v_hard'

        # plot bars:
        fig = plt.figure(figsize=(4,2))
        ax = fig.add_subplot(111)
        sns.barplot(x='variable', y='value', hue='drug', data=params, ax=ax)
        plt.tight_layout()
        fig.savefig(os.path.join(project_dir, 'figs', '{}_{}_bars.pdf'.format(exp_name, version)))