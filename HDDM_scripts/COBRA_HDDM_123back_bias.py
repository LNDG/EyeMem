#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 21:44:05 2016
conda activate hddm_env
ipcluster start
spyder
@author: kloosterman
"""
#reset
# %matplotlib inline

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import hddm
import os
import kabuki
import platform
from kabuki.analyze import gelman_rubin

if platform.system() == 'Darwin':
    path = '/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data'
elif platform.system() == 'Linux':
    path = '/home/mpib/kloosterman/MATLAB/COBRA/123back_bias_novelvsfam/data'

os.chdir(path)

# csvfile = '/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata.csv'
csvfile = os.path.join(path, 'COBRA_DDMdata_drop_lowdprime.csv')
# csvpath = '/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias/COBRA_HDDM_May26_2021_response-button_press_incsubj_stimcoding.csv'
data = hddm.load_csv(csvfile)
#data.head(10)

#%% plot RT distributions
data = data.dropna()
data = data[data.rt > 0.2] # drop too fast RT's
data = hddm.utils.flip_errors(data)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby('subj_idx'):
    subj_data.rt.hist(bins=25, histtype='step', ax=ax)

plt.savefig('RTdistribution.pdf')
plt.show()

#%% fit HDDM, multiple chains in parallel
# run chains in parallel inc. bias models
def run_biasmodel(id):
    # Why are you loading and cleaning the data again here? seems repetitive, equivalent to what you did above outside the function
    import hddm
    data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')
    data = data.dropna()
    data = data[data.rt > 0.2] # drop too fast RT's
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
                            depends_on={'v':'stim', 'a':'stim', 't':'stim', 'dc':'stim', 'z':'stim' }, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    m.sample(5000, burn=2500, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/db_bias%i'%id, db='pickle')
    return m

def run_biasmodel_early(id):
    import hddm
    import platform
    if platform.system() == 'Darwin':
        data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')
    else:
        data = hddm.load_csv('/home/mpib/kloosterman/MATLAB/COBRA/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')

    data = data.dropna()
    data = data[data.late==0]
    data = data[data.rt > 0.2] # drop too fast RT's
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
                            depends_on={'v':'stim', 'a':'stim', 't':'stim', 'dc':'stim', 'z':'stim' }, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    m.sample(1000, burn=250, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/db_bias%i'%id, db='pickle')
    # m.sample(1000, burn=250, dbname= os.path.join(path, 'db_bias%i'%id), db='pickle')
    return m

def run_biasmodel_late(id):
    import hddm
    import platform
    if platform.system() == 'Darwin':
        data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')
    else:
        data = hddm.load_csv('/home/mpib/kloosterman/MATLAB/COBRA/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')

    data = data.dropna()
    data = data[data.early==0]
    data = data[data.rt > 0.2] # drop too fast RT's
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
                            depends_on={'v':'stim', 'a':'stim', 't':'stim', 'dc':'stim', 'z':'stim' }, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    m.sample(1000, burn=250, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/db_bias%i'%id, db='pickle')
    # m.sample(1000, burn=250, dbname= os.path.join(path, 'db_bias%i'%id), db='pickle')
    return m

def run_biasmodel_middle(id):
    import hddm
    import platform
    if platform.system() == 'Darwin':
        data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')
    else:
        data = hddm.load_csv('/home/mpib/kloosterman/MATLAB/COBRA/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')

    data = data.dropna()
    data = data[(data.early==1) & (data.late==1)]
    data = data[data.rt > 0.2] # drop too fast RT's
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
                            depends_on={'v':'stim', 'a':'stim', 't':'stim', 'dc':'stim', 'z':'stim' }, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    m.sample(1000, burn=250, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/db_bias%i'%id, db='pickle')
    # m.sample(1000, burn=250, dbname= os.path.join(path, 'db_bias%i'%id), db='pickle')
    return m

def run_standardmodel(id):
    import hddm
    data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/COBRA_DDMdata_drop_lowdprime.csv')
    data = data.dropna()
    data = data[data.rt > 0.2] # drop too fast RT's
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=False, bias=False,
                            depends_on={'v':'stim', 'a':'stim', 't':'stim'}, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    m.sample(5000, burn=2500, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias_novelvsfam/data/db_standard%i'%id, db='pickle')
    return m

from IPython.parallel import Client
v = Client()[:]

# run model
jobs = v.map(run_biasmodel, range(15)) # 4 is the number of CPUs

models = jobs.get()

a = gelman_rubin(models)
b = pd.DataFrame.from_dict(a, orient='index')
b.to_csv('gelman_rubin_vals_drop_lowdprime.csv')

# Create a new model that has all traces concatenated
# of individual models.
m = kabuki.utils.concat_models(models)

#%% export data
m.save('hddmmodel_biasmodel_drop_lowdprime') # save to file

test = m.gen_stats()
test.to_csv('params__biasmodel_drop_lowdprime.csv' )

#%% plotting and model fit checks
# a = m.plot_posteriors_conditions()
# plt.savefig('plot_posteriors_conditions.pdf')
# m.plot_posteriors(['a', 't', 'v', 'dc', 'z'])

# m.plot_posterior_predictive(figsize=(100, 50), ) # bins=1000

m.plot_posterior_predictive(figsize=(100, 50), value_range= np.linspace(-1.5, 1.5, 10)) # bins=1000

# # m.plot_posterior_quantiles(samples=3, columns=3, figsize=(100, 50))
# plt.show()
# plt.savefig('modelfitsHDDMbiasmodel_late.pdf')


# TODO LATER fix statsmodels
# from kabuki.analyze import check_geweke
# import statsmodels
# print check_geweke(m)


# standard DDM
# subj_params = []
# for subj_idx, subj_data in data.groupby('subj_idx'):
#     m_subj = hddm.HDDM(subj_data, depends_on={'v':'stim', 'a':'stim', 't':'stim' }, p_outlier=0.05,)
#     # m_subj = hddm.HDDM(subj_data, depends_on={'v':'stim' }, p_outlier=0.05,)
#     subj_params.append(m_subj.optimize(method='chisquare'))
# params = pd.DataFrame(subj_params)
# # params.to_csv('params_cobra_onlyvfree.csv' )
# params.to_csv('params_cobra_optimize.csv' )

# m = hddm.HDDM(data, depends_on={'v':'stim', 'a':'stim', 't':'stim' }, bias=False,  p_outlier=0.05,) # , include='all'
# m.find_starting_values()
# m.sample(10000, burn=1000)
# m.sample(5000, burn=2500, thin=2) # Urai 2019 eLife settings
# m.sample(5000, burn=2500, dbname='traces.db', db='pickle') # Frank 2021 paper
# # m.sample(5, burn=1, dbname='traces.db', db='pickle') # test settings
# m.save('hddmmodel') # save to file
# to load: m = hddm.load('hddmmodel')


# # #TODO plot results
# #paralist=["v" "a" "t"]
# paralist=[ 'v', 'a', 't' ]
# paranamelist = [ 'drift-rate', 'boundary separation', 'non-decision time' ]
# phase = ''
# # TODO f, ax = plt.subplots(2,2)
# for idx, ipara, in enumerate(paralist):
#     v_0, v_1 = m.nodes_db.node[['{}(2-back)'.format(ipara), '{}(3-back)'.format(ipara)]]
#     hddm.analyze.plot_posterior_nodes([v_0, v_1])
#     plt.xlabel(ipara)
#     plt.ylabel('Posterior probability')
#     plt.title('Posterior of {} group means, p = {}'.format(paranamelist[idx], (v_0.trace() < v_1.trace()).mean()))
#     plt.savefig('hddm_{}_{}.pdf'.format(phase, paranamelist[idx]))
#
