#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 21:44:05 2016

@author: kloosterman
"""
#reset

import pandas as pd
import matplotlib.pyplot as plt
import hddm
import os
import kabuki
from kabuki.analyze import gelman_rubin

path = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM'
os.chdir(path)

phase = 'test'
# phase = 'study'

data = hddm.load_csv('./EyeMem1_ddm_{}.csv'.format(phase))
#data = hddm.load_csv('./EyeMem_hddm_study.csv')
#data.head(10)

# data = data.rename(columns={'response':'accuracy', 'accuracy': 'response' })

data = data.dropna()
data.loc[data.response==1, 'response'] = 0
data.loc[data.response==2, 'response'] = 1

data.loc[data.age==0, 'age'] = 'YA'
data.loc[data.age==1, 'age'] = 'OA'

data = hddm.utils.flip_errors(data)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions {}'.format(phase))
for i, subj_data in data.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)

plt.savefig('RTdistribution_{}.pdf'.format(phase))

#%% fit HDDM, multiple chains in parallel
# run chains in parallel
def run_model(id):
    import hddm
    data = hddm.load_csv('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM/EyeMem1_ddm_test.csv')
    data = data.dropna()
    # data = data[data.rt > 0.2] # drop too fast RT's
    data.loc[data.response==1, 'response'] = 0
    data.loc[data.response==2, 'response'] = 1
    
    data.loc[data.age==0, 'age'] = 'YA'
    data.loc[data.age==1, 'age'] = 'OA'
    m = hddm.HDDMStimCoding(data, stim_col='stim', split_param='v', drift_criterion=True, bias=True, 
                            depends_on={'v':'age', 'a':'age', 't':'age', 'dc':'age', 'z':'age' }, p_outlier=0.05,) # , include='all'
    m.find_starting_values()
    # m.sample(50, burn=20, dbname='/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/db%i'%id, db='pickle')
    m.sample(500, burn=250, dbname='/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM/db%i'%id, db='pickle')
    return m

from IPython.parallel import Client
v = Client()[:]
jobs = v.map(run_model, range(5)) # 4 is the number of CPUs
models = jobs.get()

a = gelman_rubin(models)
b = pd.DataFrame.from_dict(a, orient='index')
b.to_csv('gelman_rubin_vals.csv')

# Create a new model that has all traces concatenated
# of individual models.
m = kabuki.utils.concat_models(models)

#%% export data
m.save('hddmmodel') # save to file

test = m.gen_stats()
test.to_csv('params_HDDMbias.csv' )

#%% plotting and model fit checks 
m.plot_posteriors_conditions()
plt.savefig('plot_posteriors_conditions.pdf')
m.plot_posteriors(['a', 't', 'v', 'dc', 'z'])
m.plot_posterior_predictive(figsize=(100, 50), ) # bins=1000
# m.plot_posterior_quantiles(samples=3, columns=3, figsize=(100, 50))
plt.show()
plt.savefig('modelfitsHDDMbias.pdf')





#%% OLD
# # Instantiate model object passing it our data (no need to call flip_errors() before passing it).
# # This will tailor an individual hierarchical DDM around your dataset.
# #m = hddm.HDDM(data)
# # find a good starting point which helps with the convergence.
# #m.find_starting_values()
# ## start drawing 7000 samples and discarding 5000 as burn-in
# #m.sample(2000, burn=20)

# m = hddm.HDDM(data, depends_on={'v': 'age', 'a': 'age', 't': 'age'}, bias=False) # , include='all'
# m.find_starting_values()
# #m.sample(10000, burn=1000)
# m.sample(5000, burn=2500)

# #%%

# # TODO find out how to plot young and old in different figures
# m.plot_posterior_predictive(figsize=(20, 20))
# plt.savefig('hddm_subject_fits {}.pdf'.format(phase))


# #%%

# #TODO make loop across ddm para's
# #paralist=["v" "a" "t"]
# paralist=[ 'v', 'a', 't' ]
# paranamelist = [ 'drift-rate', 'boundary separation', 'non-decision time' ]
# # TODO f, ax = plt.subplots(2,2)
# for idx, ipara, in enumerate(paralist):    
#     v_0, v_1 = m.nodes_db.node[['{}(0)'.format(ipara), '{}(1)'.format(ipara)]]    
#     hddm.analyze.plot_posterior_nodes([v_0, v_1])
#     plt.xlabel(ipara)
#     plt.ylabel('Posterior probability')
#     plt.title('Posterior of {} group means, p = {}'.format(paranamelist[idx], (v_0.trace() < v_1.trace()).mean()))
#     plt.savefig('hddm_{}_{}.pdf'.format(phase, paranamelist[idx]))
    

