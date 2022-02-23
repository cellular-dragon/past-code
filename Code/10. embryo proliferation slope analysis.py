# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 16:52:37 2022

@author: dell
"""


import pandas as pd
import numpy as np
import scipy
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
#input the division timing matrix of ctr
control = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\Cell_division_timing_matrix_ctr.txt", sep="\t", index_col=0)

#input the division timing matrix of RNAi
experiment = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\Cell_division_timing_matrix_ctr.txt", sep="\t", index_col=0)

#input the the path and name of result file 
out_put=(r"C:\Users\dell\Desktop\\RNAi_slope.txt")
control_mean = np.mean(control,axis=1)
control_cell_list = list(control_mean.index.values)
ctr_lin = []
experiment_lin=[]

for i in control.columns:
    control_tmp=control[i]
    peremetre_control = stats.linregress(control_mean,control_tmp)
    
    if type(ctr_lin) == list:
        ctr_lin = pd.DataFrame(list(peremetre_control))                                  
    else:
        ctr_lin = pd.concat([ctr_lin, pd.DataFrame(list(peremetre_control))], axis=1)

ctr_lin.columns=control.columns    
ctr_lin.index=['slope', 'intercept', 'r_value', 'p_value', 'std_err']        
a=0        
for j in experiment.columns:
    experiment_tmp=experiment[j]
    
    experiment_tmp_cell_list = list(experiment_tmp.dropna().index.values)
    
    if experiment_tmp_cell_list == control_cell_list:
        peremetre_experiment = stats.linregress(control_mean,experiment_tmp)
        a+=1
        print(a)
    else: 
        
        cell_list_cross = list(set(experiment_tmp_cell_list) & set(control_cell_list))
        a+=1
        print(a)
        
        
        experiment_tmp=experiment_tmp.loc[cell_list_cross]
        control_mean_tmp = control_mean.loc[cell_list_cross]
        peremetre_experiment = stats.linregress(control_mean_tmp,experiment_tmp)
        
    if type(experiment_lin) == list:
        experiment_lin = pd.DataFrame(list(peremetre_experiment))                                  
    else:
        experiment_lin = pd.concat([experiment_lin, pd.DataFrame(list(peremetre_experiment))], axis=1)

experiment_lin.columns=experiment.columns    
experiment_lin.index=['slope', 'intercept', 'r_value', 'p_value', 'std_err']
   
ctr_slope_mean=np.mean(ctr_lin.loc['slope'])
ctr_slope_std=np.std(ctr_lin.loc['slope'])

slope_Z_score = (experiment_lin.loc['slope'] - ctr_slope_mean)/ctr_slope_std
p_value_RNAi_embryo=pd.DataFrame(1-scipy.stats.norm.cdf(abs(slope_Z_score)),index=slope_Z_score.index,columns=['slope P-value'])

q_value_bh=multipletests(list(p_value_RNAi_embryo['slope P-value'].values), method='fdr_bh')[1]
q_value_bh=pd.DataFrame(q_value_bh,index=p_value_RNAi_embryo.index,columns=['slope q-value bh'])

sign_matrix = np.sign(slope_Z_score.values)
sign_matrix_pd = pd.DataFrame(sign_matrix)
sign_matrix_pd = sign_matrix_pd.replace(0.0, 1.0)
q_value_bh = pd.DataFrame(sign_matrix_pd.values*(q_value_bh.values))

q_value_bh.index=p_value_RNAi_embryo.index
q_value_bh.columns=['slope q-value bh']
q_value_bh['z-score']=slope_Z_score
q_value_bh=pd.concat([q_value_bh,pd.DataFrame(experiment_lin.loc['slope'])],axis=1)

q_value_binary=[]
for q in q_value_bh.index:
    q_value_temp=q_value_bh.loc[q]['slope q-value bh']
    if q_value_temp <0.01 and q_value_temp>-0.01:
        q_value_binary.append(1)
    else:
        q_value_binary.append(0)    

q_value_bh['q-value binary']=q_value_binary

q_value_bh.to_csv(out_put, sep="\t")




