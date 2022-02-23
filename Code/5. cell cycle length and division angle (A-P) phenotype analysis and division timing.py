# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 20:14:41 2022

@author: dell
"""

# -*- coding: utf-8 -*-
import scipy
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import os

#######input the path and name of control file of cell cycle length or division timing or division angle (A-P)
control = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\Cell_cycle_length_matrix_ctr.txt", sep="\t", index_col=0)                 #in 输入文件 control
#######input the path and name of RNAi file of cell cycle length or division timing or division angle (A-P)
experiment = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\Cell_cycle_length_matrix_RNAi.txt", sep="\t", index_col=0)                 #in 输入文件 control

###input the phenotype name and output folder path 
Phenotype_name='RNAi cell cycle length'
phenotype_output_folder_path=os.path.abspath(r'C:\Users\dell\Desktop\output_folder\cell cycle phenotype')
if not os.path.exists(phenotype_output_folder_path):
    os.makedirs(phenotype_output_folder_path)

p_value = []
q_value_bh = []

control_std = np.std(control, ddof=1, axis=1)
control_mean = np.mean(control,axis=1)

diff=experiment.sub(control_mean,axis=0)

z_score=diff.divide(control_std,axis=0)
z_score=z_score.dropna(axis=1,how="all")
z_score.to_csv(os.path.join(phenotype_output_folder_path,Phenotype_name+"_Z-score.txt"), sep="\t") #out 输出文件 qvalue_bh

p_value=pd.DataFrame(1-scipy.stats.norm.cdf(abs(z_score)),index=z_score.index,columns=z_score.columns)

for i in p_value.columns:
    p_value_tmp=p_value[i]
    mask = np.isfinite(p_value_tmp)
    
    p_value_corrected=np.full(p_value_tmp.shape,np.nan)
    p_value_corrected[mask]=multipletests(p_value_tmp[mask], method='fdr_bh')[1]#Benjamini/Hochberg  correction
    
    if type(q_value_bh) == list:
        q_value_bh = pd.DataFrame(list(p_value_corrected))
    else:
        q_value_bh = pd.concat([q_value_bh, pd.DataFrame(list(p_value_corrected))], axis=1)

q_value_bh.columns = z_score.columns
q_value_bh.index = z_score.index
q_value_bh.to_csv(os.path.join(phenotype_output_folder_path,Phenotype_name+"_Q-value.txt"), sep="\t") 

mask_non_significant = q_value_bh >= 0.01
mask_significant = q_value_bh < 0.01

q_value_bh[mask_significant] = 1
q_value_bh[mask_non_significant] = 0

sign_matrix = np.sign(z_score.values)
sign_matrix_pd = pd.DataFrame(sign_matrix)
sign_matrix_pd = sign_matrix_pd.replace(0.0, 1.0)
q_value_bh = pd.DataFrame(sign_matrix_pd.values*(q_value_bh.values))

q_value_bh.columns = z_score.columns
q_value_bh.index = z_score.index
q_value_bh.to_csv(os.path.join(phenotype_output_folder_path,Phenotype_name+"_Q-value_binary_0.01.txt"), sep="\t")                 #out 输出文件 qvalue_bh

