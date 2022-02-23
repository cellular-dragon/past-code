# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 18:52:32 2022

@author: dell
"""

import pandas as pd
import numpy as np

#input the path and file of control CND-1 expression binary matrix
CTR_exp_matrix = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\CND-1_expression_value_matrix_binary.txt", sep="\t", index_col=0, header=0)

#input the path and file of RNAi CND-1 expression binary matrix
RNAi_exp_matrix = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\CND-1_expression_value_matrix_binary.txt", sep="\t", index_col=0, header=0)

CTR_exp=np.sum(CTR_exp_matrix,axis=1)/len(CTR_exp_matrix.columns)
CTR_exp_cell_list=list(CTR_exp[CTR_exp>=0.99].index)
CTR_noexp_cell_list=list(CTR_exp[CTR_exp<=0.01].index)
CTR_cell_list=CTR_exp_cell_list + CTR_noexp_cell_list
RNAi_exp=RNAi_exp_matrix.loc[CTR_exp_cell_list]
RNAi_noexp=RNAi_exp_matrix.loc[CTR_noexp_cell_list]
RNAi_exp_differ=RNAi_exp-1
RNAi_noexp_differ=RNAi_noexp-0
RNAi_differ= pd.concat([RNAi_exp_differ, RNAi_noexp_differ], axis=0)
sign_matrix = np.sign(RNAi_differ.values)
sign_matrix_pd = pd.DataFrame(sign_matrix)
sign_matrix_pd = sign_matrix_pd.replace(0.0, 1.0)

#input the path and file name of CND-1 expression phenotype
RNAi_differ.to_csv(r'C:\\Users\\dell\\Desktop\\CND-1_expression_phenotype.txt', sep='\t')






