# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 13:58:22 2022

@author: dell
"""


import scipy
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import os


###input the folder path of CTR divsion timing matrix used to define division asynchrony
CTR_distance_matrix_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_distance_matrix')

###input the folder path of RNAi divsion timing matrix used to define division asynchrony
RNAi_distance_matrix_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_distance_matrix')


###input the folder path of cell pair distance phenotype
output = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_different_matrix_AB256')
if not os.path.exists(output):
    os.makedirs(output)

########################################################################################################
cell_list=['ABalaaaala','ABalaaaalp','ABalaaaarl','ABalaaaarr','ABalaaapal','ABalaaapar','ABalaaappl','ABalaaappr','ABalaapaaa','ABalaapaap','ABalaapapa','ABalaapapp','ABalaappaa','ABalaappap','ABalaapppa','ABalaapppp','ABalapaaaa','ABalapaaap','ABalapaapa','ABalapaapp','ABalapapaa','ABalapapap','ABalapappa','ABalapappp','ABalappaaa','ABalappaap','ABalappapa','ABalappapp','ABalapppaa','ABalapppap','ABalappppa','ABalappppp','ABalpaaaaa','ABalpaaaap','ABalpaaapa','ABalpaaapp','ABalpaapaa','ABalpaapap','ABalpaappa','ABalpaappp','ABalpapaaa','ABalpapaap','ABalpapapa','ABalpapapp','ABalpappaa','ABalpappap','ABalpapppa','ABalpapppp','ABalppaaaa','ABalppaaap','ABalppaapa','ABalppaapp','ABalppapaa','ABalppapap','ABalppappa','ABalppappp','ABalpppaaa','ABalpppaap','ABalpppapa','ABalpppapp','ABalppppaa','ABalppppap','ABalpppppa','ABalpppppp','ABaraaaaaa','ABaraaaaap','ABaraaaapa','ABaraaaapp','ABaraaapaa','ABaraaapap','ABaraaappa','ABaraaappp','ABaraapaaa','ABaraapaap','ABaraapapa','ABaraapapp','ABaraappaa','ABaraappap','ABaraapppa','ABaraapppp','ABarapaaaa','ABarapaaap','ABarapaapa','ABarapaapp','ABarapapaa','ABarapapap','ABarapappa','ABarapappp','ABarappaaa','ABarappaap','ABarappapa','ABarappapp','ABarapppaa','ABarapppap','ABarappppa','ABarappppp','ABarpaaaaa','ABarpaaaap','ABarpaaapa','ABarpaaapp','ABarpaapaa','ABarpaapap','ABarpaappa','ABarpaappp','ABarpapaaa','ABarpapaap','ABarpapapa','ABarpapapp','ABarpappaa','ABarpappap','ABarpapppa','ABarpapppp','ABarppaaaa','ABarppaaap','ABarppaapa','ABarppaapp','ABarppapaa','ABarppapap','ABarppappa','ABarppappp','ABarpppaaa','ABarpppaap','ABarpppapa','ABarpppapp','ABarppppaa','ABarppppap','ABarpppppa','ABarpppppp','ABplaaaaaa','ABplaaaaap','ABplaaaapa','ABplaaaapp','ABplaaapaa','ABplaaapap','ABplaaappa','ABplaaappp','ABplaapaaa','ABplaapaap','ABplaapapa','ABplaapapp','ABplaappaa','ABplaappap','ABplaapppa','ABplaapppp','ABplapaaaa','ABplapaaap','ABplapaapa','ABplapaapp','ABplapapaa','ABplapapap','ABplapappa','ABplapappp','ABplappaaa','ABplappaap','ABplappapa','ABplappapp','ABplapppaa','ABplapppap','ABplappppa','ABplappppp','ABplpaaaaa','ABplpaaaap','ABplpaaapa','ABplpaaapp','ABplpaapaa','ABplpaapap','ABplpaappa','ABplpaappp','ABplpapaaa','ABplpapaap','ABplpapapa','ABplpapapp','ABplpappaa','ABplpappap','ABplpapppa','ABplpapppp','ABplppaaaa','ABplppaaap','ABplppaapa','ABplppaapp','ABplppapaa','ABplppapap','ABplppappa','ABplppappp','ABplpppaaa','ABplpppaap','ABplpppapa','ABplpppapp','ABplppppaa','ABplppppap','ABplpppppa','ABplpppppp','ABpraaaaaa','ABpraaaaap','ABpraaaapa','ABpraaaapp','ABpraaapaa','ABpraaapap','ABpraaappa','ABpraaappp','ABpraapaaa','ABpraapaap','ABpraapapa','ABpraapapp','ABpraappaa','ABpraappap','ABpraapppa','ABpraapppp','ABprapaaaa','ABprapaaap','ABprapaapa','ABprapaapp','ABprapapaa','ABprapapap','ABprapappa','ABprapappp','ABprappaaa','ABprappaap','ABprappapa','ABprappapp','ABprapppaa','ABprapppap','ABprappppa','ABprappppp','ABprpaaaaa','ABprpaaaap','ABprpaaapa','ABprpaaapp','ABprpaapaa','ABprpaapap','ABprpaappa','ABprpaappp','ABprpapaaa','ABprpapaap','ABprpapapa','ABprpapapp','ABprpappaa','ABprpappap','ABprpapppa','ABprpapppp','ABprppaaaa','ABprppaaap','ABprppaapa','ABprppaapp','ABprppapaa','ABprppapap','ABprppappa','ABprppappp','ABprpppaaa','ABprpppaap','ABprpppapa','ABprpppapp','ABprppppaa','ABprppppap','ABprpppppa','ABprpppppp','Caaaaa','Caaaap','Caaapa','Caaapp','Caapa','Caappd','Caappv','Capaaa','Capaap','Capapa','Capapp','Cappaa','Cappap','Capppa','Capppp','Cpaaaa','Cpaaap','Cpaapa','Cpaapp','Cpapaa','Cpapap','Cpappd','Cpappv','Cppaaa','Cppaap','Cppapa','Cppapp','Cpppaa','Cpppap','Cppppa','Cppppp','Daaa','Daap','Dapa','Dapp','Dpaa','Dpap','Dppa','Dppp','Ealaa','Ealap','Ealpa','Ealpp','Earaa','Earap','Earpa','Earpp','Eplaa','Eplap','Eplpa','Eplpp','Epraa','Eprap','Eprpa','Eprpp','MSaaaaaa','MSaaaaap','MSaaaapa','MSaaaapp','MSaaapaa','MSaaapap','MSaaapp','MSaapaaa','MSaapaap','MSaapapa','MSaapapp','MSaappa','MSaappp','MSapaaaa','MSapaaap','MSapaap','MSapapaa','MSapapap','MSapappa','MSapappp','MSappaa','MSappap','MSapppa','MSapppp','MSpaaaaa','MSpaaaap','MSpaaapa','MSpaaapp','MSpaapaa','MSpaapap','MSpaapp','MSpapaaa','MSpapaap','MSpapapa','MSpapapp','MSpappa','MSpappp','MSppaaaa','MSppaaap','MSppaap','MSppapaa','MSppapap','MSppappa','MSppappp','MSpppaa','MSpppap','MSppppa','MSppppp','Z2','Z3']

os.chdir(CTR_distance_matrix_path) 

cell_number=len(cell_list)
File_list=[]
embryo_list=[]
Embryo_cell_dist=[]
Embryo_cell_std=[]
Embryo_dist_sum = np.zeros(shape=(cell_number,cell_number))
Embryo_dist_unstack=[]
for file in os.listdir(CTR_distance_matrix_path):
    Embryo=str.split(file,"_distMatrix")[0]
    File_list.append(Embryo)
    Embryo_dist=pd.read_csv(file, index_col=0, sep='\t')
    Embryo_dist=Embryo_dist[cell_list].loc[cell_list]
    Embryo_dist_sum += Embryo_dist.values
    Embryo_dist=list(Embryo_dist.unstack())
    
    if type(Embryo_dist_unstack) == list:
        Embryo_dist_unstack = pd.DataFrame(list(Embryo_dist))                                  
    else:
        Embryo_dist_unstack = pd.concat([Embryo_dist_unstack, pd.DataFrame(list(Embryo_dist))], axis=1)
    Embryo_dist_std = np.std(Embryo_dist_unstack, ddof=1, axis=1)
    Embryo_dist_std=np.array(Embryo_dist_std).reshape(cell_number,cell_number)#更改细胞列表
        
Embryo_dist_mean=Embryo_dist_sum/105#平均的胚胎距离矩阵
Embryo_dist_mean=pd.DataFrame(Embryo_dist_mean)#更改细胞列表
Embryo_dist_mean.index=cell_list
Embryo_dist_mean.columns=cell_list
Embryo_dist_std=pd.DataFrame(Embryo_dist_std)#更改细胞列表
Embryo_dist_std.index=cell_list
Embryo_dist_std.columns=cell_list
        
##############################################################################        
def func(unit_list: list) -> np.array:
    """
    Generate a np.array by unit_list:
        [1, 2, 5] ->  [[nan 1   2  ]
                       [1   nan 5  ]
                       [2   5   nan]]
    """
    length = len(unit_list)
    for i in range(2*length):
        if 2*length == i*(i+1):
            n = i + 1

    list_to_array = []
    for i in range(1, n+1):
        zero = [0, ]*i
        zero.extend(unit_list[:n-i])
        list_to_array.append(zero)
        unit_list = unit_list[n-i:]

    array = np.array(list_to_array)
    array = array + array.T + np.diag((np.nan, )*n)
    return array
###############################################################################

File_list=[]
rmsd_list=[]
cell_position_difference=[]
Total_cell=pd.DataFrame(Embryo_dist_mean.index)
Total_cell=Total_cell.set_index(0)
os.chdir(RNAi_distance_matrix_path) 

a=0
for file in os.listdir(RNAi_distance_matrix_path):
    Embryo=str.split(file,"_distMatrix")[0]
    File_list.append(Embryo)
    Embryo_dist=pd.read_csv(file,index_col=0, sep='\t')
    Embryo_dist=Embryo_dist.replace(0,np.nan)#去除空值
    Embryo_dist=Embryo_dist.dropna(axis=0, how='all').dropna(axis=1, how='all')
    Embryo_dist=Embryo_dist.fillna(0)
    Embryo_dist=Embryo_dist
    coexist_cell_list=list(set(list(Embryo_dist_mean.index)).intersection(set(list(Embryo_dist.index))))
    coexist_cell_list.sort()
    if len(coexist_cell_list)>0:        
        Embryo_dist_mean_temp=Embryo_dist_mean[coexist_cell_list].loc[coexist_cell_list]
        Embryo_dist_std_temp=Embryo_dist_std[coexist_cell_list].loc[coexist_cell_list]
        Embryo_dist_temp=Embryo_dist[coexist_cell_list].loc[coexist_cell_list]
        Embryo_dist_zscore=(Embryo_dist_temp-Embryo_dist_mean_temp)/Embryo_dist_std_temp
        Embryo_dist_pvalue=pd.DataFrame(1-scipy.stats.norm.cdf(abs(Embryo_dist_zscore)),index=Embryo_dist_temp.index,columns=Embryo_dist_temp.columns)

        Embryo_dist_pvalue_triu=[Embryo_dist_pvalue.values[i,j] for i in range(len(Embryo_dist_pvalue)) for j in range(len(Embryo_dist_pvalue)) if i<j]
        Embryo_dist_qvalue_triu=multipletests(Embryo_dist_pvalue_triu, method='fdr_bh')[1]
        
        Embryo_dist_qvalue_mtrix=func(Embryo_dist_qvalue_triu)
        Embryo_dist_qvalue_mtrix=pd.DataFrame(Embryo_dist_qvalue_mtrix,
                                              index=Embryo_dist_pvalue.index,columns=Embryo_dist_pvalue.columns)#可导出每个胚胎差异细胞矩阵
                                              
        significant=Embryo_dist_qvalue_mtrix<0.01
        no_significant=Embryo_dist_qvalue_mtrix>=0.01
        Embryo_dist_qvalue_mtrix[significant]=1                                     
        Embryo_dist_qvalue_mtrix[no_significant]=0                                     
                                      
        Embryo_dist_qvalue_mtrix.to_csv(os.path.join(output,Embryo+'_AB256_comparison_Matrix.txt'), sep="\t")
        a+=1
        print(a)
        
        