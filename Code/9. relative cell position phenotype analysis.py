# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 20:31:12 2022

@author: dell
"""



import scipy
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import copy

###input the folder path of CTR divsion timing matrix used to define division asynchrony
CTR_distance_matrix_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_distance_matrix')

###input the folder path of RNAi divsion timing matrix used to define division asynchrony
RNAi_distance_matrix_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_distance_matrix')

###input the folder path of relative cell position phenotype
Relative_cell_position_phenotype_folder = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\relative_cell_position_phenotype')
if not os.path.exists(Relative_cell_position_phenotype_folder):
    os.makedirs(Relative_cell_position_phenotype_folder)
    

Stage_AB16=['ABalaa','ABalap','ABalpa','ABalpp','ABaraa','ABarap','ABarpa','ABarpp','ABplaa','ABplap','ABplpa','ABplpp','ABpraa','ABprap','ABprpa','ABprpp','Ca','Cp','D','Ea','Ep','MSaa','MSap','MSpa','MSpp','P4']
Stage_AB32=['ABalaaa','ABalaap','ABalapa','ABalapp','ABalpaa','ABalpap','ABalppa','ABalppp','ABaraaa','ABaraap','ABarapa','ABarapp','ABarpaa','ABarpap','ABarppa','ABarppp','ABplaaa','ABplaap','ABplapa','ABplapp','ABplpaa','ABplpap','ABplppa','ABplppp','ABpraaa','ABpraap','ABprapa','ABprapp','ABprpaa','ABprpap','ABprppa','ABprppp','Caa','Cap','Cpa','Cpp','Da','Dp','Eal','Ear','Epl','Epr','MSaaa','MSaap','MSapa','MSapp','MSpaa','MSpap','MSppa','MSppp','P4']
Stage_AB64=['ABalaaaa','ABalaaap','ABalaapa','ABalaapp','ABalapaa','ABalapap','ABalappa','ABalappp','ABalpaaa','ABalpaap','ABalpapa','ABalpapp','ABalppaa','ABalppap','ABalpppa','ABalpppp','ABaraaaa','ABaraaap','ABaraapa','ABaraapp','ABarapaa','ABarapap','ABarappa','ABarappp','ABarpaaa','ABarpaap','ABarpapa','ABarpapp','ABarppaa','ABarppap','ABarpppa','ABarpppp','ABplaaaa','ABplaaap','ABplaapa','ABplaapp','ABplapaa','ABplapap','ABplappa','ABplappp','ABplpaaa','ABplpaap','ABplpapa','ABplpapp','ABplppaa','ABplppap','ABplpppa','ABplpppp','ABpraaaa','ABpraaap','ABpraapa','ABpraapp','ABprapaa','ABprapap','ABprappa','ABprappp','ABprpaaa','ABprpaap','ABprpapa','ABprpapp','ABprppaa','ABprppap','ABprpppa','ABprpppp','Caaa','Caap','Capa','Capp','Cpaa','Cpap','Cppa','Cppp','Da','Dp','Eal','Ear','Epl','Epr','MSaaaa','MSaaap','MSaapa','MSaapp','MSapaa','MSapap','MSappa','MSappp','MSpaaa','MSpaap','MSpapa','MSpapp','MSppaa','MSppap','MSpppa','MSpppp','Z2','Z3']
Stage_AB128=['ABalaaaal','ABalaaaar','ABalaaapa','ABalaaapp','ABalaapaa','ABalaapap','ABalaappa','ABalaappp','ABalapaaa','ABalapaap','ABalapapa','ABalapapp','ABalappaa','ABalappap','ABalapppa','ABalapppp','ABalpaaaa','ABalpaaap','ABalpaapa','ABalpaapp','ABalpapaa','ABalpapap','ABalpappa','ABalpappp','ABalppaaa','ABalppaap','ABalppapa','ABalppapp','ABalpppaa','ABalpppap','ABalppppa','ABalppppp','ABaraaaaa','ABaraaaap','ABaraaapa','ABaraaapp','ABaraapaa','ABaraapap','ABaraappa','ABaraappp','ABarapaaa','ABarapaap','ABarapapa','ABarapapp','ABarappaa','ABarappap','ABarapppa','ABarapppp','ABarpaaaa','ABarpaaap','ABarpaapa','ABarpaapp','ABarpapaa','ABarpapap','ABarpappa','ABarpappp','ABarppaaa','ABarppaap','ABarppapa','ABarppapp','ABarpppaa','ABarpppap','ABarppppa','ABarppppp','ABplaaaaa','ABplaaaap','ABplaaapa','ABplaaapp','ABplaapaa','ABplaapap','ABplaappa','ABplaappp','ABplapaaa','ABplapaap','ABplapapa','ABplapapp','ABplappaa','ABplappap','ABplapppa','ABplapppp','ABplpaaaa','ABplpaaap','ABplpaapa','ABplpaapp','ABplpapaa','ABplpapap','ABplpappa','ABplpappp','ABplppaaa','ABplppaap','ABplppapa','ABplppapp','ABplpppaa','ABplpppap','ABplppppa','ABplppppp','ABpraaaaa','ABpraaaap','ABpraaapa','ABpraaapp','ABpraapaa','ABpraapap','ABpraappa','ABpraappp','ABprapaaa','ABprapaap','ABprapapa','ABprapapp','ABprappaa','ABprappap','ABprapppa','ABprapppp','ABprpaaaa','ABprpaaap','ABprpaapa','ABprpaapp','ABprpapaa','ABprpapap','ABprpappa','ABprpappp','ABprppaaa','ABprppaap','ABprppapa','ABprppapp','ABprpppaa','ABprpppap','ABprppppa','ABprppppp','MSaaaaa','MSaaaap','MSaaapa','MSaaapp','MSaapaa','MSaapap','MSaappa','MSaappp','MSapaaa','MSapaap','MSapapa','MSapapp','MSappaa','MSappap','MSapppa','MSapppp','MSpaaaa','MSpaaap','MSpaapa','MSpaapp','MSpapaa','MSpapap','MSpappa','MSpappp','MSppaaa','MSppaap','MSppapa','MSppapp','MSpppaa','MSpppap','MSppppa','MSppppp','Eala','Ealp','Eara','Earp','Epla','Eplp','Epra','Eprp','Caaaa','Caaap','Caapa','Caapp','Capaa','Capap','Cappa','Cappp','Cpaaa','Cpaap','Cpapa','Cpapp','Cppaa','Cppap','Cpppa','Cpppp','Daa','Dap','Dpa','Dpp','Z2','Z3']
Stage_AB256=['ABalaaaala','ABalaaaalp','ABalaaaarl','ABalaaaarr','ABalaaapal','ABalaaapar','ABalaaappl','ABalaaappr','ABalaapaaa','ABalaapaap','ABalaapapa','ABalaapapp','ABalaappaa','ABalaappap','ABalaapppa','ABalaapppp','ABalapaaaa','ABalapaaap','ABalapaapa','ABalapaapp','ABalapapaa','ABalapapap','ABalapappa','ABalapappp','ABalappaaa','ABalappaap','ABalappapa','ABalappapp','ABalapppaa','ABalapppap','ABalappppa','ABalappppp','ABalpaaaaa','ABalpaaaap','ABalpaaapa','ABalpaaapp','ABalpaapaa','ABalpaapap','ABalpaappa','ABalpaappp','ABalpapaaa','ABalpapaap','ABalpapapa','ABalpapapp','ABalpappaa','ABalpappap','ABalpapppa','ABalpapppp','ABalppaaaa','ABalppaaap','ABalppaapa','ABalppaapp','ABalppapaa','ABalppapap','ABalppappa','ABalppappp','ABalpppaaa','ABalpppaap','ABalpppapa','ABalpppapp','ABalppppaa','ABalppppap','ABalpppppa','ABalpppppp','ABaraaaaaa','ABaraaaaap','ABaraaaapa','ABaraaaapp','ABaraaapaa','ABaraaapap','ABaraaappa','ABaraaappp','ABaraapaaa','ABaraapaap','ABaraapapa','ABaraapapp','ABaraappaa','ABaraappap','ABaraapppa','ABaraapppp','ABarapaaaa','ABarapaaap','ABarapaapa','ABarapaapp','ABarapapaa','ABarapapap','ABarapappa','ABarapappp','ABarappaaa','ABarappaap','ABarappapa','ABarappapp','ABarapppaa','ABarapppap','ABarappppa','ABarappppp','ABarpaaaaa','ABarpaaaap','ABarpaaapa','ABarpaaapp','ABarpaapaa','ABarpaapap','ABarpaappa','ABarpaappp','ABarpapaaa','ABarpapaap','ABarpapapa','ABarpapapp','ABarpappaa','ABarpappap','ABarpapppa','ABarpapppp','ABarppaaaa','ABarppaaap','ABarppaapa','ABarppaapp','ABarppapaa','ABarppapap','ABarppappa','ABarppappp','ABarpppaaa','ABarpppaap','ABarpppapa','ABarpppapp','ABarppppaa','ABarppppap','ABarpppppa','ABarpppppp','ABplaaaaaa','ABplaaaaap','ABplaaaapa','ABplaaaapp','ABplaaapaa','ABplaaapap','ABplaaappa','ABplaaappp','ABplaapaaa','ABplaapaap','ABplaapapa','ABplaapapp','ABplaappaa','ABplaappap','ABplaapppa','ABplaapppp','ABplapaaaa','ABplapaaap','ABplapaapa','ABplapaapp','ABplapapaa','ABplapapap','ABplapappa','ABplapappp','ABplappaaa','ABplappaap','ABplappapa','ABplappapp','ABplapppaa','ABplapppap','ABplappppa','ABplappppp','ABplpaaaaa','ABplpaaaap','ABplpaaapa','ABplpaaapp','ABplpaapaa','ABplpaapap','ABplpaappa','ABplpaappp','ABplpapaaa','ABplpapaap','ABplpapapa','ABplpapapp','ABplpappaa','ABplpappap','ABplpapppa','ABplpapppp','ABplppaaaa','ABplppaaap','ABplppaapa','ABplppaapp','ABplppapaa','ABplppapap','ABplppappa','ABplppappp','ABplpppaaa','ABplpppaap','ABplpppapa','ABplpppapp','ABplppppaa','ABplppppap','ABplpppppa','ABplpppppp','ABpraaaaaa','ABpraaaaap','ABpraaaapa','ABpraaaapp','ABpraaapaa','ABpraaapap','ABpraaappa','ABpraaappp','ABpraapaaa','ABpraapaap','ABpraapapa','ABpraapapp','ABpraappaa','ABpraappap','ABpraapppa','ABpraapppp','ABprapaaaa','ABprapaaap','ABprapaapa','ABprapaapp','ABprapapaa','ABprapapap','ABprapappa','ABprapappp','ABprappaaa','ABprappaap','ABprappapa','ABprappapp','ABprapppaa','ABprapppap','ABprappppa','ABprappppp','ABprpaaaaa','ABprpaaaap','ABprpaaapa','ABprpaaapp','ABprpaapaa','ABprpaapap','ABprpaappa','ABprpaappp','ABprpapaaa','ABprpapaap','ABprpapapa','ABprpapapp','ABprpappaa','ABprpappap','ABprpapppa','ABprpapppp','ABprppaaaa','ABprppaaap','ABprppaapa','ABprppaapp','ABprppapaa','ABprppapap','ABprppappa','ABprppappp','ABprpppaaa','ABprpppaap','ABprpppapa','ABprpppapp','ABprppppaa','ABprppppap','ABprpppppa','ABprpppppp','Caaaaa','Caaaap','Caaapa','Caaapp','Caapa','Caappd','Caappv','Capaaa','Capaap','Capapa','Capapp','Cappaa','Cappap','Capppa','Capppp','Cpaaaa','Cpaaap','Cpaapa','Cpaapp','Cpapaa','Cpapap','Cpappd','Cpappv','Cppaaa','Cppaap','Cppapa','Cppapp','Cpppaa','Cpppap','Cppppa','Cppppp','Daaa','Daap','Dapa','Dapp','Dpaa','Dpap','Dppa','Dppp','Ealaa','Ealap','Ealpa','Ealpp','Earaa','Earap','Earpa','Earpp','Eplaa','Eplap','Eplpa','Eplpp','Epraa','Eprap','Eprpa','Eprpp','MSaaaaaa','MSaaaaap','MSaaaapa','MSaaaapp','MSaaapaa','MSaaapap','MSaaapp','MSaapaaa','MSaapaap','MSaapapa','MSaapapp','MSaappa','MSaappp','MSapaaaa','MSapaaap','MSapaap','MSapapaa','MSapapap','MSapappa','MSapappp','MSappaa','MSappap','MSapppa','MSapppp','MSpaaaaa','MSpaaaap','MSpaaapa','MSpaaapp','MSpaapaa','MSpaapap','MSpaapp','MSpapaaa','MSpapaap','MSpapapa','MSpapapp','MSpappa','MSpappp','MSppaaaa','MSppaaap','MSppaap','MSppapaa','MSppapap','MSppappa','MSppappp','MSpppaa','MSpppap','MSppppa','MSppppp','Z2','Z3']
Embryo_dist_mean_WT=pd.DataFrame()
RNAi_cellposition_zscore=pd.DataFrame()
WT_single_cell_RMSD_all=[]
RNAi_cellposition_qvalue_binary=[]
RNAi_cellposition_qvalue=[]  
RNAi_single_cell_RMSD_all=[] 
for stage in [Stage_AB16,Stage_AB32,Stage_AB64,Stage_AB128,Stage_AB256]:

    cell_list=stage
    
    os.chdir (CTR_distance_matrix_path)#control cell pair distance matrix
    
    File_list=[]
    Total_matrix=[]
    for file in os.listdir(CTR_distance_matrix_path):
        Embryo=str.split(file,"_distMatrix_")[0]
        File_list.append(Embryo)
        Embryo_dist=pd.read_csv(file,index_col=0,header=0,sep='\t')
        Embryo_dist=Embryo_dist[cell_list].loc[cell_list]
        if type(Total_matrix) == list:
            Total_matrix = Embryo_dist                                 
        else:
            Total_matrix = pd.concat([Total_matrix, Embryo_dist], axis=1)

    cell_rmsd_mean=[]#Calculate the mean and standard deviation of the RMSD per cell between wild-type embryos 
    cell_rmsd_sd=[]
    Embryo_dist_mean=[]        
    for cell in cell_list:
        cell_vector_matrix=Total_matrix[cell] 
        cell_vector_mean=np.mean(cell_vector_matrix,axis=1)       
        cell_vector_matrix=cell_vector_matrix.drop([cell],axis=0)
        cell_vector_matrix=cell_vector_matrix.T
        cell_rmsd=scipy.spatial.distance.pdist(cell_vector_matrix)/len(cell_vector_matrix.columns)**0.5#caculate rmsd       
        
        cell_rmsd = [x for x in cell_rmsd if not np.isnan(x)]    
        
        cell_rmsd_mean_temp=np.average(cell_rmsd)
        cell_rmsd_sd_temp=np.std(cell_rmsd)
        cell_rmsd_mean.append(cell_rmsd_mean_temp)
        cell_rmsd_sd.append(cell_rmsd_sd_temp)
        if type(Embryo_dist_mean)==list:
            Embryo_dist_mean=pd.DataFrame(cell_vector_mean)
        else:
            Embryo_dist_mean=pd.concat([Embryo_dist_mean,pd.DataFrame(cell_vector_mean)],axis=1)
    Embryo_dist_mean.index=cell_list
    Embryo_dist_mean.columns=cell_list 
    Embryo_dist_mean_WT=pd.concat([Embryo_dist_mean_WT,Embryo_dist_mean],axis=0)
      
    WT_single_cell_RMSD=pd.DataFrame()
    WT_single_cell_RMSD['cell_rmsd_mean']=cell_rmsd_mean
    WT_single_cell_RMSD['cell_rmsd_sd']=cell_rmsd_sd
    WT_single_cell_RMSD.index=cell_list
    
    if type(WT_single_cell_RMSD_all)==list:
        WT_single_cell_RMSD_all=WT_single_cell_RMSD
    else:
        WT_single_cell_RMSD_all=pd.concat([WT_single_cell_RMSD_all,WT_single_cell_RMSD],axis=0)
    #########################################################################################################      
    os.chdir (RNAi_distance_matrix_path)  #RNAi cell pair distance matrix
    embryo_list=[]
    single_cell_rmsd_matrix=[]
    a=0
    for file in os.listdir(RNAi_distance_matrix_path):
        Embryo=str.split(file,"_distMatrix_")[0]
        embryo_list.append(Embryo)
        Embryo_dist=pd.read_csv(file, index_col=0,header=0,sep='\t')
        Embryo_dist=Embryo_dist[cell_list].loc[cell_list]
        single_cell_RMSD=[]
        for cell in cell_list:
            cell_vactor_temp=Embryo_dist[cell]
            cell_vactor_mean=Embryo_dist_mean[cell]
            cell_vactor=pd.concat([pd.DataFrame(cell_vactor_temp),pd.DataFrame(cell_vactor_mean)],axis=1)
            cell_vactor=cell_vactor.drop([cell]).dropna()
            if len(cell_vactor)>0:
                cell_vactor=cell_vactor.T
                cell_rmsd=scipy.spatial.distance.pdist(cell_vactor)/len(cell_vactor.columns)**0.5
                single_cell_RMSD.append(cell_rmsd)
            else:
                cell_rmsd=np.array([np.nan])
                single_cell_RMSD.append(cell_rmsd)
                
            
        if type(single_cell_rmsd_matrix) == list:
            single_cell_rmsd_matrix = pd.DataFrame(single_cell_RMSD)                                  
        else:
            single_cell_rmsd_matrix = pd.concat([single_cell_rmsd_matrix, pd.DataFrame(single_cell_RMSD)], axis=1)
        a+=1
        print(a)
    
    single_cell_rmsd_matrix.index=cell_list
    single_cell_rmsd_matrix.columns=embryo_list
    if type(RNAi_single_cell_RMSD_all)==list:
        RNAi_single_cell_RMSD_all=single_cell_rmsd_matrix
    else:
        RNAi_single_cell_RMSD_all=pd.concat([RNAi_single_cell_RMSD_all,single_cell_rmsd_matrix],axis=0) 
        
    p_value = []
    q_value_bh = []
    
    control_std = WT_single_cell_RMSD['cell_rmsd_sd']
    control_mean = WT_single_cell_RMSD['cell_rmsd_mean']
    
    diff=single_cell_rmsd_matrix.sub(control_mean,axis=0)
    
    z_score=diff.divide(control_std,axis=0)
    RNAi_cellposition_zscore=pd.concat([RNAi_cellposition_zscore,pd.DataFrame(z_score)],axis=0)
    p_value=pd.DataFrame(1-scipy.stats.norm.cdf(abs(z_score)),index=z_score.index,columns=z_score.columns)
    
    for i in p_value.columns:
        p_value_tmp=p_value[i]
        mask = np.isfinite(p_value_tmp)    
        p_value_corrected=np.full(p_value_tmp.shape,np.nan)
        value_numbers=p_value_tmp[mask]
        if len(value_numbers)>0:
            p_value_corrected[mask]=multipletests(p_value_tmp[mask], method='fdr_bh')[1] 
        else:
            p_value_corrected[mask]=p_value_tmp[mask]
        if type(q_value_bh) == list:
            q_value_bh = pd.DataFrame(list(p_value_corrected))                                   
        else:
            q_value_bh = pd.concat([q_value_bh, pd.DataFrame(list(p_value_corrected))], axis=1)
    
    q_value_bh.columns = z_score.columns
    q_value_bh.index = z_score.index

    if type(RNAi_cellposition_qvalue)==list:
        RNAi_cellposition_qvalue=q_value_bh
    else:
        RNAi_cellposition_qvalue=pd.concat([RNAi_cellposition_qvalue,q_value_bh],axis=0)  

    q_value_bh_copy=copy.deepcopy(q_value_bh)     
    mask_non_significant = q_value_bh_copy >= 0.01
    mask_significant = q_value_bh_copy < 0.01
    
    q_value_bh_copy[mask_significant] = 1
    q_value_bh_copy[mask_non_significant] = 0
    
    sign_matrix = np.sign(z_score.values)
    sign_matrix_pd = pd.DataFrame(sign_matrix)
    sign_matrix_pd = sign_matrix_pd.replace(0.0, 1.0)
    q_value_bh_copy = pd.DataFrame(sign_matrix_pd.values*(q_value_bh_copy.values))
        
    q_value_bh_copy.columns = z_score.columns
    q_value_bh_copy.index = z_score.index    
    q_value_bh_copy=q_value_bh_copy.replace(-1,0)

    if type(RNAi_cellposition_qvalue_binary)==list:
        RNAi_cellposition_qvalue_binary=q_value_bh_copy
    else:
        RNAi_cellposition_qvalue_binary=pd.concat([RNAi_cellposition_qvalue_binary,q_value_bh_copy],axis=0)    

Embryo_dist_mean_WT.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'CTR_relative_cell_pair_distance_vector_mean.txt'), sep="\t")   
WT_single_cell_RMSD_all.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'CTR_relative_cell_position_RMSD_mean.txt'), sep="\t")   

RNAi_single_cell_RMSD_all=RNAi_single_cell_RMSD_all[~RNAi_single_cell_RMSD_all.index.duplicated(keep='first')]     
RNAi_single_cell_RMSD_all.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'RNAi_relative_cell_position_RMSD.txt'), sep="\t")

RNAi_cellposition_zscore=RNAi_cellposition_zscore[~RNAi_cellposition_zscore.index.duplicated(keep='first')]     
RNAi_cellposition_zscore.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'RNAi_relative_cell_position_Z-score.txt'), sep="\t")

RNAi_cellposition_qvalue=RNAi_cellposition_qvalue[~RNAi_cellposition_qvalue.index.duplicated(keep='first')]      
RNAi_cellposition_qvalue.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'RNAi_relative_cell_position_Q-value.txt'), sep="\t")

RNAi_cellposition_qvalue_binary=RNAi_cellposition_qvalue_binary[~RNAi_cellposition_qvalue_binary.index.duplicated(keep='first')]   
RNAi_cellposition_qvalue_binary.to_csv(os.path.join(Relative_cell_position_phenotype_folder,'RNAi_relative_cell_position_Q-value_binary_0.01.txt'), sep="\t")




