# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 13:54:11 2022

@author: dell
"""

import os
import pandas as pd
import numpy as np
import scipy

"""
The script runs in a large memory, so you are advised to run it on a computer with a memory of more than 16 GB.  
"""
File_Directory = r'C:\Users\dell\Desktop\raw_data'   # input raw data of control or RNAi or RNAi_add

###input the folder path of distance matrix used to define relative cell position
distance_matrix_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_distance_matrix')
if not os.path.exists(distance_matrix_path):
    os.makedirs(distance_matrix_path)

###input the folder path which store all cell XYZ coordinates of each embryo used to analyze the cell position correction in the simulated case 
XYZ_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder\RNAi_XYZ')
if not os.path.exists(distance_matrix_path):    
    os.makedirs(distance_matrix_path)    
###No need to change, store a list of embryo lengths in the parent directory of distance_matrix_path 
embryo_length_path = os.path.abspath(os.path.join(distance_matrix_path, ".."))

###input the embryo length file name
embryo_length_file_name='RNAi_embryo_length'


os.chdir(File_Directory)
raw_data_list = os.listdir(File_Directory)
Cell_list_714=['ABal', 'ABala', 'ABalaa', 'ABalaaa', 'ABalaaaa', 'ABalaaaal', 'ABalaaaala', 'ABalaaaalp', 'ABalaaaar', 'ABalaaaarl', 'ABalaaaarr', 'ABalaaap', 'ABalaaapa', 'ABalaaapal', 'ABalaaapar', 'ABalaaapp', 'ABalaaappl', 'ABalaaappr', 'ABalaap', 'ABalaapa', 'ABalaapaa', 'ABalaapaaa', 'ABalaapaap', 'ABalaapap', 'ABalaapapa', 'ABalaapapp', 'ABalaapp', 'ABalaappa', 'ABalaappaa', 'ABalaappap', 'ABalaappp', 'ABalaapppa', 'ABalaapppp', 'ABalap', 'ABalapa', 'ABalapaa', 'ABalapaaa', 'ABalapaaaa', 'ABalapaaap', 'ABalapaap', 'ABalapaapa', 'ABalapaapp', 'ABalapap', 'ABalapapa', 'ABalapapaa', 'ABalapapap', 'ABalapapp', 'ABalapappa', 'ABalapappp', 'ABalapp', 'ABalappa', 'ABalappaa', 'ABalappaaa', 'ABalappaap', 'ABalappap', 'ABalappapa', 'ABalappapp', 'ABalappp', 'ABalapppa', 'ABalapppaa', 'ABalapppap', 'ABalapppp', 'ABalappppa', 'ABalappppp', 'ABalp', 'ABalpa', 'ABalpaa', 'ABalpaaa', 'ABalpaaaa', 'ABalpaaaaa', 'ABalpaaaap', 'ABalpaaap', 'ABalpaaapa', 'ABalpaaapp', 'ABalpaap', 'ABalpaapa', 'ABalpaapaa', 'ABalpaapap', 'ABalpaapp', 'ABalpaappa', 'ABalpaappp', 'ABalpap', 'ABalpapa', 'ABalpapaa', 'ABalpapaaa', 'ABalpapaap', 'ABalpapap', 'ABalpapapa', 'ABalpapapp', 'ABalpapp', 'ABalpappa', 'ABalpappaa', 'ABalpappap', 'ABalpappp', 'ABalpapppa', 'ABalpapppp', 'ABalpp', 'ABalppa', 'ABalppaa', 'ABalppaaa', 'ABalppaaaa', 'ABalppaaap', 'ABalppaap', 'ABalppaapa', 'ABalppaapp', 'ABalppap', 'ABalppapa', 'ABalppapaa', 'ABalppapap', 'ABalppapp', 'ABalppappa', 'ABalppappp', 'ABalppp', 'ABalpppa', 'ABalpppaa', 'ABalpppaaa', 'ABalpppaap', 'ABalpppap', 'ABalpppapa', 'ABalpppapp', 'ABalpppp', 'ABalppppa', 'ABalppppaa', 'ABalppppap', 'ABalppppp', 'ABalpppppa', 'ABalpppppp', 'ABar', 'ABara', 'ABaraa', 'ABaraaa', 'ABaraaaa', 'ABaraaaaa', 'ABaraaaaaa', 'ABaraaaaap', 'ABaraaaap', 'ABaraaaapa', 'ABaraaaapp', 'ABaraaap', 'ABaraaapa', 'ABaraaapaa', 'ABaraaapap', 'ABaraaapp', 'ABaraaappa', 'ABaraaappp', 'ABaraap', 'ABaraapa', 'ABaraapaa', 'ABaraapaaa', 'ABaraapaap', 'ABaraapap', 'ABaraapapa', 'ABaraapapp', 'ABaraapp', 'ABaraappa', 'ABaraappaa', 'ABaraappap', 'ABaraappp', 'ABaraapppa', 'ABaraapppp', 'ABarap', 'ABarapa', 'ABarapaa', 'ABarapaaa', 'ABarapaaaa', 'ABarapaaap', 'ABarapaap', 'ABarapaapa', 'ABarapaapp', 'ABarapap', 'ABarapapa', 'ABarapapaa', 'ABarapapap', 'ABarapapp', 'ABarapappa', 'ABarapappp', 'ABarapp', 'ABarappa', 'ABarappaa', 'ABarappaaa', 'ABarappaap', 'ABarappap', 'ABarappapa', 'ABarappapp', 'ABarappp', 'ABarapppa', 'ABarapppaa', 'ABarapppap', 'ABarapppp', 'ABarappppa', 'ABarappppp', 'ABarp', 'ABarpa', 'ABarpaa', 'ABarpaaa', 'ABarpaaaa', 'ABarpaaaaa', 'ABarpaaaap', 'ABarpaaap', 'ABarpaaapa', 'ABarpaaapp', 'ABarpaap', 'ABarpaapa', 'ABarpaapaa', 'ABarpaapap', 'ABarpaapp', 'ABarpaappa', 'ABarpaappp', 'ABarpap', 'ABarpapa', 'ABarpapaa', 'ABarpapaaa', 'ABarpapaap', 'ABarpapap', 'ABarpapapa', 'ABarpapapp', 'ABarpapp', 'ABarpappa', 'ABarpappaa', 'ABarpappap', 'ABarpappp', 'ABarpapppa', 'ABarpapppp', 'ABarpp', 'ABarppa', 'ABarppaa', 'ABarppaaa', 'ABarppaaaa', 'ABarppaaap', 'ABarppaap', 'ABarppaapa', 'ABarppaapp', 'ABarppap', 'ABarppapa', 'ABarppapaa', 'ABarppapap', 'ABarppapp', 'ABarppappa', 'ABarppappp', 'ABarppp', 'ABarpppa', 'ABarpppaa', 'ABarpppaaa', 'ABarpppaap', 'ABarpppap', 'ABarpppapa', 'ABarpppapp', 'ABarpppp', 'ABarppppa', 'ABarppppaa', 'ABarppppap', 'ABarppppp', 'ABarpppppa', 'ABarpppppp', 'ABpl', 'ABpla', 'ABplaa', 'ABplaaa', 'ABplaaaa', 'ABplaaaaa', 'ABplaaaaaa', 'ABplaaaaap', 'ABplaaaap', 'ABplaaaapa', 'ABplaaaapp', 'ABplaaap', 'ABplaaapa', 'ABplaaapaa', 'ABplaaapap', 'ABplaaapp', 'ABplaaappa', 'ABplaaappp', 'ABplaap', 'ABplaapa', 'ABplaapaa', 'ABplaapaaa', 'ABplaapaap', 'ABplaapap', 'ABplaapapa', 'ABplaapapp', 'ABplaapp', 'ABplaappa', 'ABplaappaa', 'ABplaappap', 'ABplaappp', 'ABplaapppa', 'ABplaapppp', 'ABplap', 'ABplapa', 'ABplapaa', 'ABplapaaa', 'ABplapaaaa', 'ABplapaaap', 'ABplapaap', 'ABplapaapa', 'ABplapaapp', 'ABplapap', 'ABplapapa', 'ABplapapaa', 'ABplapapap', 'ABplapapp', 'ABplapappa', 'ABplapappp', 'ABplapp', 'ABplappa', 'ABplappaa', 'ABplappaaa', 'ABplappaap', 'ABplappap', 'ABplappapa', 'ABplappapp', 'ABplappp', 'ABplapppa', 'ABplapppaa', 'ABplapppap', 'ABplapppp', 'ABplappppa', 'ABplappppp', 'ABplp', 'ABplpa', 'ABplpaa', 'ABplpaaa', 'ABplpaaaa', 'ABplpaaaaa', 'ABplpaaaap', 'ABplpaaap', 'ABplpaaapa', 'ABplpaaapp', 'ABplpaap', 'ABplpaapa', 'ABplpaapaa', 'ABplpaapap', 'ABplpaapp', 'ABplpaappa', 'ABplpaappp', 'ABplpap', 'ABplpapa', 'ABplpapaa', 'ABplpapaaa', 'ABplpapaap', 'ABplpapap', 'ABplpapapa', 'ABplpapapp', 'ABplpapp', 'ABplpappa', 'ABplpappaa', 'ABplpappap', 'ABplpappp', 'ABplpapppa', 'ABplpapppp', 'ABplpp', 'ABplppa', 'ABplppaa', 'ABplppaaa', 'ABplppaaaa', 'ABplppaaap', 'ABplppaap', 'ABplppaapa', 'ABplppaapp', 'ABplppap', 'ABplppapa', 'ABplppapaa', 'ABplppapap', 'ABplppapp', 'ABplppappa', 'ABplppappp', 'ABplppp', 'ABplpppa', 'ABplpppaa', 'ABplpppaaa', 'ABplpppaap', 'ABplpppap', 'ABplpppapa', 'ABplpppapp', 'ABplpppp', 'ABplppppa', 'ABplppppaa', 'ABplppppap', 'ABplppppp', 'ABplpppppa', 'ABplpppppp', 'ABpr', 'ABpra', 'ABpraa', 'ABpraaa', 'ABpraaaa', 'ABpraaaaa', 'ABpraaaaaa', 'ABpraaaaap', 'ABpraaaap', 'ABpraaaapa', 'ABpraaaapp', 'ABpraaap', 'ABpraaapa', 'ABpraaapaa', 'ABpraaapap', 'ABpraaapp', 'ABpraaappa', 'ABpraaappp', 'ABpraap', 'ABpraapa', 'ABpraapaa', 'ABpraapaaa', 'ABpraapaap', 'ABpraapap', 'ABpraapapa', 'ABpraapapp', 'ABpraapp', 'ABpraappa', 'ABpraappaa', 'ABpraappap', 'ABpraappp', 'ABpraapppa', 'ABpraapppp', 'ABprap', 'ABprapa', 'ABprapaa', 'ABprapaaa', 'ABprapaaaa', 'ABprapaaap', 'ABprapaap', 'ABprapaapa', 'ABprapaapp', 'ABprapap', 'ABprapapa', 'ABprapapaa', 'ABprapapap', 'ABprapapp', 'ABprapappa', 'ABprapappp', 'ABprapp', 'ABprappa', 'ABprappaa', 'ABprappaaa', 'ABprappaap', 'ABprappap', 'ABprappapa', 'ABprappapp', 'ABprappp', 'ABprapppa', 'ABprapppaa', 'ABprapppap', 'ABprapppp', 'ABprappppa', 'ABprappppp', 'ABprp', 'ABprpa', 'ABprpaa', 'ABprpaaa', 'ABprpaaaa', 'ABprpaaaaa', 'ABprpaaaap', 'ABprpaaap', 'ABprpaaapa', 'ABprpaaapp', 'ABprpaap', 'ABprpaapa', 'ABprpaapaa', 'ABprpaapap', 'ABprpaapp', 'ABprpaappa', 'ABprpaappp', 'ABprpap', 'ABprpapa', 'ABprpapaa', 'ABprpapaaa', 'ABprpapaap', 'ABprpapap', 'ABprpapapa', 'ABprpapapp', 'ABprpapp', 'ABprpappa', 'ABprpappaa', 'ABprpappap', 'ABprpappp', 'ABprpapppa', 'ABprpapppp', 'ABprpp', 'ABprppa', 'ABprppaa', 'ABprppaaa', 'ABprppaaaa', 'ABprppaaap', 'ABprppaap', 'ABprppaapa', 'ABprppaapp', 'ABprppap', 'ABprppapa', 'ABprppapaa', 'ABprppapap', 'ABprppapp', 'ABprppappa', 'ABprppappp', 'ABprppp', 'ABprpppa', 'ABprpppaa', 'ABprpppaaa', 'ABprpppaap', 'ABprpppap', 'ABprpppapa', 'ABprpppapp', 'ABprpppp', 'ABprppppa', 'ABprppppaa', 'ABprppppap', 'ABprppppp', 'ABprpppppa', 'ABprpppppp', 'C', 'Ca', 'Caa', 'Caaa', 'Caaaa', 'Caaaaa', 'Caaaap', 'Caaap', 'Caaapa', 'Caaapp', 'Caap', 'Caapa', 'Caapp', 'Caappd', 'Caappv', 'Cap', 'Capa', 'Capaa', 'Capaaa', 'Capaap', 'Capap', 'Capapa', 'Capapp', 'Capp', 'Cappa', 'Cappaa', 'Cappap', 'Cappp', 'Capppa', 'Capppp', 'Cp', 'Cpa', 'Cpaa', 'Cpaaa', 'Cpaaaa', 'Cpaaap', 'Cpaap', 'Cpaapa', 'Cpaapp', 'Cpap', 'Cpapa', 'Cpapaa', 'Cpapap', 'Cpapp', 'Cpappd', 'Cpappv', 'Cpp', 'Cppa', 'Cppaa', 'Cppaaa', 'Cppaap', 'Cppap', 'Cppapa', 'Cppapp', 'Cppp', 'Cpppa', 'Cpppaa', 'Cpppap', 'Cpppp', 'Cppppa', 'Cppppp', 'D', 'Da', 'Daa', 'Daaa', 'Daap', 'Dap', 'Dapa', 'Dapp', 'Dp', 'Dpa', 'Dpaa', 'Dpap', 'Dpp', 'Dppa', 'Dppp', 'E', 'Ea', 'Eal', 'Eala', 'Ealaa', 'Ealap', 'Ealp', 'Ealpa', 'Ealpp', 'Ear', 'Eara', 'Earaa', 'Earap', 'Earp', 'Earpa', 'Earpp', 'Ep', 'Epl', 'Epla', 'Eplaa', 'Eplap', 'Eplp', 'Eplpa', 'Eplpp', 'Epr', 'Epra', 'Epraa', 'Eprap', 'Eprp', 'Eprpa', 'Eprpp', 'MS', 'MSa', 'MSaa', 'MSaaa', 'MSaaaa', 'MSaaaaa', 'MSaaaaaa', 'MSaaaaap', 'MSaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaap', 'MSaaapa', 'MSaaapaa', 'MSaaapap', 'MSaaapp', 'MSaap', 'MSaapa', 'MSaapaa', 'MSaapaaa', 'MSaapaap', 'MSaapap', 'MSaapapa', 'MSaapapp', 'MSaapp', 'MSaappa', 'MSaappp', 'MSap', 'MSapa', 'MSapaa', 'MSapaaa', 'MSapaaaa', 'MSapaaap', 'MSapaap', 'MSapap', 'MSapapa', 'MSapapaa', 'MSapapap', 'MSapapp', 'MSapappa', 'MSapappp', 'MSapp', 'MSappa', 'MSappaa', 'MSappap', 'MSappp', 'MSapppa', 'MSapppp', 'MSp', 'MSpa', 'MSpaa', 'MSpaaa', 'MSpaaaa', 'MSpaaaaa', 'MSpaaaap', 'MSpaaap', 'MSpaaapa', 'MSpaaapp', 'MSpaap', 'MSpaapa', 'MSpaapaa', 'MSpaapap', 'MSpaapp', 'MSpap', 'MSpapa', 'MSpapaa', 'MSpapaaa', 'MSpapaap', 'MSpapap', 'MSpapapa', 'MSpapapp', 'MSpapp', 'MSpappa', 'MSpappp', 'MSpp', 'MSppa', 'MSppaa', 'MSppaaa', 'MSppaaaa', 'MSppaaap', 'MSppaap', 'MSppap', 'MSppapa', 'MSppapaa', 'MSppapap', 'MSppapp', 'MSppappa', 'MSppappp', 'MSppp', 'MSpppa', 'MSpppaa', 'MSpppap', 'MSpppp', 'MSppppa', 'MSppppp', 'P3', 'P4', 'Z2', 'Z3']

a=0
embryo_name = []
embryo_length = []
for each_emb in raw_data_list:
    
    emb_output_ = pd.read_csv(each_emb, sep='\t')
    
    embryo_name_ = each_emb.split('_raw')[0]
    embryo_name.append(embryo_name_)

    emb_output_info = emb_output_.set_index(['cell_name'])
    
    emb_output_info = emb_output_info[['time','X','Y','Z']]
    
    emb_output_info['Z']=emb_output_info['Z']*4.5
    
    #################################################################################################################
    emb_output_temp=emb_output_info.drop(['time'],axis=1)
    cell_pair_distance_all_tp=scipy.spatial.distance.pdist(emb_output_temp)                
    cell_pair_distance_all_tp.sort()            
    dist_max=np.median(cell_pair_distance_all_tp[-50:])
    embryo_length.append(dist_max)                            
    
    cell_list=[]
    cell_distance=[]
    
    for cell in Cell_list_714:
        cell_list.append(cell)
        if cell in emb_output_info.index:
            if type(emb_output_info.loc[cell]['time'])==pd.core.series.Series:
                time_point=int(emb_output_info.loc[cell]['time'].median()+0.5)-(emb_output_info.loc[cell]['time'].min())#int()+0.5为了取偶数后面的点
                cell_info=emb_output_info.loc[cell].iloc[time_point]
                cell_position=cell_info.drop('time') 
            else:
                cell_position=emb_output_info.loc[cell].drop('time')
        else:
            cell_position=pd.DataFrame([np.nan,np.nan,np.nan],index=['X','Y','Z'],columns=[cell])
        if type(cell_distance)==list:
            cell_distance=pd.DataFrame(cell_position)
        else:
            cell_distance=pd.concat([cell_distance,pd.DataFrame(cell_position)],axis=1)
                    
    cell_distance=cell_distance.T
    
    cell_distance.to_csv(os.path.join(XYZ_path,embryo_name_+'_XYZ.txt'), sep="\t")

    
    cell_pair_distance=scipy.spatial.distance.pdist(cell_distance)
    #distance_exclude_nan = cell_pair_distance[~np.isnan(cell_pair_distance)]#保留不是nan的值
    #dist_max=distance_exclude_nan.max()#胚胎根据最长细胞间距离标准化
    cell_pair_distance=cell_pair_distance/dist_max
    cell_pair_distance=pd.DataFrame(scipy.spatial.distance.squareform(cell_pair_distance),
                                    index=cell_list,columns=cell_list)
    cell_pair_distance.to_csv(os.path.join(distance_matrix_path,embryo_name_+'_distMatrix_normalized.txt'), sep="\t")
    a+=1
    print(a)   
embryo_length=pd.DataFrame(embryo_length,columns=['embryo length'])
embryo_length.index=embryo_name
embryo_length.to_csv(os.path.join(embryo_length_path,embryo_length_file_name+'.txt'), sep='\t')