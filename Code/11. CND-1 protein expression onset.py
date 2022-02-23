# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 17:10:59 2022

@author: dell
"""

import os
import numpy as np
import pandas as pd
import math


##############################functions
def creat_track(cell_name):
    track=[]
    for j in Cell_List_AB:
        if j in cell_name:
            track.append(j)
    return track

def creat_sublineage(ternimal_cell):
    sub=ternimal_cell[0:4]
    sublineage=[]
    for cell_name in Cell_List_AB:        
        if cell_name[0:4] == sub and len(cell_name)<=len(ternimal_cell):
            sublineage.append(cell_name)
    return sublineage
##############################
Expressed_terminal_cells=['ABalpapaa','ABalpapap','ABalppapa','ABalppapp','ABarappa','ABarappp','ABplaapaa','ABplaapap','ABplpaaa','ABplpaap','ABplppapa','ABplppapp','ABpraapaa','ABpraapap','ABprpaaa','ABprpaap','ABprppapa','ABprppapp']
Cell_List_AB = ['ABal','ABala','ABalaa','ABalaaa','ABalaaaa','ABalaaaal','ABalaaaala','ABalaaaalp','ABalaaaar','ABalaaaarl','ABalaaaarr','ABalaaap','ABalaaapa','ABalaaapal','ABalaaapar','ABalaaapp','ABalaaappl','ABalaaappr','ABalaap','ABalaapa','ABalaapaa','ABalaapaaa','ABalaapaap','ABalaapap','ABalaapapa','ABalaapapp','ABalaapp','ABalaappa','ABalaappaa','ABalaappap','ABalaappp','ABalaapppa','ABalaapppp','ABalap','ABalapa','ABalapaa','ABalapaaa','ABalapaaaa','ABalapaaap','ABalapaap','ABalapaapa','ABalapaapp','ABalapap','ABalapapa','ABalapapaa','ABalapapap','ABalapapp','ABalapappa','ABalapappp','ABalapp','ABalappa','ABalappaa','ABalappaaa','ABalappaap','ABalappap','ABalappapa','ABalappapp','ABalappp','ABalapppa','ABalapppaa','ABalapppap','ABalapppp','ABalappppa','ABalappppp','ABalp','ABalpa','ABalpaa','ABalpaaa','ABalpaaaa','ABalpaaaaa','ABalpaaaap','ABalpaaap','ABalpaaapa','ABalpaaapp','ABalpaap','ABalpaapa','ABalpaapaa','ABalpaapap','ABalpaapp','ABalpaappa','ABalpaappp','ABalpap','ABalpapa','ABalpapaa','ABalpapaaa','ABalpapaap','ABalpapap','ABalpapapa','ABalpapapp','ABalpapp','ABalpappa','ABalpappaa','ABalpappap','ABalpappp','ABalpapppa','ABalpapppp','ABalpp','ABalppa','ABalppaa','ABalppaaa','ABalppaaaa','ABalppaaap','ABalppaap','ABalppaapa','ABalppaapp','ABalppap','ABalppapa','ABalppapaa','ABalppapap','ABalppapp','ABalppappa','ABalppappp','ABalppp','ABalpppa','ABalpppaa','ABalpppaaa','ABalpppaap','ABalpppap','ABalpppapa','ABalpppapp','ABalpppp','ABalppppa','ABalppppaa','ABalppppap','ABalppppp','ABalpppppa','ABalpppppp','ABar','ABara','ABaraa','ABaraaa','ABaraaaa','ABaraaaaa','ABaraaaaaa','ABaraaaaap','ABaraaaap','ABaraaaapa','ABaraaaapp','ABaraaap','ABaraaapa','ABaraaapaa','ABaraaapap','ABaraaapp','ABaraaappa','ABaraaappp','ABaraap','ABaraapa','ABaraapaa','ABaraapaaa','ABaraapaap','ABaraapap','ABaraapapa','ABaraapapp','ABaraapp','ABaraappa','ABaraappaa','ABaraappap','ABaraappp','ABaraapppa','ABaraapppp','ABarap','ABarapa','ABarapaa','ABarapaaa','ABarapaaaa','ABarapaaap','ABarapaap','ABarapaapa','ABarapaapp','ABarapap','ABarapapa','ABarapapaa','ABarapapap','ABarapapp','ABarapappa','ABarapappp','ABarapp','ABarappa','ABarappaa','ABarappaaa','ABarappaap','ABarappap','ABarappapa','ABarappapp','ABarappp','ABarapppa','ABarapppaa','ABarapppap','ABarapppp','ABarappppa','ABarappppp','ABarp','ABarpa','ABarpaa','ABarpaaa','ABarpaaaa','ABarpaaaaa','ABarpaaaap','ABarpaaap','ABarpaaapa','ABarpaaapp','ABarpaap','ABarpaapa','ABarpaapaa','ABarpaapap','ABarpaapp','ABarpaappa','ABarpaappp','ABarpap','ABarpapa','ABarpapaa','ABarpapaaa','ABarpapaap','ABarpapap','ABarpapapa','ABarpapapp','ABarpapp','ABarpappa','ABarpappaa','ABarpappap','ABarpappp','ABarpapppa','ABarpapppp','ABarpp','ABarppa','ABarppaa','ABarppaaa','ABarppaaaa','ABarppaaap','ABarppaap','ABarppaapa','ABarppaapp','ABarppap','ABarppapa','ABarppapaa','ABarppapap','ABarppapp','ABarppappa','ABarppappp','ABarppp','ABarpppa','ABarpppaa','ABarpppaaa','ABarpppaap','ABarpppap','ABarpppapa','ABarpppapp','ABarpppp','ABarppppa','ABarppppaa','ABarppppap','ABarppppp','ABarpppppa','ABarpppppp','ABpl','ABpla','ABplaa','ABplaaa','ABplaaaa','ABplaaaaa','ABplaaaaaa','ABplaaaaap','ABplaaaap','ABplaaaapa','ABplaaaapp','ABplaaap','ABplaaapa','ABplaaapaa','ABplaaapap','ABplaaapp','ABplaaappa','ABplaaappp','ABplaap','ABplaapa','ABplaapaa','ABplaapaaa','ABplaapaap','ABplaapap','ABplaapapa','ABplaapapp','ABplaapp','ABplaappa','ABplaappaa','ABplaappap','ABplaappp','ABplaapppa','ABplaapppp','ABplap','ABplapa','ABplapaa','ABplapaaa','ABplapaaaa','ABplapaaap','ABplapaap','ABplapaapa','ABplapaapp','ABplapap','ABplapapa','ABplapapaa','ABplapapap','ABplapapp','ABplapappa','ABplapappp','ABplapp','ABplappa','ABplappaa','ABplappaaa','ABplappaap','ABplappap','ABplappapa','ABplappapp','ABplappp','ABplapppa','ABplapppaa','ABplapppap','ABplapppp','ABplappppa','ABplappppp','ABplp','ABplpa','ABplpaa','ABplpaaa','ABplpaaaa','ABplpaaaaa','ABplpaaaap','ABplpaaap','ABplpaaapa','ABplpaaapp','ABplpaap','ABplpaapa','ABplpaapaa','ABplpaapap','ABplpaapp','ABplpaappa','ABplpaappp','ABplpap','ABplpapa','ABplpapaa','ABplpapaaa','ABplpapaap','ABplpapap','ABplpapapa','ABplpapapp','ABplpapp','ABplpappa','ABplpappaa','ABplpappap','ABplpappp','ABplpapppa','ABplpapppp','ABplpp','ABplppa','ABplppaa','ABplppaaa','ABplppaaaa','ABplppaaap','ABplppaap','ABplppaapa','ABplppaapp','ABplppap','ABplppapa','ABplppapaa','ABplppapap','ABplppapp','ABplppappa','ABplppappp','ABplppp','ABplpppa','ABplpppaa','ABplpppaaa','ABplpppaap','ABplpppap','ABplpppapa','ABplpppapp','ABplpppp','ABplppppa','ABplppppaa','ABplppppap','ABplppppp','ABplpppppa','ABplpppppp','ABpr','ABpra','ABpraa','ABpraaa','ABpraaaa','ABpraaaaa','ABpraaaaaa','ABpraaaaap','ABpraaaap','ABpraaaapa','ABpraaaapp','ABpraaap','ABpraaapa','ABpraaapaa','ABpraaapap','ABpraaapp','ABpraaappa','ABpraaappp','ABpraap','ABpraapa','ABpraapaa','ABpraapaaa','ABpraapaap','ABpraapap','ABpraapapa','ABpraapapp','ABpraapp','ABpraappa','ABpraappaa','ABpraappap','ABpraappp','ABpraapppa','ABpraapppp','ABprap','ABprapa','ABprapaa','ABprapaaa','ABprapaaaa','ABprapaaap','ABprapaap','ABprapaapa','ABprapaapp','ABprapap','ABprapapa','ABprapapaa','ABprapapap','ABprapapp','ABprapappa','ABprapappp','ABprapp','ABprappa','ABprappaa','ABprappaaa','ABprappaap','ABprappap','ABprappapa','ABprappapp','ABprappp','ABprapppa','ABprapppaa','ABprapppap','ABprapppp','ABprappppa','ABprappppp','ABprp','ABprpa','ABprpaa','ABprpaaa','ABprpaaaa','ABprpaaaaa','ABprpaaaap','ABprpaaap','ABprpaaapa','ABprpaaapp','ABprpaap','ABprpaapa','ABprpaapaa','ABprpaapap','ABprpaapp','ABprpaappa','ABprpaappp','ABprpap','ABprpapa','ABprpapaa','ABprpapaaa','ABprpapaap','ABprpapap','ABprpapapa','ABprpapapp','ABprpapp','ABprpappa','ABprpappaa','ABprpappap','ABprpappp','ABprpapppa','ABprpapppp','ABprpp','ABprppa','ABprppaa','ABprppaaa','ABprppaaaa','ABprppaaap','ABprppaap','ABprppaapa','ABprppaapp','ABprppap','ABprppapa','ABprppapaa','ABprppapap','ABprppapp','ABprppappa','ABprppappp','ABprppp','ABprpppa','ABprpppaa','ABprpppaaa','ABprpppaap','ABprpppap','ABprpppapa','ABprpppapp','ABprpppp','ABprppppa','ABprppppaa','ABprppppap','ABprppppp','ABprpppppa','ABprpppppp']

all_track=[]
for i in Expressed_terminal_cells:
    track=[]
    for j in Cell_List_AB:
        if j in i:
            track.append(j)
    if type(all_track)==list:
        all_track=pd.DataFrame(track)
    else:
        all_track=pd.concat([all_track,pd.DataFrame(track)],axis=1)        
all_track.columns=Expressed_terminal_cells

all_sublineage=[]
for ternimal_cell in all_track.columns:
    sub=ternimal_cell[0:4]
    sublineage=[]
    for cell_name in Cell_List_AB:        
        if cell_name[0:4] == sub:
            sublineage.append(cell_name)
    if type(all_sublineage)==list:
        all_sublineage=pd.DataFrame(sublineage)
    else:
        all_sublineage=pd.concat([all_sublineage,pd.DataFrame(sublineage)],axis=1)        
all_sublineage.columns=all_track.columns
###############################

File_Directory = r'C:\Users\dell\Desktop\raw_data'   # input raw data of control or RNAi or RNAi_add
# input the corresponding CND-1 binary expression matrix of control or RNAi or RNAi_add from script "3. extract relative cell position distance matrix"
exp_matrix_binary = pd.read_csv(r"C:\Users\dell\Desktop\output_folder\CND-1_expression_value_matrix_binary.txt", index_col= 0,header=0,sep='\t')

Output_Folder = r'C:\Users\dell\Desktop\output_folder'   # assign the output folder path
output_file_name='CND-1_expression_onset_matrix' #input expression value file name

os.chdir(File_Directory)
raw_data_list = os.listdir(File_Directory)

lineage_start_time=[]
a=0
for terminal_cell_name in Expressed_terminal_cells:
    terminal_cell_name='ABplaapaa'
    track_selected=list(all_track[terminal_cell_name].dropna())
     
    embryo_name = []
    start_time=[]
    
    for each_emb in raw_data_list:
        each_emb='ctr_emb1_raw_data.txt'
        emb_output_ = pd.read_csv(each_emb, sep='\t')
        
        embryo_name_ = each_emb.split('_raw')[0]
        embryo_name.append(embryo_name_)
    
        emb_output_['expression'] = emb_output_ ['raw-expression'] -  emb_output_['blot-correction']
                        
        emb_output_ = emb_output_.set_index(['cell_name'])
        
        emb_output_ = emb_output_[['time','expression']]
        
    #################################################################################################根据表达矩阵找到最开始表达的细胞
        embryo_cell_exp=exp_matrix_binary[embryo_name_]                 
        sub_=embryo_cell_exp.loc[track_selected]                
        sub_=sub_.dropna()
        
        if sub_.max()==1:
            first_exp_cell=sub_[sub_==1].index[0]    
                        
            track_temp=creat_track(first_exp_cell)
            #sulineage_temp=creat_sublineage(first_exp_cell)
    ################################################################################################产生cell track校正的表达值信息                                             
            
            emb_output_2=emb_output_.loc[track_temp]#主要track的细胞列表
            time_point_min=emb_output_2['time'].min()
            emb_output_2['time']=emb_output_2['time']-time_point_min
            expression_value=emb_output_2.loc[first_exp_cell]['expression']
            expression_start=np.nan
            times=int(len(expression_value)-2)
            for i in range(times):                            
                exp_value_3=expression_value.iloc[i:i+3]
                if exp_value_3.iloc[0] >=1500 and exp_value_3.iloc[1] >=1500 and exp_value_3.iloc[2] >=1500:
                    expression_start=i
                    break
            if math.isnan(expression_start):
                start_time.append(expression_start)
            else:
                time_point=emb_output_2.loc[first_exp_cell]['time']
                expression_start_time=time_point.iloc[expression_start]                        
                start_time.append(expression_start_time)
        else:
            start_time.append(np.nan)

    if type(lineage_start_time)==list:
        lineage_start_time=pd.DataFrame(start_time)
    else:
        lineage_start_time=pd.concat([lineage_start_time,pd.DataFrame(start_time)],axis=1) 
    a+=1
    print(a)
    
lineage_start_time.index=embryo_name
lineage_start_time.columns=Expressed_terminal_cells
lineage_start_time.to_csv(os.path.join(Output_Folder,output_file_name+'.txt'), sep="\t")                













