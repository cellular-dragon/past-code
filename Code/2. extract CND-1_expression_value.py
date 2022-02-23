# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 13:06:36 2022

@author: dell
"""

import os
import pandas as pd
import numpy as np

#output_path = os.path.abspath(r'C:\Users\dell\Desktop\output_folder')###目的总文件夹
#if not os.path.exists(output_path):
    #os.makedirs(output_path)

File_Directory = r'C:\Users\dell\Desktop\Raw_data_Control'   # input raw data of control or RNAi or RNAi_add
Output_Folder = r'C:\Users\dell\Desktop\output_folder'   # assign the output folder path
output_gene_expression_value_file_name='CND-1_expression_value_matrix' #input expression value file name
expression_cut_off=1500 #input expression cut off, for example "expression: >= 1500","no expression: < 1500"
output_gene_expression_value_binary_file_name=output_gene_expression_value_file_name+'_binary'

os.chdir(File_Directory)
raw_data_list = os.listdir(File_Directory)

a=0
embryo_name = []
output=[]
for each_emb in raw_data_list:
    
    emb_output_ = pd.read_csv(each_emb, sep='\t')
    
    embryo_name_ = each_emb.split('_raw')[0]
    embryo_name.append(embryo_name_)

    emb_output_['expression'] = emb_output_ ['raw-expression'] -  emb_output_['blot-correction']
    
    emb_output_ = emb_output_.set_index(['time', 'cell_name'])
    
    emb_output_ = emb_output_[['expression']]
    emb_output_mean = emb_output_.mean(level='cell_name')
            
    if type(output) == list:
        output = emb_output_mean
    else:
        output = pd.concat([output, emb_output_mean], axis=1)
    a+=1
    print(a)                
output.columns = embryo_name

Cell_list_714=['ABal', 'ABala', 'ABalaa', 'ABalaaa', 'ABalaaaa', 'ABalaaaal', 'ABalaaaala', 'ABalaaaalp', 'ABalaaaar', 'ABalaaaarl', 'ABalaaaarr', 'ABalaaap', 'ABalaaapa', 'ABalaaapal', 'ABalaaapar', 'ABalaaapp', 'ABalaaappl', 'ABalaaappr', 'ABalaap', 'ABalaapa', 'ABalaapaa', 'ABalaapaaa', 'ABalaapaap', 'ABalaapap', 'ABalaapapa', 'ABalaapapp', 'ABalaapp', 'ABalaappa', 'ABalaappaa', 'ABalaappap', 'ABalaappp', 'ABalaapppa', 'ABalaapppp', 'ABalap', 'ABalapa', 'ABalapaa', 'ABalapaaa', 'ABalapaaaa', 'ABalapaaap', 'ABalapaap', 'ABalapaapa', 'ABalapaapp', 'ABalapap', 'ABalapapa', 'ABalapapaa', 'ABalapapap', 'ABalapapp', 'ABalapappa', 'ABalapappp', 'ABalapp', 'ABalappa', 'ABalappaa', 'ABalappaaa', 'ABalappaap', 'ABalappap', 'ABalappapa', 'ABalappapp', 'ABalappp', 'ABalapppa', 'ABalapppaa', 'ABalapppap', 'ABalapppp', 'ABalappppa', 'ABalappppp', 'ABalp', 'ABalpa', 'ABalpaa', 'ABalpaaa', 'ABalpaaaa', 'ABalpaaaaa', 'ABalpaaaap', 'ABalpaaap', 'ABalpaaapa', 'ABalpaaapp', 'ABalpaap', 'ABalpaapa', 'ABalpaapaa', 'ABalpaapap', 'ABalpaapp', 'ABalpaappa', 'ABalpaappp', 'ABalpap', 'ABalpapa', 'ABalpapaa', 'ABalpapaaa', 'ABalpapaap', 'ABalpapap', 'ABalpapapa', 'ABalpapapp', 'ABalpapp', 'ABalpappa', 'ABalpappaa', 'ABalpappap', 'ABalpappp', 'ABalpapppa', 'ABalpapppp', 'ABalpp', 'ABalppa', 'ABalppaa', 'ABalppaaa', 'ABalppaaaa', 'ABalppaaap', 'ABalppaap', 'ABalppaapa', 'ABalppaapp', 'ABalppap', 'ABalppapa', 'ABalppapaa', 'ABalppapap', 'ABalppapp', 'ABalppappa', 'ABalppappp', 'ABalppp', 'ABalpppa', 'ABalpppaa', 'ABalpppaaa', 'ABalpppaap', 'ABalpppap', 'ABalpppapa', 'ABalpppapp', 'ABalpppp', 'ABalppppa', 'ABalppppaa', 'ABalppppap', 'ABalppppp', 'ABalpppppa', 'ABalpppppp', 'ABar', 'ABara', 'ABaraa', 'ABaraaa', 'ABaraaaa', 'ABaraaaaa', 'ABaraaaaaa', 'ABaraaaaap', 'ABaraaaap', 'ABaraaaapa', 'ABaraaaapp', 'ABaraaap', 'ABaraaapa', 'ABaraaapaa', 'ABaraaapap', 'ABaraaapp', 'ABaraaappa', 'ABaraaappp', 'ABaraap', 'ABaraapa', 'ABaraapaa', 'ABaraapaaa', 'ABaraapaap', 'ABaraapap', 'ABaraapapa', 'ABaraapapp', 'ABaraapp', 'ABaraappa', 'ABaraappaa', 'ABaraappap', 'ABaraappp', 'ABaraapppa', 'ABaraapppp', 'ABarap', 'ABarapa', 'ABarapaa', 'ABarapaaa', 'ABarapaaaa', 'ABarapaaap', 'ABarapaap', 'ABarapaapa', 'ABarapaapp', 'ABarapap', 'ABarapapa', 'ABarapapaa', 'ABarapapap', 'ABarapapp', 'ABarapappa', 'ABarapappp', 'ABarapp', 'ABarappa', 'ABarappaa', 'ABarappaaa', 'ABarappaap', 'ABarappap', 'ABarappapa', 'ABarappapp', 'ABarappp', 'ABarapppa', 'ABarapppaa', 'ABarapppap', 'ABarapppp', 'ABarappppa', 'ABarappppp', 'ABarp', 'ABarpa', 'ABarpaa', 'ABarpaaa', 'ABarpaaaa', 'ABarpaaaaa', 'ABarpaaaap', 'ABarpaaap', 'ABarpaaapa', 'ABarpaaapp', 'ABarpaap', 'ABarpaapa', 'ABarpaapaa', 'ABarpaapap', 'ABarpaapp', 'ABarpaappa', 'ABarpaappp', 'ABarpap', 'ABarpapa', 'ABarpapaa', 'ABarpapaaa', 'ABarpapaap', 'ABarpapap', 'ABarpapapa', 'ABarpapapp', 'ABarpapp', 'ABarpappa', 'ABarpappaa', 'ABarpappap', 'ABarpappp', 'ABarpapppa', 'ABarpapppp', 'ABarpp', 'ABarppa', 'ABarppaa', 'ABarppaaa', 'ABarppaaaa', 'ABarppaaap', 'ABarppaap', 'ABarppaapa', 'ABarppaapp', 'ABarppap', 'ABarppapa', 'ABarppapaa', 'ABarppapap', 'ABarppapp', 'ABarppappa', 'ABarppappp', 'ABarppp', 'ABarpppa', 'ABarpppaa', 'ABarpppaaa', 'ABarpppaap', 'ABarpppap', 'ABarpppapa', 'ABarpppapp', 'ABarpppp', 'ABarppppa', 'ABarppppaa', 'ABarppppap', 'ABarppppp', 'ABarpppppa', 'ABarpppppp', 'ABpl', 'ABpla', 'ABplaa', 'ABplaaa', 'ABplaaaa', 'ABplaaaaa', 'ABplaaaaaa', 'ABplaaaaap', 'ABplaaaap', 'ABplaaaapa', 'ABplaaaapp', 'ABplaaap', 'ABplaaapa', 'ABplaaapaa', 'ABplaaapap', 'ABplaaapp', 'ABplaaappa', 'ABplaaappp', 'ABplaap', 'ABplaapa', 'ABplaapaa', 'ABplaapaaa', 'ABplaapaap', 'ABplaapap', 'ABplaapapa', 'ABplaapapp', 'ABplaapp', 'ABplaappa', 'ABplaappaa', 'ABplaappap', 'ABplaappp', 'ABplaapppa', 'ABplaapppp', 'ABplap', 'ABplapa', 'ABplapaa', 'ABplapaaa', 'ABplapaaaa', 'ABplapaaap', 'ABplapaap', 'ABplapaapa', 'ABplapaapp', 'ABplapap', 'ABplapapa', 'ABplapapaa', 'ABplapapap', 'ABplapapp', 'ABplapappa', 'ABplapappp', 'ABplapp', 'ABplappa', 'ABplappaa', 'ABplappaaa', 'ABplappaap', 'ABplappap', 'ABplappapa', 'ABplappapp', 'ABplappp', 'ABplapppa', 'ABplapppaa', 'ABplapppap', 'ABplapppp', 'ABplappppa', 'ABplappppp', 'ABplp', 'ABplpa', 'ABplpaa', 'ABplpaaa', 'ABplpaaaa', 'ABplpaaaaa', 'ABplpaaaap', 'ABplpaaap', 'ABplpaaapa', 'ABplpaaapp', 'ABplpaap', 'ABplpaapa', 'ABplpaapaa', 'ABplpaapap', 'ABplpaapp', 'ABplpaappa', 'ABplpaappp', 'ABplpap', 'ABplpapa', 'ABplpapaa', 'ABplpapaaa', 'ABplpapaap', 'ABplpapap', 'ABplpapapa', 'ABplpapapp', 'ABplpapp', 'ABplpappa', 'ABplpappaa', 'ABplpappap', 'ABplpappp', 'ABplpapppa', 'ABplpapppp', 'ABplpp', 'ABplppa', 'ABplppaa', 'ABplppaaa', 'ABplppaaaa', 'ABplppaaap', 'ABplppaap', 'ABplppaapa', 'ABplppaapp', 'ABplppap', 'ABplppapa', 'ABplppapaa', 'ABplppapap', 'ABplppapp', 'ABplppappa', 'ABplppappp', 'ABplppp', 'ABplpppa', 'ABplpppaa', 'ABplpppaaa', 'ABplpppaap', 'ABplpppap', 'ABplpppapa', 'ABplpppapp', 'ABplpppp', 'ABplppppa', 'ABplppppaa', 'ABplppppap', 'ABplppppp', 'ABplpppppa', 'ABplpppppp', 'ABpr', 'ABpra', 'ABpraa', 'ABpraaa', 'ABpraaaa', 'ABpraaaaa', 'ABpraaaaaa', 'ABpraaaaap', 'ABpraaaap', 'ABpraaaapa', 'ABpraaaapp', 'ABpraaap', 'ABpraaapa', 'ABpraaapaa', 'ABpraaapap', 'ABpraaapp', 'ABpraaappa', 'ABpraaappp', 'ABpraap', 'ABpraapa', 'ABpraapaa', 'ABpraapaaa', 'ABpraapaap', 'ABpraapap', 'ABpraapapa', 'ABpraapapp', 'ABpraapp', 'ABpraappa', 'ABpraappaa', 'ABpraappap', 'ABpraappp', 'ABpraapppa', 'ABpraapppp', 'ABprap', 'ABprapa', 'ABprapaa', 'ABprapaaa', 'ABprapaaaa', 'ABprapaaap', 'ABprapaap', 'ABprapaapa', 'ABprapaapp', 'ABprapap', 'ABprapapa', 'ABprapapaa', 'ABprapapap', 'ABprapapp', 'ABprapappa', 'ABprapappp', 'ABprapp', 'ABprappa', 'ABprappaa', 'ABprappaaa', 'ABprappaap', 'ABprappap', 'ABprappapa', 'ABprappapp', 'ABprappp', 'ABprapppa', 'ABprapppaa', 'ABprapppap', 'ABprapppp', 'ABprappppa', 'ABprappppp', 'ABprp', 'ABprpa', 'ABprpaa', 'ABprpaaa', 'ABprpaaaa', 'ABprpaaaaa', 'ABprpaaaap', 'ABprpaaap', 'ABprpaaapa', 'ABprpaaapp', 'ABprpaap', 'ABprpaapa', 'ABprpaapaa', 'ABprpaapap', 'ABprpaapp', 'ABprpaappa', 'ABprpaappp', 'ABprpap', 'ABprpapa', 'ABprpapaa', 'ABprpapaaa', 'ABprpapaap', 'ABprpapap', 'ABprpapapa', 'ABprpapapp', 'ABprpapp', 'ABprpappa', 'ABprpappaa', 'ABprpappap', 'ABprpappp', 'ABprpapppa', 'ABprpapppp', 'ABprpp', 'ABprppa', 'ABprppaa', 'ABprppaaa', 'ABprppaaaa', 'ABprppaaap', 'ABprppaap', 'ABprppaapa', 'ABprppaapp', 'ABprppap', 'ABprppapa', 'ABprppapaa', 'ABprppapap', 'ABprppapp', 'ABprppappa', 'ABprppappp', 'ABprppp', 'ABprpppa', 'ABprpppaa', 'ABprpppaaa', 'ABprpppaap', 'ABprpppap', 'ABprpppapa', 'ABprpppapp', 'ABprpppp', 'ABprppppa', 'ABprppppaa', 'ABprppppap', 'ABprppppp', 'ABprpppppa', 'ABprpppppp', 'C', 'Ca', 'Caa', 'Caaa', 'Caaaa', 'Caaaaa', 'Caaaap', 'Caaap', 'Caaapa', 'Caaapp', 'Caap', 'Caapa', 'Caapp', 'Caappd', 'Caappv', 'Cap', 'Capa', 'Capaa', 'Capaaa', 'Capaap', 'Capap', 'Capapa', 'Capapp', 'Capp', 'Cappa', 'Cappaa', 'Cappap', 'Cappp', 'Capppa', 'Capppp', 'Cp', 'Cpa', 'Cpaa', 'Cpaaa', 'Cpaaaa', 'Cpaaap', 'Cpaap', 'Cpaapa', 'Cpaapp', 'Cpap', 'Cpapa', 'Cpapaa', 'Cpapap', 'Cpapp', 'Cpappd', 'Cpappv', 'Cpp', 'Cppa', 'Cppaa', 'Cppaaa', 'Cppaap', 'Cppap', 'Cppapa', 'Cppapp', 'Cppp', 'Cpppa', 'Cpppaa', 'Cpppap', 'Cpppp', 'Cppppa', 'Cppppp', 'D', 'Da', 'Daa', 'Daaa', 'Daap', 'Dap', 'Dapa', 'Dapp', 'Dp', 'Dpa', 'Dpaa', 'Dpap', 'Dpp', 'Dppa', 'Dppp', 'E', 'Ea', 'Eal', 'Eala', 'Ealaa', 'Ealap', 'Ealp', 'Ealpa', 'Ealpp', 'Ear', 'Eara', 'Earaa', 'Earap', 'Earp', 'Earpa', 'Earpp', 'Ep', 'Epl', 'Epla', 'Eplaa', 'Eplap', 'Eplp', 'Eplpa', 'Eplpp', 'Epr', 'Epra', 'Epraa', 'Eprap', 'Eprp', 'Eprpa', 'Eprpp', 'MS', 'MSa', 'MSaa', 'MSaaa', 'MSaaaa', 'MSaaaaa', 'MSaaaaaa', 'MSaaaaap', 'MSaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaap', 'MSaaapa', 'MSaaapaa', 'MSaaapap', 'MSaaapp', 'MSaap', 'MSaapa', 'MSaapaa', 'MSaapaaa', 'MSaapaap', 'MSaapap', 'MSaapapa', 'MSaapapp', 'MSaapp', 'MSaappa', 'MSaappp', 'MSap', 'MSapa', 'MSapaa', 'MSapaaa', 'MSapaaaa', 'MSapaaap', 'MSapaap', 'MSapap', 'MSapapa', 'MSapapaa', 'MSapapap', 'MSapapp', 'MSapappa', 'MSapappp', 'MSapp', 'MSappa', 'MSappaa', 'MSappap', 'MSappp', 'MSapppa', 'MSapppp', 'MSp', 'MSpa', 'MSpaa', 'MSpaaa', 'MSpaaaa', 'MSpaaaaa', 'MSpaaaap', 'MSpaaap', 'MSpaaapa', 'MSpaaapp', 'MSpaap', 'MSpaapa', 'MSpaapaa', 'MSpaapap', 'MSpaapp', 'MSpap', 'MSpapa', 'MSpapaa', 'MSpapaaa', 'MSpapaap', 'MSpapap', 'MSpapapa', 'MSpapapp', 'MSpapp', 'MSpappa', 'MSpappp', 'MSpp', 'MSppa', 'MSppaa', 'MSppaaa', 'MSppaaaa', 'MSppaaap', 'MSppaap', 'MSppap', 'MSppapa', 'MSppapaa', 'MSppapap', 'MSppapp', 'MSppappa', 'MSppappp', 'MSppp', 'MSpppa', 'MSpppaa', 'MSpppap', 'MSpppp', 'MSppppa', 'MSppppp', 'P3', 'P4', 'Z2', 'Z3']

output=output.loc[Cell_list_714]

mask_exp = output <= 0
output[mask_exp] = 0

output.to_csv(os.path.join(Output_Folder,output_gene_expression_value_file_name+'.txt'), sep='\t')

mask_exp = output >= 1500
mask_non_exp = output < 1500

output[mask_exp] = 1
output[mask_non_exp] = 0

output.to_csv(os.path.join(Output_Folder,output_gene_expression_value_binary_file_name+'.txt'), sep='\t')

