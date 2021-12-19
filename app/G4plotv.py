"""
Created on Tue Oct 12 21:58:17 2021
@author: priyasharma16
"""
import matplotlib.pyplot as plt
import pandas as pd
import os, json

def drawG4Plot(folderPath, seq_len_list, seq_1_dict, seq_2_dict):
    ID_1, pos_lis, length_lis, content_lis = getPositionsData(seq_1_dict)
    ID_2, pos_lis_2, length_lis_2, content_lis_2 = getPositionsData(seq_2_dict)

    unique_content_lis, unique_pos_lis, unique_length_lis = getUniqueData(content_lis, content_lis_2, pos_lis, length_lis, ID_1, folderPath)
    unique_content_lis2, unique_pos_lis2, unique_length_lis2 = getUniqueData(content_lis_2, content_lis, pos_lis_2, length_lis_2,ID_2, folderPath)
    
    getCommonData(content_lis, content_lis_2, pos_lis, ID_1, folderPath)
    getCommonData(content_lis_2, content_lis, pos_lis_2, ID_2, folderPath)
    if unique_content_lis and unique_content_lis2:
        pos_lis = unique_pos_lis
        length_lis = unique_length_lis
        pos_lis2 = unique_pos_lis2
        length_lis2 = unique_length_lis2
                        
        xmin = min(pos_lis)
        xmax = max(pos_lis)
        xmin2 = min(pos_lis2)
        xmax2 = max(pos_lis2)
                        
        df_1 = pd.DataFrame({'positions': pos_lis,'lengths': length_lis})
        df_2 = pd.DataFrame({'positions': pos_lis2,'lengths': length_lis2})
                        
        fig, ax = plt.subplots(figsize=(16,10), dpi= 80)
        #for plot 1      
        ax.hlines(0, xmin, xmax, color ="black")
        ax.vlines(x=df_1.positions, ymin=0, ymax=df_1.lengths, color='red', alpha=0.7, linewidth=1,label = f'Quadruplex in {ID_1}')
        ax.scatter(x=df_1.positions, y=df_1.lengths, s=50, color='red', alpha=0.7)
        for i in pos_lis:
            ax.text(i,0, s = i, fontsize=7, fontweight = 'bold')
        text_lis = [str(tex).lstrip("-") for tex in length_lis]
        df_1['text'] = text_lis
                    
        for row in df_1.itertuples():
            ax.text(row.positions, row.lengths+.5, s=(int(row.text), row.positions), horizontalalignment= 'center', rotation = 'vertical', fontweight = 'bold', verticalalignment='bottom', fontsize=7)
            
        ax.hlines(0, xmin2, xmax2, color ="black")
        ax.vlines(x=df_2.positions, ymin=0, ymax=df_2.lengths, color='blue', alpha=0.7, linewidth=1, label = f'Quadruplex in {ID_2}')
        ax.scatter(x=df_2.positions, y=df_2.lengths, s=50, color='blue', alpha=0.7)
        for j in pos_lis2:
            ax.text(j,0, s = j, fontsize=7, fontweight = 'bold')
    
        ax.set_xticks([])
        ax.set_yticks([])
        plt.ylabel("Length of Quadruplexes")
        plt.xlabel("Positions of Quadruplexes")
        text_lis2 = [str(tex).lstrip("-") for tex in length_lis2]
        df_2['text'] = text_lis2
        plt.box(False)
                    
        for row in df_2.itertuples():
            ax.text(row.positions, row.lengths+.5, s=(int(row.text), row.positions), horizontalalignment= 'center', rotation = 'vertical', fontweight = 'bold', verticalalignment='bottom', fontsize=7)
        ax.legend(prop = {'style': 'italic'}, loc='upper left', bbox_to_anchor=(1, 0.5))
        plt.savefig(fname = os.path.join(folderPath, 'plot_over.png'), pad_inches = 0.1, bbox_inches = "tight")
        plt.close()           

def getPositionsData(seq_dict):
    for ID, coordinates in seq_dict.items():
        pos_lis = []
        length_lis = []
        content_lis = []
        for values in coordinates:            
            pos_lis.append(int(values.split('\t')[0]))
            length_lis.append(int(values.split('\t')[1]))
            content_lis.append(str(values.split('\t')[2]))
    return (ID, pos_lis, length_lis, content_lis)

def getUniqueData(content_list_1, content_list_2, pos_lis, length_lis, query_id, folderPath):
    unique_content_lis = list(set(content_list_1) - set(content_list_2))
    unique_pos_lis = []
    unique_length_lis = []
    for i in unique_content_lis:
        unique_pos_lis.append(pos_lis[content_list_1.index(i)])
        unique_length_lis.append(length_lis[content_list_1.index(i)])

    writable_dict = motif_file_header()
    position_dict = {}
    for ele in unique_content_lis:
        position_dict[unique_pos_lis[ unique_content_lis.index(ele)]] = [ele]
    if position_dict:
        position_dict = { key: position_dict[key] for key in sorted(position_dict)}
        writable_dict.update({ query_id : position_dict})
        file_name_write = os.path.join(folderPath, 'output_' + query_id + '_unique_data.json')
        with open(file_name_write, 'w') as fp_write:
            json.dump(writable_dict, fp_write)
    return (unique_content_lis, unique_pos_lis, unique_length_lis)

def getCommonData(content_list_1, content_list_2, pos_lis, query_id, folderPath):
    common_from_content_lis = list(set(content_list_1) & set(content_list_2))
    common_pos_lis_1=[]
    for i in common_from_content_lis:
        common_pos_lis_1.append(pos_lis[content_list_1.index(i)])
    
    writable_dict = motif_file_header()
    position_dict = {}
    for ele in common_from_content_lis:
        position_dict[common_pos_lis_1[ common_from_content_lis.index(ele)]] = [ele]
    if position_dict:
        position_dict = { key: position_dict[key] for key in sorted(position_dict)}
        writable_dict.update({ query_id : position_dict})
        file_name_write = os.path.join(folderPath, 'output_' + query_id + '_common_data.json')
        with open(file_name_write, 'w') as fp_write:
            json.dump(writable_dict, fp_write)

def motif_file_header():
    return { 
        "header": ["Sequence ID", "G4 sequence", "Length", "Start Position", "End Position"]
    }