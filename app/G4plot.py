"""
Created on Tue Oct 12 21:58:17 2021
@author: priyasharma16
"""
import matplotlib.pyplot as plt
import os

def drawG4Plot(folderPath, seq_len_list, positiveDict, negativeDict):
    main_dict = {}
    for keys,values in positiveDict.items():
        for key, value in values.items():
            for item in value:
                data_list = main_dict.get(keys, [])
                data_list.append(f'{key}\t{len(item)}')
                main_dict[keys] = data_list

    for keys,values in negativeDict.items():
        for key, value in values.items():
            for item in value:
                data_list = main_dict.get(keys, [])
                data_list.append(f'{key}\t{-len(item)}')
                main_dict[keys] = data_list

    for ID, coordinates in main_dict.items():
        pos_lis = []
        length_lis = []
    
        for values in coordinates:            
            pos_lis.append(int(values.split('\t')[0]))
            length_lis.append(int(values.split('\t')[1]))

        text_lis = [str(tex).lstrip("-") for tex in length_lis]
        tex_lis = [int(x) for x in text_lis]
        xmin = min(pos_lis)
        xmax = max(pos_lis)
        
        fig, ax = plt.subplots(figsize=(16,10), dpi= 80)      
        ax.hlines(0, xmin, xmax, color ="black")
        ax.vlines(x=pos_lis, ymin=0, ymax=length_lis, color=['red' if x < 0 else 'blue' for x in length_lis], alpha=0.7, linewidth=0.8)
        ax.scatter(x=pos_lis, y=length_lis, s=25, color=['red' if x < 0 else 'blue' for x in length_lis], alpha=0.7)
        ax.set_xticks([])
        ax.set_yticks([]) 
        plt.ylabel("Length of Quadruplexes")
        plt.xlabel("Positions of Quadruplexes")
        plt.box(False)
        for x,y,tex in zip(pos_lis, length_lis, tex_lis):
            ax.text(x, y+.5, s=(tex, x), horizontalalignment= 'right' if y < 0 else 'center', rotation = 'vertical', fontweight = 'bold', verticalalignment='bottom', fontsize=7)
            
        #ax.legend(prop = {'style': 'italic'}, loc='upper left', bbox_to_anchor=(1, 0.5))
        png_file_path = os.path.join(folderPath, f'{ID}.png')
        plt.savefig(fname = png_file_path, pad_inches = 0.1, format='png', bbox_inches = "tight", dpi=200)   
        plt.clf()   
        plt.close()
