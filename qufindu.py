"""
Created on Sun Sep 12 20:20:43 2021
@author: priyasharma16
"""
import re
import json, os, shutil
import G4plot as plt

class qufindu():
    def __init__(self, outputfolderPath,positive_dict,negative_dict, search_model, nucleobase, min_stem_size, max_stem_size, min_loop, max_loop, isPositiveChecked, isNegativeChecked, bulge_size=0, mismatch_num=0, CpGisland=False, cpg_model='nonoverlapping'):
        self.outputfolderPath = outputfolderPath
        self.nucleobase = nucleobase
        self.search_model = search_model
        self.cpg_model = cpg_model
        self.min_stem = min_stem_size
        self.max_stem = max_stem_size
        self.loopmin = min_loop
        self.loopmax = max_loop
        self.CpG = CpGisland
        self.bulge_size = bulge_size
        self.mismatch_num = mismatch_num
        self.plus_strands = positive_dict
        self.minus_strands  = negative_dict
        self.base, self.inv_base = self.make_base()
        reg_ex = self.regex_func() 

        self.positive_motif_dict, self.negative_motif_dict = {}, {}
        if isPositiveChecked:
            file_name = 'output_positive'
            CpG_dict, quad_dict = {},{}
            if self.CpG:
                CpG_dict, quad_dict = self.CpG_island(self.plus_strands, reg_ex)
            else:
                quad_dict = self.quadruplex_finder(self.plus_strands, reg_ex)
            
            CpG_dict = self.processDictionaryData(CpG_dict)
            quad_dict = self.processDictionaryData(quad_dict)
            
            # if CpG_dict and quad_dict:
            #     generatePositionFigure(CpG_dict, quad_dict, self.outputfolderPath, 'positive')
            
            if CpG_dict:
                CpG_dict = self.sortResult(CpG_dict)
                writable_dict = cpg_file_header()
                writable_dict.update(CpG_dict)
                file_name_write = os.path.join(self.outputfolderPath, file_name + '_cpg.json')
                with open(file_name_write, 'w') as fp_write:
                    json.dump(writable_dict, fp_write)

            if quad_dict:
                quad_dict = self.sortResult(quad_dict)
                self.positive_motif_dict = quad_dict
                writable_dict = motif_file_header()
                writable_dict.update(quad_dict)
                file_name_write = os.path.join(self.outputfolderPath, file_name + '_motif.json')
                with open(file_name_write, 'w') as fp_write:
                    json.dump(writable_dict, fp_write)

        if isNegativeChecked:
            file_name = 'output_negative'
            CpG_dict, quad_dict = {}, {}
            if self.CpG:
               CpG_dict, quad_dict = self.CpG_island(self.minus_strands, reg_ex)
            else:
                quad_dict = self.quadruplex_finder(self.minus_strands, reg_ex) 
            
            CpG_dict = self.processDictionaryData(CpG_dict)
            quad_dict = self.processDictionaryData(quad_dict)

            # if CpG_dict and quad_dict:
            #     generatePositionFigure(CpG_dict, quad_dict, self.outputfolderPath, 'negative')

            if CpG_dict:
                CpG_dict = self.sortResult(CpG_dict)
                writable_dict = cpg_file_header()
                writable_dict.update(CpG_dict)
                file_name_write = os.path.join(self.outputfolderPath, file_name + '_cpg.json')
                with open(file_name_write, 'w') as fp_write:
                    json.dump(writable_dict, fp_write)
            if quad_dict:
                quad_dict = self.sortResult(quad_dict)
                self.negative_motif_dict = quad_dict
                writable_dict = motif_file_header()
                writable_dict.update(quad_dict)
                file_name_write = os.path.join(self.outputfolderPath, file_name + '_motif.json')
                with open(file_name_write, 'w') as fp_write:
                    json.dump(writable_dict, fp_write)
        
        if self.positive_motif_dict or self.negative_motif_dict:
            if self.positive_motif_dict:
                seq_len_list = [len(value) for value in positive_dict.values()] 
            else:
                seq_len_list = [len(value) for value in negative_dict.values()]
            plt.drawG4Plot(self.outputfolderPath, seq_len_list, self.positive_motif_dict, self.negative_motif_dict)

    def processDictionaryData(self, input_dict):
        output_dict = {}
        for key, value in input_dict.items():
            if value:
                output_dict[key] = value
        return output_dict

    def sortResult(self, input_dict):
        final_dict = {}
        for key, value in input_dict.items():
            final_dict[ key ] = { i: value[i] for i in sorted(value.keys())}
        return final_dict
    
    def make_base(self):
        if self.nucleobase == 'G':
            invert_base = 'C'
        elif self.nucleobase == 'A':
            invert_base = 'T'
        elif self.nucleobase == 'T':
            invert_base = 'A'
        else:
            invert_base = 'G'
        return self.nucleobase, invert_base
    
    def regex_func(self):
        if self.bulge_size > 0 and self.mismatch_num == 0 and (self.base == 'G' or self.base == 'C'):
            
            regx = '((({b}{{1,}}[AT{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[AT{i}U]{{{bs}}}{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{3}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{1,}}[AT{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[AT{i}U]{{{bs}}}{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{2}}({b}{{1,}}[AT{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[AT{i}U]{{{bs}}}{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}({b}{{1,}}[AT{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[AT{i}U]{{{bs}}}{b}{{1,}})))'.format(b = self.base, i = self.inv_base,  bs = self.bulge_size, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        elif self.bulge_size > 0 and self.mismatch_num == 0 and (self.base == 'A' or self.base == 'T'):
            
            regx = '((({b}{{1,}}[GC{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[GC{i}U]{{{bs}}}{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{3}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{1,}}[GC{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[GC{i}U]{{{bs}}}{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{2}}({b}{{1,}}[GC{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[GC{i}U]{{{bs}}}{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}({b}{{1,}}[GC{i}U]{{{bs}}}{b}{{2,3}}|{b}{{2,3}}[GC{i}U]{{{bs}}}{b}{{1,}})))'.format(b = self.base, i = self.inv_base,  bs = self.bulge_size, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)

        elif self.mismatch_num == 1 and self.bulge_size == 0 and (self.base == 'G' or self.base == 'C'):
            regx = '((({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{3}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{2}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})))'.format(b = self.base, i = self.inv_base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)

        elif self.mismatch_num == 1 and self.bulge_size == 0 and (self.base == 'A' or self.base == 'T'):
            regx = '((({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{3}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{2}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})))'.format(b = self.base, i = self.inv_base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)

        elif self.mismatch_num == 2 and self.bulge_size == 0 and (self.base == 'G' or self.base == 'C'):
            regx = '((({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}}))|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[AT{i}U]{b}{{1,}})))'.format(b = self.base, i = self.inv_base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
    
        elif self.mismatch_num == 2 and self.bulge_size == 0 and (self.base == 'A' or self.base == 'T'):    
            regx = '((({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}})|(({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|(({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})([ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}){{2}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}}))|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}})|({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})[ATCGUN]{{{l},{lx}}}{b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}({b}{{2,}}|{b}{{1,}}[GC{i}U]{b}{{1,}})))'.format(b = self.base, i = self.inv_base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        else:
            regx = '(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}{b}{{{m},{x}}})'.format(b = self.base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        return regx    
    
    def quadruplex_finder(self, strand_dict, pattern):
        quadruplex_dict = {}
        for sequence_id, sequence in strand_dict.items():
            data_dict = self.quadruplex_finder_inner_function(sequence_id, sequence, pattern)
            for key, value in data_dict.items():
                data_dict[key] = [value]
            quadruplex_dict[ sequence_id ] = data_dict
        return quadruplex_dict
    
    def quadruplex_finder_inner_function(self, sequence_id, sequence ,pattern, QuadfromCpG = False):
        data_dict = {}
        if self.search_model == 'overlapping':
            start = 0
            original_seq = sequence
            seq_length = len(sequence)
            while start < seq_length:
                quadruplex = re.search(pattern, sequence)
                if not quadruplex:
                    break
                pos = original_seq.find(quadruplex.group(), start)
                sequence = original_seq[pos + 1 :]
                CpGisland_start = int(sequence_id) if QuadfromCpG else 0
                start_position = CpGisland_start + pos
                data_dict[ start_position ] = f"{quadruplex.group()}"
                start = pos + 1
        else:
            reg_obj = re.finditer(pattern, sequence)
            for match in reg_obj:
                CpGisland_start = int(sequence_id) if QuadfromCpG else 0
                start_position = CpGisland_start + match.start()
                data_dict[ start_position ] = f"{match.group()}"
        return data_dict

    def CpG_island(self, strand_dict, pattern):
        seq_dict_data = {} 
        quad_data_dict = {}
        for sequence_id, sequence in strand_dict.items():
            start = 0
            CpG_seq_dict = {}
            position_sequence_dict = {}
            while start < len(sequence): 
                x1 = start
                x2 = start + 200
                window = sequence[x1:x2]
                if len(window) < 200 :
                    break
                gc_content = 100 * (window.count("C") + window.count("G"))/len(window)
                if self.cpg_model == 'overlapping':
                    no_of_CpG = window.count("CG") or window.count("GC")
                else:
                    no_of_CpG = window.count("CG")
                CpG_ratio = no_of_CpG * len(window)/(window.count("C") * window.count("G"))
                if gc_content > 50 and CpG_ratio > 0.6:
                    position_sequence_dict[ x1+1 ] = window
                    CpG_seq_dict[ x1+1 ] = [x2, gc_content, CpG_ratio] 
                    if self.cpg_model != 'overlapping':
                        start += 200
                        continue
                start += 1
            
            if CpG_seq_dict:
                seq_dict_data[sequence_id] = CpG_seq_dict
            if position_sequence_dict:
                QuadfromCpG = True
                quad_dict = {}
                for sequence_id_, sequence in position_sequence_dict.items():
                    response_dict = self.quadruplex_finder_inner_function(sequence_id_, sequence, pattern, QuadfromCpG)
                    for key, value in response_dict.items():
                        key_set = quad_dict.get(key, set())
                        key_set.add(value)
                        quad_dict[key] = key_set
                if quad_dict:
                    #convert inner set to list to make it json writable
                    for key,value in quad_dict.items():
                        quad_dict[ key ] = list(value)
                    quad_data_dict[sequence_id] = quad_dict
        
        return seq_dict_data, quad_data_dict

def motif_file_header():
    return { 
        "header": ["Sequence ID", "G4 sequence", "Length", "Start Position", "End Position"]
    }

def cpg_file_header():
    return {
        "header": ["Sequence ID", "Position", "GC Content", "CpG Ratio"]
    }

'''
def generatePositionFigure(CpG_dict, quad_dict, outputfolderPath, strand):
    outputfolderPath = os.path.join(outputfolderPath, strand)
    if os.path.exists(outputfolderPath):
        shutil.rmtree(outputfolderPath)
    os.mkdir(outputfolderPath)
    motif_main_dict = {}
    cpg_main_dict= {}
    if quad_dict:
        for keys, values in quad_dict.items():                          
            for key, value in values.items():
                for item in value:
                    data_list = motif_main_dict.get(keys, [])
                    data_list.append(f'{key}\t{int(key)+len(item)}')
                    motif_main_dict[keys] = data_list
    if CpG_dict:
        for keys, values in CpG_dict.items(): 
            for key2, value2 in values.items():
                data_list = cpg_main_dict.get(keys, [])
                data_list.append(f'{key2}\t{value2[0]}')
                cpg_main_dict[keys]= data_list
    
    if motif_main_dict and cpg_main_dict:
        for coordinates in cpg_main_dict.values():
            cgi = []    
            for item in coordinates:        
                cgi.append([int(item.split('\t')[0]),int(item.split('\t')[1])])

        for coordinates2 in motif_main_dict.values():
            motif_lis = []
            for item in coordinates2:
                motif_lis.append([int(item.split('\t')[0]),int(item.split('\t')[1])])
        import matplotlib.pyplot as plt
        for ele in cgi:
            g4_start = [] 
            g4_end = []
            for ele2 in motif_lis:
                count = 0
                if ((ele2[0] >= ele[0]) and (ele2[1] <= ele[1])):  
                    g4_start.append(ele2[0])
                    g4_end.append(ele2[1])
                    y = range(0, len(g4_start))
                    plt.figure(figsize=(4,3))
                    plt.hlines(y, g4_start, g4_end, color = 'red', label = 'Quadruplexes')            
                    for x1,x2 in zip(g4_start, g4_end):
                        plt.text(x1 + 0.5, y[count], s = (x1,x2), horizontalalignment = 'right', fontweight = 'bold', verticalalignment='bottom',fontsize=4, color = 'blue')
                        count = count + 1
                    plt.xlim(ele)
                    plt.yticks([])
                    plt.xticks(ele)
                    plt.xlabel("Positions of Quadruplexes", fontsize = 6)
                    plt.legend(prop = {'style': 'italic'}, loc='upper left', bbox_to_anchor=(1, 0.5))
                    plt.savefig(fname = os.path.join(outputfolderPath,f'{ele[0]}_{ele[1]}_plot.png'), pad_inches = 0.1,quality = 95,bbox_inches = "tight")
                    plt.close()

'''