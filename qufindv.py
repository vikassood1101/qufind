"""
Created on Sun Sep 12 20:20:43 2021
@author: priyasharma16
"""
import re
import json, os
import G4plotv as plt

class qufindv():
    def __init__(self, outputfolderPath, seq_1_dict, seq_2_dict, search_model, nucleobase, min_stem_size, max_stem_size, min_loop, max_loop, bulge_size=0, mismatch_num=0):
        self.outputfolderPath = outputfolderPath
        self.nucleobase = nucleobase
        self.search_model = search_model
        self.min_stem = min_stem_size
        self.max_stem = max_stem_size
        self.loopmin = min_loop
        self.loopmax = max_loop
        self.bulge_size = bulge_size
        self.mismatch_num = mismatch_num

        self.seq_1_dict = seq_1_dict
        self.seq_2_dict  = seq_2_dict
        
        self.base, self.inv_base = self.make_base()
        reg_ex = self.regex_func() 

        self.seq_motif_dict_1, self.seq_motif_dict_2 = {}, {}

        quad_dict_1 = self.quadruplex_finder(self.seq_1_dict, reg_ex)
        if quad_dict_1:
            quad_dict_1 = self.processDictionaryData(quad_dict_1)
            quad_dict = self.sortResult(quad_dict_1)
            self.seq_motif_dict_1 = quad_dict
            writable_dict = motif_file_header()
            writable_dict.update(quad_dict)
            file_name_write = os.path.join(self.outputfolderPath, 'output_seq_1_motif.json')
            with open(file_name_write, 'w') as fp_write:
                json.dump(writable_dict, fp_write)

        quad_dict_2 = self.quadruplex_finder(self.seq_2_dict, reg_ex)
        if quad_dict_2:
            quad_dict_2 = self.processDictionaryData(quad_dict_2)
            quad_dict = self.sortResult(quad_dict_2)
            self.seq_motif_dict_2 = quad_dict
            writable_dict = motif_file_header()
            writable_dict.update(quad_dict)
            file_name_write = os.path.join(self.outputfolderPath, 'output_seq_2_motif.json')
            with open(file_name_write, 'w') as fp_write:
                json.dump(writable_dict, fp_write)
        
        if self.seq_motif_dict_1 and self.seq_motif_dict_2:
            final_dict_1 = {}
            for key, value in self.seq_motif_dict_1.items():
                final_dict_1[key] = []
                for i_key, i_value in value.items():
                    for ele in i_value:
                        final_dict_1.get(key, []).append(f'{i_key}\t{len(ele)}\t{ele}')
            
            final_dict_2 = {}
            for key, value in self.seq_motif_dict_2.items():
                final_dict_2[key] = []
                for i_key, i_value in value.items():
                    for ele in i_value:
                        final_dict_2.get(key, []).append(f'{i_key}\t{-len(ele)}\t{ele}')
            
            seq_len_list = [len(value) for value in seq_1_dict.values()] 
            seq_len_list.append(len(value) for value in seq_2_dict.values())
            plt.drawG4Plot(self.outputfolderPath, seq_len_list, final_dict_1, final_dict_2)

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
        
        elif self.mismatch_num == 0 and self.bulge_size == 0 and (self.base == 'G' or self.base == 'C'):
            regx = '(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}{b}{{{m},{x}}})'.format(b = self.base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        
        elif self.mismatch_num == 0 and self.bulge_size == 0 and (self.base == 'A' or self.base == 'T'):
            regx = '(({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}{b}{{{m},{x}}})'.format(b = self.base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        else:
            regx = '([ATC]{{0,5}}({b}{{{m},{x}}}[ATCGUN]{{{l},{lx}}}){{3}}{b}{{{m},{x}}}[ATC]{{0,5}})'.format(b = self.base, m = self.min_stem, x = self.max_stem, l = self.loopmin, lx = self.loopmax)
        return regx    
    
    def quadruplex_finder(self, strand_dict, pattern):
        quadruplex_dict = {}
        for sequence_id, sequence in strand_dict.items():
            data_dict_result = {}
            data_dict = self.quadruplex_finder_inner_function(sequence, pattern)
            for key, value in data_dict.items():
                if re.search('((A){10,}|(T){10,})',value) is None:               
                    data_dict_result[key] = [value]
            
            if data_dict_result:
                quadruplex_dict[ sequence_id ] = data_dict_result
        return quadruplex_dict
    
    def quadruplex_finder_inner_function(self, sequence ,pattern):
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
                data_dict[ pos+1 ] = f"{quadruplex.group()}"
                start = pos + 1
        else:
            reg_obj = re.finditer(pattern, sequence)
            for match in reg_obj:
                data_dict[ match.start()+1 ] = f"{match.group()}"
        return data_dict

def motif_file_header():
    return { 
        "header": ["Sequence ID", "G4 sequence", "Length", "Start Position", "End Position"]
    }


    
    