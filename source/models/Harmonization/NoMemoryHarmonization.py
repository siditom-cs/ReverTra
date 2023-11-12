import pandas as pd
import json
import os
from Bio.Seq import Seq, reverse_complement, transcribe, back_transcribe, translate

def init_codons():
    DNA_alphabet = ['T', 'C', 'A', 'G']
    codons = []
    for c1 in DNA_alphabet:
        for c2 in DNA_alphabet:
            for c3 in DNA_alphabet:
                codons.append(c1+c2+c3)
    return codons

class Harmonizer(object):
    def __init__(self, dataset_path='.data/datasets/processed_data_SCPECBS3', 
            species_list=['S_cerevisiae', 'S_pombe', 'E_coli', 'B_subtilis'], 
            cai_type='.0.nt.fasta.cai_idx.json'):

        super(Harmonizer,self).__init__()
        self.dataset_path = dataset_path
        self.codons = init_codons()
        self.cai_species_index = self.init_cai_translation_dict(dataset_path,species_list, cai_type)



    def init_CAI(self, path):
        with open(path,"r") as handle:
            cai_data = json.load(handle)

        AAs = set(Seq.translate("".join(self.codons)))
        cai_index = {}
        for aa in AAs:

            #Get codons for this AA
            codons_for_aa_codons = [(c,cai_data['CAI_probs'][c]) for c in self.codons if translate(c)==aa]
            codons_for_aa_codons.sort(key=lambda x: x[1], reverse=True)
            cai_index[aa] = {c[0]:i for i,c in enumerate(codons_for_aa_codons)}

        return cai_index

    def init_cai_translation_dict(self, dataset_path, species_list, cai_type):
        d = {}
        for s in species_list:
            d[s] = self.init_CAI(os.path.join(os.path.join(dataset_path,s),s+cai_type))
        return d


    def calc_harmonic_acc(self,origin_seq, dest_seq, origin_species='S_pombe', destination_species='S_cerevisiae'):
        orig_index = self.cai_species_index[origin_species]
        dest_index = self.cai_species_index[destination_species]

        dest_pred_seq = []
        dest_aaidentity = []

        for i in range(len(dest_seq)):
            if dest_seq[i] == '<gap>':
                # There is a gap in the alignment - there is no AA/Codon for S_cerevisiae seq in this location.
                continue
            dest_aa = translate(dest_seq[i])
            orig_codon = origin_seq[i]
            if 'mask' in orig_codon:
                # There is a gap in the alignment - there is no codon in s_pombe sequence - taking the most frequent codon for S_cerevisiae
                dest_pred_codon = list(dest_index[dest_aa].keys())[0]
                dest_pred_seq.append(dest_pred_codon)
                continue

            orig_aa = translate(orig_codon)
            orig_codon_loc = orig_index[orig_aa][orig_codon]

            if orig_codon_loc>len(dest_index[dest_aa].keys())-1:
                dest_pred_codon = list(dest_index[dest_aa].keys())[-1]
            else:
                dest_pred_codon = list(dest_index[dest_aa].keys())[orig_codon_loc]
            if orig_aa == dest_aa:
                dest_aaidentity.append(1)
            else:
                dest_aaidentity.append(0)
    
            dest_pred_seq.append(dest_pred_codon)

        # Calculate Accuracy:
        dest_seq = list(filter(lambda a: a != '<gap>', dest_seq))
        accuracy = sum([x==y for (x,y) in zip(dest_seq,dest_pred_seq)])/len(dest_seq)
        """
        print(sc_aaidentity)
        accuracy = sum([x == y for (x, y, z) in zip(s_cerevisiae_seq, sc_pred_seq, sc_aaidentity) if z == 1]) / sum(
            sc_aaidentity)
        """
        return accuracy, ' '.join(dest_pred_seq)

if __name__=="__main__":
    print(cai_species_index)
    print(cai_species_index.keys())





