from Bio.Seq import Seq, translate
from Bio import SeqIO
import json
import os
import torch
from transformers import AutoTokenizer
import pickle

def init_codons():
    DNA_alphabet = ['T', 'C', 'A', 'G']
    codons = []
    for c1 in DNA_alphabet:
        for c2 in DNA_alphabet:
            for c3 in DNA_alphabet:
                codons.append(c1+c2+c3)
    return codons

def load_tokenizer(tokenizer_path):
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)
    masking_dict = {}
    with open(os.path.join(tokenizer_path,"masking_dict.json"),"r") as f:
        masking_dict = json.load(f)
    restriction_dict = {}
    with open(os.path.join(tokenizer_path,"mask_restrict_dict.json"),"r") as f:
        restriction_dict = json.load(f)
    #print('Tokenizer: ', tokenizer.vocab)
    #print("Masking_dict: ",masking_dict)
    #print("Restriction_dict: ",restriction_dict)
    return tokenizer, masking_dict, restriction_dict

class CAI(object):
    def __init__(self, training_data_path):
        super(CAI,self).__init__()
        self.training_data_path = training_data_path
        self.codons = init_codons()
        self.codons2idx = {c:i for i,c in enumerate(self.codons)}
        self.codons_translated_to_aa = self.split(Seq(''.join(self.codons)).translate())
        #self.unique_aa  =  ''.join(set(self.codons_translated_to_aa))
        self.unique_aa = [aa for aa in set(self.codons_translated_to_aa)]
        self.unique_aa.sort()
        self.codons_for_aa = self.init_AA_codons()
        self.cai_index_path = self.training_data_path + '.cai_idx.json'
        if os.path.exists(self.cai_index_path):
            with open(self.cai_index_path,'r') as f:
                d = json.load(f)
                self.CAI_index = d['CAI_index']
                self.CAI_probs = d['CAI_probs']
        else:
            CAI_index, CAI_probs, CAI_count = self.init_CAI_dict()
            d = dict()
            d['CAI_index'] = CAI_index
            d['CAI_probs'] = CAI_probs
            d['CAI_count'] = CAI_count
            self.CAI_index = CAI_index
            self.CAI_probs = CAI_probs
            with open(self.cai_index_path,'w') as f:
                json.dump(d,f)
        print(self.CAI_index)

    def get_cai_vectorized_representations(self, tokenizer):
        cai_vec_rep = {aa: torch.ones(1,len(tokenizer))*float("inf")*(-1) for aa in self.unique_aa}
        for aa in self.unique_aa:
            for ci in self.codons_for_aa[aa]:
                cai_vec_rep[aa][0,tokenizer.convert_tokens_to_ids(self.codons[ci])] = self.CAI_probs[self.codons[ci]]
        return cai_vec_rep
    
    def split(self, word):
        return [char for char in word]
    
    def init_AA_codons(self):
        #preparing a dictionary with the relevant codons for every aa:
        codons_for_aa = dict.fromkeys(self.unique_aa)
        for aa in self.unique_aa:
            codons_for_aa[aa] = [x for x,i in enumerate(self.codons_translated_to_aa) if i == aa]
        return codons_for_aa

    def normalize_codons_count(self, codons_count):
        probs = dict()
        for a in self.unique_aa:
            sum_AA = sum([codons_count[self.codons[ci]] for ci in self.codons_for_aa[a]])
            for ci in self.codons_for_aa[a]:
                probs[self.codons[ci]] = codons_count[self.codons[ci]]/sum_AA
        return probs

    def init_CAI_dict(self):
        dna_O = self.training_data_path
        recordsO = list(SeqIO.parse(dna_O, "fasta"))

        codons_countO = dict.fromkeys(self.codons, 0)
        for record in recordsO:
            words = [record.seq[i:i+3] for i in range(0, len(record.seq), 3)]    
            for codon in self.codons:
                codons_countO[codon] += words.count(codon)
        print(codons_countO)
        cai_codonO = dict.fromkeys(self.unique_aa)
        codons_probsO = self.normalize_codons_count(codons_countO)
        for a in self.unique_aa:
            count = 0
            for ci in self.codons_for_aa[a]:
                if codons_countO[self.codons[ci]] > count:
                    cai_codonO[a] = self.codons[ci]
                    count = codons_countO[self.codons[ci]]
        return cai_codonO, codons_probsO, codons_countO
    
    def getCAIfromDNASeq(self, dna_seq):
        s1 = dna_seq
        cai_codon_s1 = [self.CAI_index[v] for v in list(Seq(s1).translate())]
        return cai_codon_s1

    def calc_cai_acc(self, dna_seq):
        s1 = dna_seq
        true_codon_s1 = [s1[j:j+3] for j in range(0, len(s1), 3)]
        cai_codon_s1 = self.getCAIfromDNASeq(dna_seq)
        acc = sum(x == y for x, y in zip(true_codon_s1, cai_codon_s1))/len(true_codon_s1)

        cai_bases = "".join(cai_codon_s1)
        acc_bases = sum(x == y for x, y in zip(dna_seq, cai_bases))/len(dna_seq)
        return acc, acc_bases

    def calc_seq_cai(self, tokenizer, dna_seq):
        
        cai_vec_rep = self.get_cai_vectorized_representations(tokenizer)
        aa_seq = Seq(dna_seq).translate()
        cai_logits = torch.ones((len(aa_seq),len(tokenizer)))*float("inf")*(-1)
        for aai in range(len(aa_seq)):
            cai_logits[aai,:] = cai_vec_rep[aa_seq[aai]]
        cai_logits = cai_logits.log_softmax(dim=-1)
        cai_seq = tokenizer(" ".join(self.getCAIfromDNASeq(dna_seq))).input_ids[1:-1]
        return cai_seq, cai_logits

