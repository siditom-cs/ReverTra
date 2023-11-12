import numpy as np
import torch
from Bio.Seq import Seq


def predict(example, mask_restriction_dict, tokenizer, model):
        query_seq = example['query_dna_seq']
        subject_seq = example['subject_dna_seq']
        query_aa_seq = example['qseq'] 
        spaced_seq = ' '.join([c for c in example["query_dna_seq"] if c != '<gap>'])
        true_vals = [c for c in query_seq if c!='<gap>']

        accuracy, harmonized_seq = model.calc_harmonic_acc(subject_seq, query_seq, example['subject_species'], example['query_species']) 

        res = dict()
        res['num_of_correct_predicted_codons'] = sum([1 for x,y in zip(true_vals, harmonized_seq.split(" ")) if (x==y) ])
        res['prot_len'] = len(true_vals)
        res['query_codons'] = spaced_seq
        res['subject_codons'] = " ".join(example['subject_dna_seq'])
        res['pred_codons'] = harmonized_seq
        res['cross_entropy_loss'] = np.nan
        res['entropy'] = np.nan
        res['accuracy'] = res['num_of_correct_predicted_codons'] / res['prot_len']
        
        return res
