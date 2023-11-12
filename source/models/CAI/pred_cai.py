import numpy as np
import torch
from Bio.Seq import Seq


def predict(example, mask_restriction_dict, tokenizer, model):
        query_seq = example['query_dna_seq']
        subject_seq = example['subject_dna_seq']
        query_aa_seq = example['qseq'] 
        spaced_seq = ' '.join([c for c in example["query_dna_seq"] if c != '<gap>'])
        true_vals = tokenizer(spaced_seq, return_tensors="pt").input_ids[:,1:-1]
        true_vals = true_vals.tolist()[0] 

        cai_seq, cai_logits = model.calc_seq_cai(tokenizer, spaced_seq.replace(' ',''))
        ce = torch.nn.CrossEntropyLoss()

        res = dict()
        res['num_of_correct_predicted_codons'] = sum([1 for x,y in zip(true_vals, cai_seq) if (x==y) ])
        res['prot_len'] = len(true_vals)
        res['query_codons'] = spaced_seq
        res['subject_codons'] = " ".join(example['subject_dna_seq'])
        res['pred_codons'] = tokenizer.decode(cai_seq)
        res['cross_entropy_loss'] = ce(cai_logits.cuda(), torch.tensor(true_vals).cuda()).item()
        res['entropy'] = (-torch.nan_to_num(torch.exp(cai_logits)*cai_logits,nan=0.0).sum(dim=-1)).mean().item()
        res['accuracy'] = res['num_of_correct_predicted_codons'] / res['prot_len']
        
        return res
