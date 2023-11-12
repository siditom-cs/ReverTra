import numpy as np
import torch
import torch.nn.functional as F
from ..CAI.CAI import CAI

def dna2codons(example):
    seq = example['dna_seq']
    return {'codon_seq': ' '.join([seq[i:i+3] for i in range(0,len(seq),3)]) }

def split_codons(example):
    return {'codon_seq': example['codon_seq'].upper().split(" ") }

def translation_query2subject_codons(example):
    if example['query_species']==example['subject_species']:
        return {'subject_dna_seq': example['query_dna_seq']}
    else:
        return {'subject_dna_seq': example['subject_dna_seq']}

def translation_AllPair_mask_codons(example):
    if example['qseqid']==example['sseqid']:
        return {'subject_dna_seq': ['<mask_'+aa+'>' for aa in example['sseq']]}
    else:
        return {'subject_dna_seq': example['subject_dna_seq']}

def dna2bases(example):
    seq = example['dna_seq']
    return {'bases': ' '.join([seq[i] for i in range(0, len(seq))]) }

def aa2masked_aa(example):
    masked_aa_seq = []
    for aa in example['aa_seq']:
        if aa == '_':
            masked_aa_seq.append('<gap>')
        else:
            masked_aa_seq.append('<mask_'+aa+'>') 
    masked_aa_seq = ' '.join(masked_aa_seq)
    return {'masked_aa_seq': masked_aa_seq}

def get_tokenize(tokenizer):
    def tokenize(example):
        seq = example['codon_seq']
        return tokenizer(seq)
    return tokenize

def get_tokenize_codons(tokenizer):
    def tokenize(example):
        seqcodons = example['codon_seq']
        return tokenizer(seqcodons)
    return tokenize

def get_tokenize_codons_wSpeciesFlag(tokenizer):
    def tokenize(example):
        seqcodons = example['codon_seq']
        species_flag = "<"+example['species']+">"
        return tokenizer(seqcodons,species_flag)
    return tokenize

def get_tokenize_bases(tokenizer):
    def tokenize(example):
        seqbases = example['bases']
        return tokenizer(seqbases)
    return tokenize

def cai_metric(tokenizer, path=None):
    cai = CAI(path)

    def calc_cai(masked_labels):
        total_dna = ''.join(tokenizer.decode(masked_labels)).replace(' ','')
        acc, acc_bases = cai.calc_cai_acc(total_dna)
        return acc, acc_bases

    return calc_cai, cai

def get_metrics(tokenizer, training_path=None, special_token_th=31):

    compute_cai, cai_model = cai_metric(tokenizer, training_path)
    def compute_perplexity(predictions, labels, mask):
        predictions = predictions.transpose(2,1)
        loss = F.cross_entropy(predictions, labels, ignore_index=-100)
        #loss = torch.masked_select(loss, mask)
        #loss = torch.sum(loss)
        perplexity = torch.exp(loss)
        return perplexity

    def compute_metrics_combined(eval_preds):
        logits, labels = eval_preds
        predictions = torch.tensor(np.argmax(logits[0], axis=-1))
        labels = torch.tensor(labels)
        masked_index = torch.ne(labels, -100)
        special_token_index = torch.gt(labels, special_token_th)

        mask = torch.logical_and(masked_index, special_token_index)

        mask_codons = torch.lt(labels, 154)
        mask_bases = torch.gt(labels, 153)

        mask_codons =  torch.logical_and(mask, mask_codons)
        mask_bases =  torch.logical_and(mask, mask_bases)


 
        acc = torch.tensor(predictions == labels, dtype=torch.float32)
        acc_all = torch.masked_select(acc,mask)
        acc_codons = torch.masked_select(acc,mask_codons)
        acc_bases = torch.masked_select(acc,mask_bases)

        cai_score, cai_bases_score = compute_cai(torch.masked_select(labels,mask_codons))
        pp = compute_perplexity(torch.tensor(logits[0]), labels, mask)
        #pred_cai_score = compute_cai(torch.masked_select(predictions,mask))
        #return {"acc": torch.mean(acc), "cai": cai_score, "pred-cai":pred_cai_score} 
        return {"acc": torch.mean(acc_all), "acc_bases": torch.mean(acc_bases),"acc_codons": torch.mean(acc_codons), "cai": cai_score, "cai_bases":cai_bases_score,"pp":pp, "acc_diff":(torch.mean(acc_codons)-cai_score)} 

    def compute_metrics(eval_preds):
        logits, labels = eval_preds
        predictions = torch.tensor(np.argmax(logits[0], axis=-1))
        labels = torch.tensor(labels)
        masked_index = torch.ne(labels, -100)
        special_token_index = torch.gt(labels, special_token_th)

        mask = torch.logical_and(masked_index, special_token_index)
        acc = torch.tensor(predictions == labels, dtype=torch.float32)
        acc = torch.masked_select(acc,mask)
        acc_mean = torch.mean(acc)

        cai_score, cai_bases_score = compute_cai(torch.masked_select(labels,mask))
        pp = compute_perplexity(torch.tensor(logits[0]), labels, mask)
        #pred_cai_score = compute_cai(torch.masked_select(predictions,mask))
        #return {"acc": torch.mean(acc), "cai": cai_score, "pred-cai":pred_cai_score} 
        return {"acc": acc_mean, "cai": cai_score, "cai_bases":cai_bases_score,"pp":pp, "acc_diff":(acc_mean-cai_score)} 
    return compute_metrics


