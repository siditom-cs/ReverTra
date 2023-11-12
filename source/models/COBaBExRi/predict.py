import json
import os
import sys
import numpy as np
import pandas as pd

from transformers import LogitsProcessor,LogitsProcessorList, BartForConditionalGeneration, AutoTokenizer
from datasets import load_from_disk, load_metric
import torch
from Bio.Seq import Seq

from .utils import translation_AllPair_mask_codons 
from ..CAI.CAI import CAI

from .restrict import RestrictToAaLogitsWarper
import argparse
parser = argparse.ArgumentParser(description='BART for Codon Optimization')
parser.add_argument('--config_path', type=str, default='./source/evaluation/config_example.json',
                    help='location of the app')
args = parser.parse_args()


def get_expr_token(expr):
    if expr==None:
        expr = 'expr_unk'
    return "<"+expr+"> "
def prepare_inputs(example, tokenizer, config):
        query_seq = example['query_dna_seq']
        subject_seq = example['subject_dna_seq']
        query_aa_seq = example['qseq']
        query_species = example['query_species']
        subject_species = example['subject_species']

        #Rename
        qaaseq, query_species, sw_aa_size = example['qseq'], example['query_species'],  config['sw_aa_size']
        if config['inference_type'] == 'mimic':
            sseq, subject_species = example['subject_dna_seq'].split(" "), example['subject_species']
        else:
            sseq, subject_species = ['<mask_'+aa+'>' if aa!='-' else '<gap>' for aa in qaaseq], example['query_species']

        #Prepare fixed-sized windows
        query_aa_wins = [qaaseq[i:i+sw_aa_size] for i in range(0,max(1,len(qaaseq)-sw_aa_size+1))]
        subject_dna_wins = [" ".join(sseq[i:i+sw_aa_size]) for i in range(0,max(1,len(qaaseq)-sw_aa_size+1))]
        mask_aa_wins = ["<"+query_species+"> "+" ".join(['<mask_'+aa+'>' if aa!='-' else '<gap>' for aa in wseq]) for wseq in query_aa_wins]
        query_aa_wins = ["<"+query_species+"> "+get_expr_token(example['expr'])+' '.join(['<mask_'+aa+'>' if aa!='-' else '<gap>' for aa in wseq]) for wseq in query_aa_wins]
        subject_dna_wins = ["<"+subject_species+"> "+wseq for wseq in subject_dna_wins]

        #Encode windows
        input_ids = tokenizer(query_aa_wins, subject_dna_wins, return_tensors="pt", padding='max_length', max_length=sw_aa_size*2+3).input_ids
        masked_ids = tokenizer(mask_aa_wins, return_tensors="pt").input_ids[:,1:-1]
        return input_ids,masked_ids


def generate_outputs(input_ids, masked_ids, mask_restriction_dict, model, sw_aa_size):
        masked_restriction = masked_ids
        
        logits_processor = LogitsProcessorList(
                [RestrictToAaLogitsWarper(masked_restriction, mask_restriction_dict)])
        maxlen = min((sw_aa_size+3),masked_restriction.shape[-1]+2)
           
        if masked_ids.shape[0]<1000:
            input_ids = input_ids.cuda()
            masked_ids = masked_ids.cuda()#TODO
            model = model.cuda()#TODO
        else:
            masked_ids = masked_ids.cpu()#TODO
            model = model.cpu()#TODO

        outputs = model.generate(input_ids, do_sample=False, output_scores = True, return_dict_in_generate = True, renormalize_logits = True, logits_processor=logits_processor, max_length=maxlen)#max_new_token=maxlen)
        outputs = torch.stack(outputs['scores'][:sw_aa_size+1],1)
        return outputs

def calc_combined_gen_from_sliding_windows_logits(sw_logits, seqlen, sw_aa_size):
        sw_logits = sw_logits.cuda()
        collect_logits = torch.zeros([seqlen, sw_logits.shape[-1]]).cuda()
        counts = torch.zeros([1,seqlen]).cuda()
        most_freq_pred = torch.zeros([seqlen,1])
        # This segment aggregates (sums) the logits of the different windows. Only the relevant codons (restricted by AA) are sumed.
        for i in range(sw_logits.shape[0]): # window num 
            for j in range(min(sw_aa_size, seqlen)): # sequence len - codon index
                collect_logits[i+j, :] += torch.exp(sw_logits[i, 1+j, :])
                counts[0,i+j] += 1
        
        
        for i in range(seqlen):
            collect_logits[i,:] /= counts[0,i]
    
        collect_logits = torch.log(collect_logits)
        collect_logits = collect_logits.log_softmax(dim=-1)

        for i in range(seqlen):
            most_freq_pred[i] = torch.argmax(collect_logits[i,:]).item()
        #add later, the probabilities of each prediction (from freq)
        return collect_logits, most_freq_pred


def predict(config, example, mask_restriction_dict, tokenizer, model):

    input_ids, masked_ids = prepare_inputs(example, tokenizer, config)
    outputs = generate_outputs(input_ids, masked_ids, mask_restriction_dict, model, config['sw_aa_size'])
    logits, most_freq_pred = calc_combined_gen_from_sliding_windows_logits(outputs, len(example['qseq']), config['sw_aa_size'])

    ce = torch.nn.CrossEntropyLoss()
    most_freq_pred=most_freq_pred.clone().detach().reshape((1,-1))


    #print("decode: ", tokenizer.decode(most_freq_pred.numpy().astype(int)[0]))
    #print("truevals: ", tokenizer.decode(true_vals))
    res = dict()

    res['prot_len'] = len(example['qseq'])
    res['prot_AAs'] = example['qseq']
    res['pred_codons'] = tokenizer.decode(most_freq_pred.numpy().astype(int)[0])
    res['entropy'] = (-torch.nan_to_num(torch.exp(logits)*logits,nan=0.0).sum(dim=-1)).mean().item()

    assert(res['prot_len']==len(res['pred_codons'].split(" ")))

    if config['calc_stats'] and 'query_dna_seq' in example.keys():
        true_vals = tokenizer(example['query_dna_seq'], return_tensors="pt").input_ids[:,1:-1]
        mask = true_vals > config['special_token_th'] #special tokens threshold
        true_vals = true_vals.tolist()[0]
        masked_most_freq_pred = most_freq_pred.masked_select(mask).numpy().astype(int)
        masked_true_vals = torch.tensor(true_vals).masked_select(mask).numpy().astype(int)

        res['subject_codons'] = example['subject_dna_seq']
        res['num_of_correct_predicted_codons'] = sum([int(x==y) for x,y in zip(masked_true_vals, masked_most_freq_pred)])
        res['query_codons'] = example['query_dna_seq']      
        res['cross_entropy_loss'] = ce(logits, torch.tensor(true_vals)).item()
        res['perplexity'] = np.exp(res['cross_entropy_loss'])
        res['accuracy'] = res['num_of_correct_predicted_codons'] / res['prot_len']
    #print(example['qseqid'], example['sseqid'],res['cross_entropy_loss'], res['entropy'],res['accuracy'])
    return res

