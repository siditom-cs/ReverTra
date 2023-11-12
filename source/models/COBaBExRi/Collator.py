import json
import os
import numpy as np
from transformers import BartConfig, BartModel, BartTokenizer
from transformers import AutoTokenizer
from datasets import load_dataset
import torch

class CodonTranslationCollator(object):
    """
    Description: Collator for COnTRA package. This collator enables codon masking to aa for the input sequence (default = 30%). 
                 Additionally, a noise mechanism of the decoder input ids is optional for introducing mistakes in the decoder for next token prediction.

    Note1: The subject can be a homologous codon sequence (from alignment), or the original input sequence.
    Note2: The subject_dna_seq key in features for single species differs in pretraining and fintuning. The former comrises from the original codon sequence (as in query_dna_seq), and the latter is already 100% masked.
           Thus, though 0% codon_mask_perc is the default for fintuning, it is 100% mask for single species, and 0% mask for homologous sequences.

    @Assumptions: 
        # input_ids:   [subject - species] [subject - codon seq] [query - species] [query - AA seq]
        # decoder_ids: <s> [query - species] [subject - codon seq][:-1]
        # labels:      [query - species] [subject - codon seq]

    """
    def __init__(self, tokenizer, masking_dict, restriction_dict, codon_mask_perc=0.3, noise_perc=0.0, rnd_seed=42, sw_aa_size=100, decoder_start_token_id=2):
        self.tokenizer = tokenizer
        self.masking_dict = masking_dict
        self.codon_mask_perc = codon_mask_perc
        np.random.seed(rnd_seed)
        self.sw_aa_size=sw_aa_size
        self.noise_perc = noise_perc
        self.restriction_dict = restriction_dict
        self.sw_aa_size = sw_aa_size
        self.decoder_start_token_id = decoder_start_token_id

    # input_ids:   [subject - species] [subject - codon seq] [query - species] [query - AA seq]
    # decoder_ids: <s> [query - species] [subject - codon seq][:-1]
    # labels:      [query - species] [subject - codon seq]
    def __call__(self, features):

        batch, mask_lens = self.get_batch(features)

        batch['decoder_input_ids'] = self.get_decoder_ids(batch['labels'])

        if self.codon_mask_perc > 0.0:
            mask_aa = self.prepare_mask(batch['input_ids'].shape, self.codon_mask_perc, batch['attention_mask'], mask_lens)
            batch['labels'] = torch.where(mask_aa[:,1:batch['labels'].shape[-1]+1],batch['labels'],-100)#the +1 is because the expr label is not included in 'labels'. 
            batch['input_ids'] = self.mask_input(batch['input_ids'], torch.where(mask_aa))
        if self.noise_perc > 0.0:
            mask_noise = self.prepare_mask(batch['decoder_input_ids'].shape, self.codon_mask_perc)
            batch['decoder_input_ids'] = self.mask_decoder_ids(batch['decoder_input_ids'], torch.where(mask_noise))
        print('InputIds: ',self.tokenizer.decode(batch['input_ids'][0,:]))
        print('Labels: ',batch['labels'][0,:])
        return batch

    def prepare_mask(self,shape, perc, attention_mask = None, per_vector_len_limits = None):
        mask = torch.tensor(np.random.binomial(1, perc, shape), dtype=torch.int64)
        if per_vector_len_limits != None:
            mask *= per_vector_len_limits
        if attention_mask != None:
            mask = mask*attention_mask
        mask[:,0] = 0
        mask = mask==1
        return mask

    def mask_input(self, inpt, listed_indices):
        for i,j in zip(listed_indices[0].numpy(), listed_indices[1].numpy()):  
            inpt[i,j] = self.masking_dict[str(inpt[i,j].item())]
        return inpt

    def mask_decoder_ids(self, decoder_input_ids, listed_indices):

        for i,j in zip(listed_indices[0].numpy(), listed_indices[1].numpy()):  
            decoder_input_ids[i,j] = self.masking_dict[str(decoder_input_ids[i,j].item())]
            k = np.random.choice(len(self.restriction_dict[str(decoder_input_ids[i,j].item())]))
            decoder_input_ids[i,j] = self.restriction_dict[str(decoder_input_ids[i,j].item())][k]
        return decoder_input_ids

    def get_expr_token(self, expr):
        if expr==None:
            expr = 'expr_unk'
        return "<"+expr+"> "

    def get_batch(self, features):
        lens = np.array([len(f['query_dna_seq']) for f in features])
        sinx = np.random.randint(0, lens)
        sinx = np.where(lens>self.sw_aa_size,sinx,0)
        einx = [min(si+self.sw_aa_size,l) for (si,l) in zip(sinx, lens)]
        query_win_codons = ["<"+f['query_species']+"> "+' '.join(f['query_dna_seq'][si:ei]) for (f,si,ei) in zip(features,sinx,einx)]
        subject_win_codons = ["<"+f['subject_species']+"> "+' '.join(f['subject_dna_seq'][si:ei]) for (f,si,ei) in zip(features,sinx,einx)]
        query_win_input_codons = ["<"+f['query_species']+"> "+self.get_expr_token(f['expr'])+ ' '.join(f['query_dna_seq'][si:ei]) for (f,si,ei) in zip(features,sinx,einx)]
        #query_win_aas = ["<"+f['query_species']+"> "+self.get_expr_token(f['expr'])+
        #        ' '.join(['<mask_'+aa+'>' if aa!='-' else '<gap>' for aa in f['qseq'][si:ei]]) for (f,si,ei) in zip(features,sinx,einx)]
        batch = self.tokenizer(query_win_input_codons, subject_win_codons, padding='longest', return_tensors='pt')
        batch.pop('token_type_ids')
        wins_lens = torch.tensor([len(w.split(" ")) for w in query_win_codons])
        mask_lens = torch.arange(0,batch['input_ids'].shape[-1]).expand(batch['input_ids'].shape) < wins_lens.unsqueeze(1)
        batch['labels'] = self.tokenizer(query_win_codons, padding='longest', return_tensors='pt')['input_ids'][:,1:-1]         
        return batch, mask_lens

    def get_decoder_ids(self, labels):  
        decoder_input_ids = labels.new_zeros(labels.shape)
        decoder_input_ids[:, 1:] = labels[:, :-1].clone()
        decoder_input_ids[:,0] = self.decoder_start_token_id
        return decoder_input_ids


