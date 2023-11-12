import json
import os
import sys
import pandas as pd
import warnings
from transformers import BartForConditionalGeneration
from datasets import load_from_disk
from .models import COBaBExRi, CAI, Harmonization
import argparse

parser = argparse.ArgumentParser(description='BART for Codon Optimization')
parser.add_argument('--eval_type', type=str, default='cai',
        help='Type of model evaluation: cai | hrm (harmonization) | model')
parser.add_argument('--config_path', type=str, default='./source/evaluation/config_example.json',
                    help='location of the app')
parser.add_argument('--out_filename', type=str, default='res_fixed.csv',
        help='')

args = parser.parse_args()

#import sys
def translation_AllPair_mask_codons(example):
    if example['qseqid']==example['sseqid']:
        return {'subject_dna_seq': ['<mask_'+aa+'>' for aa in example['sseq']]}
    else:
        return {'subject_dna_seq': example['subject_dna_seq']}


def single_species_mask(example):
    return {'subject_dna_seq':['<mask_'+aa+'>' if aa!='-' else '<gap>' for aa in example['qseq']]}


print(args)
def load_dataset(config):
    codon_dataset = load_from_disk(config['dataset_path'])
    del codon_dataset['train']
    del codon_dataset['val']
    print(codon_dataset)
    if config['mask_all']: #if single species leave only AA seq of query as input.
        print('mask all alignments')
        codon_dataset['test'] = codon_dataset['test'].map(single_species_mask)
    else: # mask only full sequences.
        print('only full seq')
        codon_dataset['test'] = codon_dataset['test'].map(translation_AllPair_mask_codons)
    print(codon_dataset)
    print(codon_dataset['test'][0])
    if config['eval_type'] == 'cai':
        codon_dataset = codon_dataset.filter(lambda example: example['query_species']==config['cai_query_species'])
        print(codon_dataset)
        #codon_dataset = codon_dataset.filter(lambda example: example['qseqid']==example['sseqid'])
    return codon_dataset

def load_model(config):
    if config['eval_type']=='cai':
        model = CAI.CAI(config['cai_refference_path']) 
    elif config['eval_type']=='hrm':
        model = Harmonization.NoMemoryHarmonization(config['hrm_refference_path']) 
    else:
        model = None
        if not config['checkpoint_flag']:
            print("Must enter a checkpoint model.")
            exit()
        model = BartForConditionalGeneration.from_pretrained(config['checkpoint_path'])

    return model

def get_example_data(config, example):
    return {k:example[k] for k in config['orig_dict']}

def dump_res(config, evalset, out_filename):
    df = pd.DataFrame(evalset)
    if not os.path.exists(config['outdir']):
        os.mkdir(config['outdir'])
    df.to_csv(os.path.join(config['outdir'], config['outfile']),index=False)
    if config['eval_type']!='cai':
        with open(os.path.join(config['outdir'], "eval_config.json"), "w") as handle:
            json.dump(config, handle, indent=4)

def check_dict(config, res):
    missing_values_flag = False
    key = None
    keys2check = [*config['out_dict'],*config['orig_dict']]
    for k in keys2check:
        if not k in res:
            missing_values_flag = True
            key = k
            break
    if missing_values_flag:
        warnings.warn("WARNING: Missing Value "+k+" In Result Dict")

if __name__ == '__main__':
    with open(args.config_path,"r") as f:
        print(args.config_path)
        config = json.load(f)
    print(config)
    module = [m for m in sys.modules.keys() if m.endswith(config['model_type'])]
    module = sys.modules[module[0]]
    
    codon_dataset = load_dataset(config)
    num_entries_in_ds = codon_dataset['test'].num_rows
    model = load_model(config)
    tokenizer, masking_dict, restriction_dict = module.load_tokenizer(config['tokenizer_path'])
    
    evalset = []
    for i in range(num_entries_in_ds):
        res = {'inx':i, 'eval_type':config['eval_type']}
        res.update(get_example_data(config, codon_dataset['test'][i]))
        if config['eval_type'] == 'cai':
            res.update(CAI.predict_cai(codon_dataset['test'][i], restriction_dict, tokenizer, model))
        elif config['eval_type'] == 'hrm':
            res.update(Harmonization.predict_harmonize(codon_dataset['test'][i], None, None, model))
        else: # args.eval_type == 'model':
            res.update(module.predict_model(config, codon_dataset['test'][i], restriction_dict, tokenizer, model))

        if config['debug']:
            check_dict(config, res)
        evalset.append(res)
        if i%500==0: #save every 1000 steps
            dump_res(config, evalset, args.out_filename)
    dump_res(config, evalset, args.out_filename)

