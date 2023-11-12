import json
import os
import torch
from torch import nn
import wandb
import transformers
from transformers import BartConfig, BartForConditionalGeneration
from transformers import AutoTokenizer, TrainingArguments, Trainer
from datasets import load_from_disk, load_metric

from .Collator import (
        CodonTranslationCollator,
        )

from .utils import (
        get_tokenize_codons,
        split_codons,
        aa2masked_aa,
        translation_AllPair_mask_codons,
        get_metrics,
    )

import argparse
parser = argparse.ArgumentParser(description='BART for Codon Optimization')
parser.add_argument('--config_path', type=str, default='./source/models/config.json',
                    help='location of the app')
args = parser.parse_args()


def load_dataset(dataset_path, single_species_flag=False, finetune_flag=False):
    codon_dataset = load_from_disk(dataset_path)
    print("Codon Dataset: ",codon_dataset)
    if single_species_flag:
        #Fileter all trainslation entries from the dataset. For each entry in the dataset, input and output for the model are only from the same species.
        codon_dataset = codon_dataset.filter(lambda example: example['sseqid']==example['qseqid'])
        #print(codon_dataset["train"])

    #Finetune flag - for single species entries, the codon sequence input is 100% masked.
    codon_dataset = codon_dataset.map(translation_AllPair_mask_codons)
        #print(codon_dataset)
    return codon_dataset

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

def load_model(config, checkpoint_flag=False, finetune_flag=False):
    model = None
    if (not checkpoint_flag) and finetune_flag:
        print("Must have a pretrained checkpoint for translation fintuning.")
        exit()     
    elif (not checkpoint_flag) and (not finetune_flag):
        #Pretraining - initializing random weights to network.
        model_config = config['model_config']
        configuration = BartConfig(vocab_size = model_config['vocab_size'], d_model = model_config["d_model"], 
                encoder_layers=model_config["attention_layers"], decoder_layers=model_config["attention_layers"], 
                encoder_attention_heads = model_config["attention_heads"], decoder_attention_heads = model_config["attention_heads"] , 
                encoder_ffn_dim = model_config["ffn_dim"], decoder_ffn_dim = model_config["ffn_dim"], 
                max_position_embeddings = model_config['max_position_embeddings'] )
        model = BartForConditionalGeneration(configuration)
    elif checkpoint_flag:
        #Finetuning or second pretrainings - loading model weights from a checkpoint.
        model = BartForConditionalGeneration.from_pretrained(config['checkpoint_path'])
        
    #print("Model - Num of parameters: ", model.num_parameters())
    return model

def init_training_args(train_config, model_path, dataset_size):
    batch_size = train_config['batch_size']
    logging_steps = train_config['logging_steps']
    label_smoothing = train_config['label_smoothing_factor']

    steps_in_epoch  = dataset_size/batch_size
    number_of_steps = train_config['steps']
    number_of_epochs = int(number_of_steps/steps_in_epoch)

    training_args = TrainingArguments(
            report_to="wandb",
            #report_to=None,
            output_dir=f"{model_path}",
            warmup_steps=train_config["warmup_steps"],
            overwrite_output_dir=True,
            evaluation_strategy="steps",
            #weight_decay=0.01,
            learning_rate=train_config["starting_lr"],
            per_device_train_batch_size=batch_size,
            per_device_eval_batch_size=batch_size,
            push_to_hub=False,
            fp16=True,
            logging_steps=logging_steps,
            num_train_epochs= number_of_epochs,
            label_smoothing_factor=label_smoothing,
            load_best_model_at_end=True,
            save_steps=logging_steps,
            save_total_limit=3,
            logging_first_step=True,
            disable_tqdm=False,
            remove_unused_columns=False,
            seed=train_config['seed'],
            #metric_for_best_model='eval_acc_diff'
        )
    return training_args

def init_trainer(training_args, tokenizer, dataset, masking_dict, restriction_dict, config):
    metrics = get_metrics(tokenizer, config['cai_refference_path'],special_token_th=config['special_token_th'])
    training_collator = CodonTranslationCollator(tokenizer,
            masking_dict,
            restriction_dict,
            codon_mask_perc=config['train_config']['mask_perc'], #Should be zero - 100% mask on single species (in the dataset), and 0% on translation (homologs). 
            noise_perc = config['train_config']['decoder_noise'],
            rnd_seed=config['train_config']['seed'], 
            sw_aa_size=config['sw_aa_size'])

    if not os.path.exists(training_args.output_dir):
        os.mkdir(training_args.output_dir)

    with open(os.path.join(training_args.output_dir,"configs.log"),"w") as handle:
        json.dump(config, handle)
    
    btrainer = Trainer(model=model, args=training_args,
                   train_dataset=dataset["train"],
                   eval_dataset=dataset["val"],
                   data_collator=training_collator,
                   tokenizer=tokenizer,
                   compute_metrics=metrics,
                   )
    return btrainer


""" Training Script:
1. Load Data
2. Augment Data 
   2.a Single Species Training
   2.b Finetuning Augmentation - replace subject with mask for single species
3. Load Tokenizer and additional dicts
4. Define/Load Model 
5. Define Training Configurations
6. Define Huggingface Trainer 
7. Train Model
"""
if __name__ == '__main__':
    #Load traingin configuration
    if not os.path.exists(args.config_path):
        print("No training configuration.")
        exit()
    with open(args.config_path,"r") as f:
        config = json.load(f)
    print(config)
    wandb.init(project=config['project_name'], name=config['model_name'])#tags=['overfit_test', 'small_data'])
    dataset = load_dataset(config['dataset_path'], single_species_flag=config['dataset_single_species_flag'], finetune_flag=config['finetune_flag'])
    tokenizer, masking_dict, restriction_dict = load_tokenizer(config['tokenizer_path'])
    model = load_model(config, checkpoint_flag=config['checkpoint_flag'], finetune_flag=config['finetune_flag'])
    train_args = init_training_args(config['train_config'], os.path.join(config['model_outpath'],config['model_name']), len(dataset['train']))
    trainer = init_trainer(train_args, tokenizer, dataset, masking_dict, restriction_dict, config)
    trainer.train()
    trainer.save_model(os.path.join(train_args.output_dir,'best_model'))
