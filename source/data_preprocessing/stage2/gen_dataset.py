import argparse
import os
import math
import numpy as np

from datasets import load_dataset
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

parser = argparse.ArgumentParser(description='Generate Arrow DataBase for Training.')
parser.add_argument('--data_path', type=str, default='/home/tomer/COnTRA/data/datasets/processed_data_SCPECBS3/homologs',
                    help='location of the partition fasta files.')
parser.add_argument('--dataset_name', type=str, default='SCPECBS3',
                    help='the name of the first species (blast query).')
parser.add_argument('--win_size', type=int, default='75',
                    help='.')


args=parser.parse_args()



def qdna2tokens(example):
    seq = example['query_dna_seq']
    return {'query_dna_seq': seq.split(" ") }

def sdna2tokens(example):
    seq = example['subject_dna_seq']
    return {'subject_dna_seq': seq.split(" ") }

    

def gendb(data_path):
    homologs_train_path = os.path.join(args.data_path,args.dataset_name+'_'+args.dataset_name+'_'+str(0)+'_homologs.csv')
    homologs_test_path = os.path.join(args.data_path,args.dataset_name+'_'+args.dataset_name+'_'+str(2)+'_homologs.csv')
    #homologs_val_path = os.path.join(args.data_path,args.species1+'_'+args.species2+'_'+str(1)+'_homologs.csv')
    homologs_val_path = os.path.join(args.data_path,args.dataset_name+'_'+args.dataset_name+
            '_'+str(1)+'_homologs.nonoverlaping_win'+str(args.win_size)+'.csv') 

    dataset_output_path = os.path.join(args.data_path,args.dataset_name+'_'+args.dataset_name+"_"+str(args.win_size)) 
    data_columns = ['qseqid', 'sseqid', 'query_species', 'subject_species', 'query_dna_seq', 'qseq', "subject_dna_seq", "sseq" ] 
     
    """
    df = pd.read_csv(homologs_test_path)
    df = df[[data_columns]]
    df.to_csv(, index=False)
    """

    data_files = {"train": homologs_train_path, "val": homologs_val_path, 'test': homologs_test_path}
    #data_files = {'test': homologs_test_path}
    
    codon_dataset = load_dataset("csv", data_files=data_files, delimiter=",")
    #codon_dataset = codon_dataset.remove_columns(['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    codon_dataset = codon_dataset.map(qdna2tokens)
    codon_dataset = codon_dataset.map(sdna2tokens)
    print(codon_dataset)

    codon_dataset.save_to_disk(dataset_output_path)
    #codon_dataset.save_to_disk(dataset_output_path+".eval")

        
if __name__ == "__main__":
    gendb(args.data_path)

