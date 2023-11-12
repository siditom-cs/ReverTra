import argparse
import os
from Bio import SeqIO
import pandas as pd
from utils import *


parser = argparse.ArgumentParser(description='Data partition for codon translation')
parser.add_argument('--rawdata_path', type=str, default='/home/tomer/ReverTra/data/raw_data',
                    help='location of the raw rna and protein abundance data.')
parser.add_argument('--data_path', type=str, default='/home/tomer/ReverTra/data/processed_data',
                    help='location of the processed data.')
parser.add_argument('--log_path', type=str, default='/home/tomer/ReverTra/data/logs',
                    help='location of the processed data.')
parser.add_argument('--agg_seqs_path', type=str, default='/home/tomer/ReverTra/data/logs/agg_seqs.temp.fasta',
                    help='location of the translated AA sequences of all sequences (for cd-hit).')
parser.add_argument('--cdhit_exec', type=str, default='/home/tomer/lib/cdhit/cd-hit',
                    help='location of the cd-hit program execution file.')
parser.add_argument('--cdhit_log', type=str, default='/home/tomer/ReverTra/data/logs/cdhit.log',
                    help='location of the cd-hit output log for the partition.')

args = parser.parse_args()

seq_dbs = [
        {'path':os.path.join(args.rawdata_path,"S_cerevisiae.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"S_cerevisiae.PA.dat"), 
            'species': "S_cerevisiae"},
        {'path':os.path.join(args.rawdata_path,"S_pombe.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"S_pombe.PA.dat"),
            'species': "S_pombe"}
        ]
seq_dbs = [
        {'path':os.path.join(args.rawdata_path,"E_coli.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"E_coli.csv"),
            'species': "E_coli"},
        {'path':os.path.join(args.rawdata_path,"B_subtilis.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"B_subtilis.csv"), 
            'species': "B_subtilis"},
        {'path':os.path.join(args.rawdata_path,"S_pombe.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"S_pombe.csv"),
            'species': "S_pombe"},
        {'path':os.path.join(args.rawdata_path,"S_cerevisiae.fasta"), 
            'pa_path':os.path.join(args.rawdata_path,"S_cerevisiae.csv"),
            'species': "S_cerevisiae"}


        ] 


##### Secion 1,2, &3: aggregate, and translate to aa seqs. 
seq_dict, trans_records = aggregate_species_records(seq_dbs)


SeqIO.write(trans_records, args.agg_seqs_path, "fasta")

print("Finished section 1, 2, & 3.")

##### Section 4: Run CD-HIT and load its results.
#Note: CD-HIT bak-file df is with columns: [ cls index | seq len | entryid | precent identity to the representative if exists]
os.system(args.cdhit_exec+" -d 10000 -i "+args.agg_seqs_path+" -c 0.7 "+" -o " + args.cdhit_log+" -bak 1")
#cls_df = pd.read_csv(args.cdhit_log+".bak.clstr",delimiter=r'\t| ',header=None,engine='python', names=['cls', 'len', 'entry', '-', 'precent_indentity'])
cls_df = pd.read_csv(args.cdhit_log+".bak.clstr",delimiter=r'\t| ',header=None,engine='python', names=[0,1,2,3,4])
for i in range(len(cls_df)):
    cls_df.iloc[i,2] = cls_df.iloc[i,2][7:-3]
print("Finished section 4 - finished cd-hit run, and uploaded cluster data.")

##### Section 5 - partition by clusters
nos, partition = split_data_by_cls(cls_df)
print("Finished section 5 - finished separting sequences by clusters.")

##### Section 6 - save records by partition and species
processed_data = process_data(seq_dict, # The aggregated index of all rna sequences
        cls_df, # The cluster index of each seq_id
        partition, # cluster indices partition  
        seq_dbs)
save_files(processed_data, args.data_path)

print("Finished section 6 - saved partition files by species.")

##### Section 7 - combine training and validation 
for seq_db in seq_dbs:
    combine_fasta_files(args.data_path, seq_db['species'])

print("Finished section 7 - combine training and validation .")
"""
##### Section 8 - add expression level for all sequences
for seq_db in seq_dbs:
    add_expr_level(args.data_path, seq_db['species'], seq_db['pa_path'])
print("Finished section 8 - add expression level.")
"""

