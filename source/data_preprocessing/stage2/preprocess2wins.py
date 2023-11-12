import argparse
import os
import math
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

parser = argparse.ArgumentParser(description='Homologs windows processing')
parser.add_argument('--data_path', type=str, default='/home/tomer/CodOpTRM/data/homologs_data/processed_test',
                    help='location of the partition fasta files.')
parser.add_argument('--log_path', type=str, default='/home/tomer/CodOpTRM/data/homologs_data/logs',
                    help='location of the partition blast output files.')
parser.add_argument('--species1', type=str, default='S_cerevisiae',
                    help='the name of the first species (blast query).')
parser.add_argument('--species2', type=str, default='S_pombe',
                    help='the name of the second species (blast db).')
parser.add_argument('--win_size', type=int, default='75',
                    help='sw - window size of the query homolog sequence.')
parser.add_argument('--sw_jump', type=int, default='5',
                    help='sw - window size of the query homolog sequence.')
parser.add_argument('--partition_id', type=int, default='0',
                    help='sw - window size of the query homolog sequence.')
parser.add_argument('--with_gaps', type=int, default='1',
                    help='Whether to leave gaps in the aliged windows or not.')
parser.add_argument('--mode', type=str, default='homologs',
        help='Script mode: homologs or augmented.')
parser.add_argument('--cut_perc', type=float, default='0.5',
                    help='What precentage of the sequence will be masked in augmented mode. cut_perc=0.0 means sournce=target.')
parser.add_argument('--msk_mean', type=int, default='3',
                    help='What is the gap mean size of the augmented sequence in mode=augmented. msk_mean=0 means source=target.')


args=parser.parse_args()



def alignment2wins(query_seq_dna, query_seq_aa, subject_seq_dna, subject_seq_aa):
    aa_ws = args.win_size
    ws = min(aa_ws, len(query_seq_dna)) 
    #ws is index in codon string
    seq = query_seq_dna
    nwins = int(np.ceil(len(seq)/ws))
    win_dict = {"query_dna_seq":[], "qseq":[],"subject_dna_seq":[], "sseq":[]}
    for i in range(0,nwins):
        si = i*ws
        win_dict["query_dna_seq"].append(" ".join(query_seq_dna[i*ws:min((i+1)*ws,len(query_seq_dna))]))
        win_dict["qseq"].append(query_seq_aa[int(i*ws):min(int((i+1)*ws),len(query_seq_aa))])
        win_dict["subject_dna_seq"].append(" ".join(subject_seq_dna[i*ws:min((i+1)*ws,len(subject_seq_dna))]))
        win_dict["sseq"].append(subject_seq_aa[int(i*ws):min(int((i+1)*ws),len(subject_seq_aa))])
    return win_dict,nwins

def preprocess_homologs(data_path):
    homologs_pairs = pd.read_csv(os.path.join(args.data_path,args.species1+'_'+args.species2+'_'+str(args.partition_id)+'_homologs.csv'))
    
    #meta_data_columns = ['qseqid', 'sseqid', 'query_species', 'subject_species'] 
    gen_columns = ['query_dna_seq', 'qseq', "subject_dna_seq", "sseq"]
    meta_data_columns = set(homologs_pairs.columns.tolist())-set(gen_columns)
    sub_columns = [*meta_data_columns, *gen_columns ]

    win_dict = {k:[] for k in sub_columns}
    print(win_dict)
    for i in range(len(homologs_pairs)):
        idata = homologs_pairs.iloc[i,:][homologs_pairs.columns.tolist()].to_dict()
        aligned_wins, nwins = alignment2wins(idata['query_dna_seq'].split(" "),idata['qseq'],idata['subject_dna_seq'].split(" "),idata['sseq'])
        idata_ = {k:[idata[k]]*nwins for k in meta_data_columns}
        idata_ = {**idata_, **aligned_wins}
        for k in win_dict.keys():
            win_dict[k].extend(idata_[k])
    df = pd.DataFrame.from_dict(win_dict)
    df.to_csv(os.path.join(args.data_path,args.species1+'_'+args.species2+'_'+str(args.partition_id)+'_homologs.nonoverlaping_win'+str(args.win_size)+'.csv'), index=False)

        

if __name__ == "__main__":
    if args.mode=='homologs':
        preprocess_homologs(args.data_path)
    else:
        print("Error: unrecognized mode argument: {}.", args.mode)

