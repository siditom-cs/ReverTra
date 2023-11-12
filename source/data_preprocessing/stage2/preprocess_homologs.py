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


def generate_dictionaries(data_path, partition_id):
    protO1 = os.path.join(data_path,args.species1+'.'+str(partition_id)+'.aa.fasta')
    protO2 = os.path.join(data_path,args.species2+'.'+str(partition_id)+'.aa.fasta')
    dnaO1 = os.path.join(data_path,args.species1+'.'+str(partition_id)+'.nt.fasta')
    dnaO2 = os.path.join(data_path,args.species2+'.'+str(partition_id)+'.nt.fasta')

    seq_dict1 = {rec.id.replace('trans_','') : rec.seq for rec in SeqIO.parse(protO1, "fasta")}
    seq_dict2 = {rec.id.replace('trans_','') : rec.seq for rec in SeqIO.parse(protO2, "fasta")}
    seq_dict_aa = {**seq_dict1, **seq_dict2}

    seq_dict1 = {rec.id : rec.seq for rec in SeqIO.parse(dnaO1, "fasta")}
    seq_dict2 = {rec.id : rec.seq for rec in SeqIO.parse(dnaO2, "fasta")}
    seq_dict_dna = {**seq_dict1, **seq_dict2}
    return seq_dict_aa, seq_dict_dna



def get_codons_for_aa(dna1, s1r, dna2, s2r):
    """
        Gets the alignment and dna seqs and returns the aligned codon sequences.
        Note: 
            1. Gaps in sequences 1 - are replaced with <mask_AA>, where the AA is the amino-acid in the subject (aka. the target) sequence.
            2. AA/Codons in sequence 1 with gaps in sequence 2 are removed from the alignment.
            3. Sequence 1 indices, with different AA from sequence 2, are replaced with <mask_AA>, where the AA of sequence 2 from the same index.
            4. As in 2, gaps in sequence 2 are removed from the codons alignment.
    """
    sp1 = dna1.translate().find(s1r.replace('-',''))
    #print('dna1:    ', dna1.translate()[(sp1):(sp1+len(s1r))])
    #print('s1r:     ', s1r)
    codons1 = []
    di=0
    for i in range(len(s1r)):
        if not s1r[i] == '-':
            codons1.append(str(dna1[3*(sp1+di):3*(sp1+di+1)]))
            di+=1
        else:
            #codons1.append('<mask_'+str(s2r[i])+'>')
            codons1.append('<gap>')
    #print('codons1: ', codons1.translate())
    sp2 = dna2.translate().find(s2r.replace('-',''))
    codons2 = []
    di=0
    for i in range(len(s2r)):
        if not s2r[i] == '-':
            codons2.append(str(dna2[3*(sp2+di):3*(sp2+di+1)]))
            di+=1
        else:
            codons2.append('<mask_'+str(s1r[i])+'>')

        #else:
        #    codons2 += '---'
    #print('dna2:    ', dna2.translate()[(sp2):(sp2+len(s2r))])
    #print('s2r:     ', s2r)
    #print('codons2: ', codons2.translate())
    assert(len(codons1)==len(codons2))
    codons1 = ' '.join(codons1)
    codons2 = ' '.join(codons2)
    return str(codons1), str(codons2)

def init_species_dict(species=['B_subtilis','E_coli','S_cerevisiae','S_pombe']):
    species_db_path=os.path.join(args.data_path,'..')
    d = {s:[] for s in species}
    for s in species:
        species_path = os.path.join(species_db_path, s)
        species_path = os.path.join(species_path, s+"."+str(args.partition_id)+'.nt.fasta')
        for rec in SeqIO.parse(species_path, "fasta"):
            d[s].append(rec.id)

        print(species_path)
    print(d)
    return d

def get_id_origin_species(rec_id, species_dict):
    for s in species_dict.keys():
        if rec_id in species_dict[s]:
            return s
    return False

def preprocess_homologs(data_path, partition_id, sw_size, sw_jump):
    seq_dict_aa, seq_dict_dna = generate_dictionaries(data_path, partition_id)
    blast_res_path= os.path.join(args.log_path,'blast_'+args.species1+'_'+args.species2+'_'+str(args.partition_id)+'_tbl')
    aligned_sequences = pd.read_csv(blast_res_path, delimiter='\t', header=None)
    aligned_sequences.columns =['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq']
    aligned_sequences['qseqid'] = [s.replace('trans_','') for s in aligned_sequences['qseqid']]
    aligned_sequences['sseqid'] = [s.replace('trans_','') for s in aligned_sequences['sseqid']]

    print(aligned_sequences.head())
    homologs_outfile=os.path.join(args.data_path,args.species1+'_'+args.species2+'_'+str(args.partition_id)+'_homologs.csv')


    print("Number of aligned sequences: ",len(aligned_sequences))
    print("Homologs windows output file: ", homologs_outfile)
    query_dna_seqs = []
    subject_dna_seqs = []
    query_species = []
    subject_species = []
    species_dict = init_species_dict()
    for i in range(len(aligned_sequences)):
        name1 = aligned_sequences.iloc[i]['qseqid']
        name2 = aligned_sequences.iloc[i]['sseqid']
        s1 = aligned_sequences.iloc[i]['qseq']
        s2 = aligned_sequences.iloc[i]['sseq']
        query_dna,subject_dna = get_codons_for_aa(seq_dict_dna[name1], s1, seq_dict_dna[name2], s2)
        query_dna_seqs.append(query_dna)
        subject_dna_seqs.append(subject_dna)
        query_species.append(get_id_origin_species(name1, species_dict))
        subject_species.append(get_id_origin_species(name2, species_dict))

    aligned_sequences['query_dna_seq']=query_dna_seqs
    aligned_sequences['subject_dna_seq']=subject_dna_seqs
    aligned_sequences['query_species'] = query_species
    aligned_sequences['subject_species'] = subject_species
    print(aligned_sequences.head())
    aligned_sequences.to_csv(homologs_outfile, index=False)    

if __name__ == "__main__":
    if args.mode=='homologs':
        preprocess_homologs(args.data_path, args.partition_id, args.win_size, args.sw_jump)
    else:
        print("Error: unrecognized mode argument: {}.", args.mode)

