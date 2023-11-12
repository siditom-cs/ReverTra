import os
import random
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    seq_rec = SeqRecord(
        seq=nuc_record.seq.translate(),
        id="trans_" + nuc_record.id,
        description="translation", 
    )
    seq_rec.species = nuc_record.species
    return seq_rec


def check_refdbs_partition(cls_members, ref_dbs):
    is_cls_in_prev_partition=False
    inx = -1
    for sid in cls_members:
        for ref_db in ref_dbs:
            for i in range(3):
                if sid in ref_db[str(i)]:
                    inx = i
                    is_cls_in_prev_partition = True
                    return is_cls_in_prev_partition, inx
    return is_cls_in_prev_partition, inx

def split_data_by_previous_partition(cls_df, ref_dbs, prec=[0.7,0.1,0.2]):
    cls_inx = cls_df[0].unique() 
    random.shuffle(cls_inx)

    tnos = len(cls_df) # total number of sequences

    nos = [0,0,0]
    partition = [list(),list(),list()]
    for icls in range(0,len(cls_inx)):
        noe_cls = sum(cls_df[0]==cls_inx[icls])
        cls_members = cls_df[2][cls_df[0]==cls_inx[icls]].values
        is_cls_in_prev_partition, ipartition = check_refdbs_partition(cls_members, ref_dbs)
        if is_cls_in_prev_partition:
            partition[ipartition].append(cls_inx[icls])
            nos[ipartition] = nos[ipartition] + noe_cls
        elif (nos[0]/tnos <prec[0]):
          partition[0].append(cls_inx[icls])
          nos[0] = nos[0] + noe_cls
        elif (nos[1]/tnos <prec[1]):
          partition[1].append(cls_inx[icls])
          nos[1] = nos[1] + noe_cls
        else:
          partition[2].append(cls_inx[icls])
          nos[2] = nos[2] + noe_cls

    return nos, partition  

def split_data_by_cls(cls_df, prec=[0.7,0.1,0.2]):
    cls_inx = cls_df[0].unique() 
    random.shuffle(cls_inx)

    tnos = len(cls_df) # total number of sequences

    nos = [0,0,0]
    partition = [list(),list(),list()]
    for icls in range(0,len(cls_inx)):
        noe_cls = sum(cls_df[0]==cls_inx[icls])
        if (nos[0]/tnos <prec[0]):
          partition[0].append(cls_inx[icls])
          nos[0] = nos[0] + noe_cls
        elif (nos[1]/tnos <prec[1]):
          partition[1].append(cls_inx[icls])
          nos[1] = nos[1] + noe_cls
        else:
          partition[2].append(cls_inx[icls])
          nos[2] = nos[2] + noe_cls

    return nos, partition  

def aggregate_species_records(seq_dbs):
    seq_dict = {}
    trans_records = []
    for i in range(len(seq_dbs)):
        records = list(SeqIO.parse(seq_dbs[i]['path'], "fasta")) 
        records = list(filter(lambda x: ((len(x)>3) & ((len(x)%3) == 0)), records))
        for rec in records:
            rec.species = seq_dbs[i]["species"]

        seq_dict.update({rec.id : rec for rec in records})
        trans_records.extend([make_protein_record(nuc_rec) for nuc_rec in records])
    return seq_dict, trans_records

def process_data(seq_dict, cls_df, partition, seq_dbs):
    #init
    processed_data = {}
    for i in range(len(partition)):
        processed_data[i] = dict()
        for j in range(len(seq_dbs)):
            species = seq_dbs[j]['species']
            processed_data[i][species] = list()

    for i in range(len(cls_df)):
        rec_id = cls_df.iloc[i,2]
        rec_cls_id = cls_df.iloc[i,0]
        partition_id = np.where([rec_cls_id in set_i for set_i in partition])[0][0]
        species = seq_dict[rec_id].species
        processed_data[partition_id][species].append(seq_dict[rec_id])

    return processed_data

def save_files(processed_data, data_path):
    for partition_id in processed_data:
        for species in processed_data[partition_id]:
            file_name = species+"."+str(partition_id)
            
            data_dir = os.path.join(data_path,species)
            if not os.path.exists(data_dir):
                os.mkdir(data_dir)

            dst_path = os.path.join(data_dir,file_name)

            SeqIO.write(processed_data[partition_id][species],dst_path+'.nt.fasta','fasta')

            SeqIO.write((make_protein_record(rec) for rec in processed_data[partition_id][species]), dst_path+".aa.fasta",'fasta')


def combine_fasta_files(data_path, species_name):
    data_path = os.path.join(data_path, species_name)
    train = os.path.join(data_path, species_name+".0.nt.fasta")
    val = os.path.join(data_path, species_name+".1.nt.fasta")
    trainval = os.path.join(data_path, species_name+".3.nt.fasta")
    records_train = list(SeqIO.parse(train, "fasta")) 
    records_val = list(SeqIO.parse(val, "fasta")) 
    records = [*records_train, *records_val]
    SeqIO.write(records,trainval,'fasta')


def create_subset_with_expr(records, expr_dict):
    for rec in records:
        print(expr_dict[rec.id])

def add_expr_level(data_path, species_name, pa_path):
    expr_df = pd.read_csv(pa_path)
    expr_df = {expr_df.iloc[i]['ids']:expr_df.iloc[i]['expr'] for i in range(len(expr_df))}
    data_path = os.path.join(data_path, species_name)
    train = os.path.join(data_path, species_name+".0.nt.fasta")
    val = os.path.join(data_path, species_name+".1.nt.fasta")
    test = os.path.join(data_path, species_name+".2.nt.fasta")
     
    records_train = list(SeqIO.parse(train, "fasta")) 
    create_subset_with_expr(records_train, expr_df)
    records_val = list(SeqIO.parse(val, "fasta")) 
    create_subset_with_expr(records_val, expr_df)
    records_test = list(SeqIO.parse(test, "fasta")) 
    create_subset_with_expr(records_test, expr_df)

