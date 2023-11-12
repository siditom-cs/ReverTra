#!/bin/sh

s1="SCPECBS3"
partition_id=1
s1_seqs="${s1}.${partition_id}.aa.fasta"
s1_dna_seqs="${s1}.${partition_id}.nt.fasta"
log_path="./logs/"
#data_path="./processed_test"
data_path="/home/tomer/CodOpTRM/data/datasets/processed_data_SCPECBS3"
data_path="/home/tomer/COnTRA/data/datasets/processed_data_SCPECBS3"
homologs_path="${data_path}/homologs"


python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path} --win_size=10
python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path} --win_size=30
python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path} --win_size=50
python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path} --win_size=100
python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path} --win_size=150


#python gen_dataset.py --win_size=10
#python gen_dataset.py --win_size=30
#python gen_dataset.py --win_size=50
#python gen_dataset.py --win_size=75
#python gen_dataset.py --win_size=100
#python gen_dataset.py --win_size=150

