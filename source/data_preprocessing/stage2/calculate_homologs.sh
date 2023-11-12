#!/bin/sh

s1="SCPECBS3"
partition_id=$1
s1_seqs="${s1}.${partition_id}.aa.fasta"
s1_dna_seqs="${s1}.${partition_id}.nt.fasta"
log_path="./logs/"
#data_path="./processed_test"
data_path="/home/tomer/CodOpTRM/data/datasets/processed_data_SCPECBS3"

mkdir "${data_path}/homologs"

homologs_path="${data_path}/homologs"
mkdir ${homologs_path}

cp ${data_path}/${s1}/${s1_seqs} ${homologs_path}
cp ${data_path}/${s1}/${s1_dna_seqs} ${homologs_path}

makeblastdb -in  ${homologs_path}/${s1_seqs} -parse_seqids  -dbtype prot

blastp -query ${homologs_path}/${s1_seqs} -db ${homologs_path}/${s1_seqs} -evalue 0.01 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" > ${homologs_path}/blast_${s1}_${s1}_${partition_id}_tbl


python preprocess_homologs.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path}
#python preprocess2wins.py --partition_id=${partition_id} --species1=${s1} --species2=${s1} --data_path=${homologs_path} --log_path=${homologs_path}
rm ${homologs_path}/blast_${s1}_${partition_id}_tbl

rm ${homologs_path}/*.fasta

