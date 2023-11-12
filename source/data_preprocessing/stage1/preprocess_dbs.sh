#!/bin/sh

python preprocess.py --rawdata_path=../../data/raw_data/SCPECBS \
	--data_path=../../data/datasets/processed_data_SCPECBS3 \
	--log_path=../..//data/logs/SCPECBS3_logs/ \
	--agg_seqs_path=../..//data/logs/SCPECBS3_logs/agg_seqs.temp.fasta \
	--cdhit_log=../../logs/SCPECBS3_logs/cdhit.log 
	
