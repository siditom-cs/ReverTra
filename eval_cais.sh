#!/bin/bash

#OLD: #CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path="./source/evaluation/cai_configs/$s/CAI.${th}.${s}.json" --out_filename="${s}.${th}.csv"
s="S_cerevisiae"
s="S_pombe"
s="E_coli"
s="B_subtilis"
#th="0.05"
#CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path="./source/evaluation/cai_configs/$s/CAI.${th}.${s}.json"
#mv CAI_Evals/res.csv CAI_Evals/${s}_${th}.csv
th="0.1"
CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path="./source/evaluation/cai_configs/$s/CAI.${th}.${s}.json"
mv CAI_Evals/res.csv CAI_Evals/${s}_${th}.csv
#th="0.15"
#CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path="./source/evaluation/cai_configs/$s/CAI.${th}.${s}.json"
mv CAI_Evals/res.csv CAI_Evals/${s}_${th}.csv
CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path="./source/evaluation/cai_configs/$s/CAI.all.${s}.json"
mv CAI_Evals/res.csv CAI_Evals/${s}_all.csv
