#!/bin/sh
#CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.eval --config_path=./source/evaluation/evaluate.py  
CUDA_VISIBLE_DEVICES=2 python -m source.evaluation.evaluate --config_path=./source/evaluation/models_configs/COBaBEx/hrm_config.json 
#CUDA_VISIBLE_DEVICES=0 python -m source.evaluation.evaluate --config_path=./source/evaluation/models_configs/COBaBEx/finetune2_homologs.json 

