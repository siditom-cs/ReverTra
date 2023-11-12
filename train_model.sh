#!/bin/sh
winsize=75
CUDA_VISIBLE_DEVICES=0 python -m source.models.COBaBExRi.train --config_path=./configs/wins_configs/win${winsize}/Pretrain_config.json
CUDA_VISIBLE_DEVICES=0 python -m source.models.COBaBExRi.train --config_path=./configs/wins_configs/win${winsize}/Finetune_oSingleSpecies.json
CUDA_VISIBLE_DEVICES=0 python -m source.models.COBaBExRi.train --config_path=./configs/wins_configs/win${winsize}/Finetune2_homologs.json
