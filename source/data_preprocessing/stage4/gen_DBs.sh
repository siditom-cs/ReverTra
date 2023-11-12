#!/bin/sh

python gen_dataset.py --win_size=10 
python gen_dataset.py --win_size=30
python gen_dataset.py --win_size=50
python gen_dataset.py --win_size=75
python gen_dataset.py --win_size=100
python gen_dataset.py --win_size=150

