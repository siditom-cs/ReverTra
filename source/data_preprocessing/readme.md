### Data Processing - README####
Here we provide a description about the data processing scripts.

Note1: the outputs of the preprocessing are provided in this project - except for the stage4 that generates the Databases for the training scripts.
Note2: Each stage scripts should be activated from its folder.

stage1 - download and partition data:
------------------------------------
Prerequisites: CD-Hit installed.

1. download_orfs.sh - (Optional; the raw data files are provided) download the mRNA sequences of the different species.  
2. preprocess_dbs.sh - partition with cd-hit


stage2 - blast: all-against-all for each set:
---------------------------------------------
Prerequisites: blast installed (specifically, blastp  and makeblastdb).

1. Blast all vs all for traning, validation, and test sets (0,1,and 2, respectively) - apply in stage2 folder:
	./calculate_homologs.sh 0  
	./calculate_homologs.sh 1
	./calculate_homologs.sh 2

2. Create validation fixed-window dataset:
	./calculate_win_homologs.sh

stage3 - combine sequence and blast data with expression level data:
-------------------------------------------------------------------
1. Add categorized expression level data to the .csv data files generated in stage2:
	./python add_expr2datasets.py 

stage4 - generate arrow database for training:
---------------------------------------------
1. Generate Arrow databases for each window size (they differ only in the validation set):
	./gen_DBs.sh


cai_option:
-----------	
Prepares subsets of fasta files from each species accordin to expression level thresholds. These subsets are then used to create difference frequency models.

