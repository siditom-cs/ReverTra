
### Download - ORF sequences

dir_path="../../data/raw_data"

# Download - protein abundance data - origin wevbsite: PAXdb

wget https://pax-db.org/downloads/4.1/datasets/4932/4932-WHOLE_ORGANISM-integrated.txt -O $dir_path/S_cerevisiae.PA.dat
wget https://pax-db.org/downloads/4.1/datasets/4896/4896-WHOLE_ORGANISM-integrated.txt -O $dir_path/S_pombe.PA.dat

#SGD
# Download - mRNA sequences of s_cerevisiae from the oficial research site of s_cereviciae: http://sgd-archive.yeastgenome.org/
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz -O $dir_path/S_cerevisiae.fasta.gz
gunzip $dir_path/S_cerevisiae.fasta.gz

#S_Pombe
# Download - mRNA sequences of s_pombe from the oficial research site of s_pombe: https://www.pombase.org/
wget https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz -O $dir_path/S_pombe.fasta.gz
gunzip $dir_path/S_pombe.fasta.gz




#SGD - Optional
# Download - additional options for yeast species from the oficial research site of s_cereviciae: http://sgd-archive.yeastgenome.org/
wget http://sgd-archive.yeastgenome.org/sequence/fungi/S_castellii/archive/WashU/orf_dna/orf_genomic.fasta.gz -O $dir_path/S_castellii.fasta.gz
gunzip $dir_path/S_castellii.fasta.gz

wget http://sgd-archive.yeastgenome.org/sequence/fungi/S_kudriavzevii/archive/WashU/orf_dna/orf_genomic.fasta.gz -O $dir_path/S_kudriavzevii.fasta.gz
gunzip $dir_path/S_kudriavzevii.fasta.gz

wget http://sgd-archive.yeastgenome.org/sequence/fungi/S_mikatae/archive/WashU/orf_dna/orf_genomic.fasta.gz -O $dir_path/S_mikatae.fasta.gz
gunzip $dir_path/S_mikatae.fasta.gz

wget http://sgd-archive.yeastgenome.org/sequence/fungi/S_paradoxus/archive/MIT/orf_dna/orf_genomic.fasta.gz -O $dir_path/S_paradoxus.fasta.gz
gunzip $dir_path/S_paradoxus.fasta.gz

wget http://sgd-archive.yeastgenome.org/sequence/fungi/S_uvarum/archive/WashU/orf_dna/orf_genomic.fasta.gz -O $dir_path/S_uvarum.fasta.gz
gunzip $dir_path/S_uvarum.fasta.gz

wget http://sgd-archive.yeastgenome.org/sequence/fungi/L_kluyveri/archive/WashU/orf_dna/orf_genomic.fasta.gz -O $dir_path/L_kluyveri.fasta.gz
gunzip $dir_path/L_kluyveri.fasta.gz



