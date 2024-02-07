# ReverTra


Tool based on ReverTra can be found here: https://www.aa2codons.info/
This project contains the code, raw data, and additional information to reproduce the training and evaluation of the ReverTra project.
The project contains the following:
1. notebooks - folder for jupyter notebooks for downloading finetuned models and using them for inference (input: amino-acid sequence; output: species specific codon sequence).
2. tokenizers - folder of custom huggingface tokenizers for codon sequence and species identifiers.
3. configs - folder with the configuration files that were used for pre-training and finetuning the models.
4. data:
   a. The raw data (mrna sequences of S.cerevisiae, S.pombe, E.coli, B.subtilis) with protein expression level files.
   b. The sequence partition (by identifiers) for training/validation/test sets.
   c. The prediction results on the test set for all trained models.
5. source:
   a. data preprocessing
   b. models
   b. 
## Assets in this project:
1. Script to produce a codon tokenizer
