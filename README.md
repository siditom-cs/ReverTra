# ReverTra

This project corresponds to the paper _"Predicting gene sequences with AI to study evolutionarily
selected codon usage patterns" by Tomer Sidi, Shir Bahiri-Elitzur, Tamir Tuller, and Rachel Kolodny._

Here you can find the code, raw data, and additional information needed to reproduce the data processing, model training and model evaluation of the project. 
- The trained models, processed datasets, and tokenizer can be found at huggingface hub, at the link: https://huggingface.co/siditom.
- Inference examples and downloading the benchmark dataset used in the research can be found at the **notebooks** folder.
- A demo and codon optimization tool based on ReverTra can be found here: https://www.aa2codons.info/
- Note that the models and the code were generated and executed with a linux ubuntu 22.04 server with A100 GPUs. Nonetheless, the trained transformer models have ~6M parameters, which can potentially be trained on less advanced systems.

## Project structure: 

1. notebooks: folder for jupyter notebooks for downloading finetuned models and using them for inference (input - amino-acid sequence; output - species specific codon sequence).
2. tokenizers: folder of custom huggingface tokenizers for codon sequence and species identifiers.
3. configs: folder with the configuration files that were used for pre-training and finetuning the models.
4. data:
   - The raw data (mrna sequences of S.cerevisiae, S.pombe, E.coli, B.subtilis) with protein expression level files.
   - The sequence partition (by identifiers) for training/validation/test sets.
   - The prediction results on the test set for all trained models.
5. source:
   - data preprocessing: code to recreate the data partition and alignments. Dependancies: CD-Hit, Blast.
   - models: contains three folders, one for each model. Frequency model (CAI), Harmonization model, and Transformers models.
6. At the main folder there are scripts for training and evaluation of the models.

### Model training and evaluation:

For training the model, use the script _train_model.sh_.  <br>
The script includes three consecutive training runs (pre-training, finetuning, and second finetuning) and depicted in the paper. The configuration files for each run is in the configs folder. The default training window is 75, for a different one - edit the script.
<br><br>
For model evaluation, use the script _eval_model.sh_. <br>
The evaluation options are configured using the config files:
- ***eval_mimic.json*** : for models with two finetuning runs (mimic).
- ***eval_sis.json*** : for models with a single finetuning run (mask).
- ***eval2_mask.json***: for models with two finetuning runs, yet with mask inference.

Evaluation of the frequency model: ***eval_cais.sh***



