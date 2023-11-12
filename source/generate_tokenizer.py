import os 
from Bio.Seq import translate
import json
from transformers import AutoTokenizer
from tokenizers import (
            models,
            pre_tokenizers,
            processors,
            Tokenizer,
            )
from transformers import PreTrainedTokenizerFast

#tokenizer_path="mbart_species_codon_tokenizer"
#tokenizer_path="mbart_translation_FungiDataset2023_codon_tokenizer"
tokenizer_path="contra_tokenizer_gen_exprrefined"

if os.path.exists(tokenizer_path):
    print("Error: Tokenizer dirctory already exists.")
    exit()

def get_dna_vocab():
    DNA_alphabet = ['T', 'C', 'A', 'G']
    codons = []
    for c1 in DNA_alphabet:
        for c2 in DNA_alphabet:
            for c3 in DNA_alphabet:
                codons.append(c1+c2+c3)

    unique_aa = set([translate(c) for c in codons])
    mask_tokens = ['<mask_'+aa+'>' for aa in unique_aa]
    special_tokens = ["<pad>","<gap>","<s>","</s>","<unk>","<cls>","<sep>"]
    expr_tokens = ["<expr_top10>","<expr_pre75_90>","<expr_pre50_75>","<expr_pre25_50>","<expr_low25>","<expr_unk>"]
    seqpos_tokens = ["<seqpos_begin>","<seqpos_end>","<seqpos_mid>"]
    similarity_tokens = ["<sim_high>","<sim_low>"]
    special_tokens = [*special_tokens, *expr_tokens, *seqpos_tokens, *similarity_tokens]
    species_tokens = ["<S_cerevisiae>","<S_pombe>","<E_coli>","<B_subtilis>"]
    #species_tokens = [ "<S_cerevisiae>", "<S_pombe>", "<L_kluyveri>", "<S_castellii>", "<S_kudriavzevii>", "<S_mikatae>", "<S_paradoxus>", "<S_uvarum>", "<Yarrowia_Lipolytica>", "<Ashbya_Gossypii>", "<Aspergillus_Acristatulus>", "<Cryptococcus_Neoformans>", "<Debaryomyces_Hansenii>", "<Fusarium_Graminearum>", "<Kluyveromyces_Lactis>", "<Magnaporthe_Grisea>"]
    return codons, mask_tokens, special_tokens, species_tokens

codons, mask_tokens, special_tokens, species_tokens = get_dna_vocab()

tokenizer = Tokenizer(models.WordLevel(unk_token="<unk>"))
tokenizer.pre_tokenizer = pre_tokenizers.Whitespace()
tokenizer.add_special_tokens(special_tokens)
tokenizer.add_special_tokens(species_tokens)
tokenizer.add_special_tokens(mask_tokens)
tokenizer.add_tokens(codons)
cls_token_id = tokenizer.token_to_id("<cls>")
sep_token_id = tokenizer.token_to_id("<sep>")

mask_translator_dict = {tokenizer.token_to_id(c):tokenizer.token_to_id('<mask_'+translate(c)+'>') for c in codons}
mask_translator_dict = {**mask_translator_dict, **{tokenizer.token_to_id(c):tokenizer.token_to_id(c) for c in mask_tokens}}
mask_translator_dict = {**mask_translator_dict, **{tokenizer.token_to_id(c):tokenizer.token_to_id(c) for c in special_tokens}}
mask_translator_dict = {**mask_translator_dict, **{tokenizer.token_to_id(c):tokenizer.token_to_id(c) for c in species_tokens}}

tokenizer.post_processor = processors.TemplateProcessing(
        single=f"<cls>:0 $A:0 <sep>:0",
        pair=f"$A:0 $B:0",
        special_tokens=[("<cls>", cls_token_id), ("<sep>", sep_token_id)],
        )


wrapped_tokenizer = PreTrainedTokenizerFast(
    tokenizer_object=tokenizer,
    
    bos_token="<s>",
    eos_token="</s>",
    unk_token="<unk>",
    pad_token="<pad>",
    cls_token="<cls>",
    sep_token="<sep>",
    mask_token="<msk>",
    model_max_length=512, 
    padding_side="right"
    )
wrapped_tokenizer.save_pretrained(tokenizer_path)
with open(os.path.join(tokenizer_path,"masking_dict.json"),"w") as f:
    json.dump(mask_translator_dict,f)



tokenizer=wrapped_tokenizer

mask_restriction_dict={}
for mask in mask_tokens:
    mask_restriction_dict[tokenizer.vocab[mask]] = []

for c in codons:
    maskid = tokenizer.vocab["<mask_"+translate(c)+">"]
    #restrict the masked aa only to codons that code it
    mask_restriction_dict[maskid] += [tokenizer.vocab[c]]
    #restrict a real codon only to itself
    mask_restriction_dict[tokenizer.vocab[c]] = [tokenizer.vocab[c]]

for s in special_tokens:
    mask_restriction_dict[tokenizer.vocab[s]] = [tokenizer.vocab[s]]

for s in species_tokens:
    mask_restriction_dict[tokenizer.vocab[s]] = [tokenizer.vocab[s]]

#this makes sure only a single copy in each list (sets)
for k in mask_restriction_dict.keys():
    mask_restriction_dict[k] = list(set(mask_restriction_dict[k]))


with open(os.path.join(tokenizer_path,'mask_restrict_dict.json'), 'w') as f:
    json.dump(mask_restriction_dict, f)



tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)
print(tokenizer)
print(mask_translator_dict)
seq = ["AAA AAG ACG CGC GGG AAA AAA AAA", "AGG ACG ATT"]
species = ["<S_pombe>", "<S_cerevisiae>"]

#print(tokenizer.encode(seq[0]).tokens)
print(wrapped_tokenizer.vocab)
print(wrapped_tokenizer(seq, species, return_tensors ='pt', padding=True).input_ids)
ids = wrapped_tokenizer(seq, species, return_tensors ='pt', padding=True).input_ids
print(wrapped_tokenizer.decode(ids[1]))
        
