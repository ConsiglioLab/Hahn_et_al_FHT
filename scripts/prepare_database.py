# load .py file with dictionary
# "data/GSEA_sets/BTM_plus.py"
# and save it as .tsv file

import os
import pandas as pd
import numpy as np
import pickle


# load dictionary
with open('../data/GSEA_sets/BTM_plus.py', 'rb') as f:
    BTM_Plus = f.read()
    exec(BTM_Plus)

# find the maximum length of the 'genes' lists
max_len = max(len(gene_set['genes']) for gene_set in BTM_Plus)

# create a DataFrame, padding shorter lists with None
df = pd.DataFrame({gene_set['name']: gene_set['genes'] + [None] * (max_len - len(gene_set['genes'])) for gene_set in BTM_Plus})

# save to tsv
df.to_csv('../data/GSEA_sets/BTM_plus.tsv', sep='\t', index=False)



