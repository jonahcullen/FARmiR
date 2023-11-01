#!/usr/bin/env python

import sys
from utils import prep_candidate_mirna
import pandas as pd
#from upsetplot import UpSet
#from upsetplot import from_contents
#from matplotlib import pyplot as plt 

# get novel mirna csv files
tissues = sys.argv[1:-1]
print(tissues)

x1 = prep_candidate_mirna(tissues[0])
print(x1.shape)

x1 = prep_candidate_mirna(tissues[0])['mat_seq'].tolist()
print(len(x1))

# get mature mirna sequence lists to dictionary
d = {}
for i in tissues:
    

#spleen = from_contents({'1_10Y': x1_seqs, '2_9Y': x2_seqs, '3_10Y': x3_seqs})

#def prep_candidate_mirna(novel_mirna):
#    '''
#        parse mirdeep2 novel mirna output to return filtered dataframe
#    '''
#    # read mirdeep2 novel mirna file, convert to dna
#    novels = []
#    with open(novel_mirna,"r") as f:
#        for line in f:
#            line = line.strip()
#            if line.startswith("chr"):
#                keep = [line.split("\t")[i] for i in [0,1,4,8,13,16]]
#                keep[-2] = keep[-2].upper().replace("U","T")
#                novels.append(keep)
#            elif line.startswith("mature"):
#                break
#    # read into dataframe            
#    cols = ["prov_id","score","read_count","sig_fold_p","mat_seq","precur_coord"]        
#    df = pd.DataFrame(novels,columns = cols)
#    # filter by score, randfold p-value and read count
#    df = df.loc[
#        (pd.to_numeric(df['score']) >= 1) \
#        & (df['sig_fold_p'] == "yes") \
#        & (pd.to_numeric(df['read_count']) >= 10)
#    ]
#    # group by mature sequence and keep oen with max mirdeep2 score
#    # cand_mirna = df.loc[df.groupby('mat_seq')['score'].idxmax()]
#    idx = df.groupby(['mat_seq'])['score'].transform(max) == df['score']
#        
#    return df[idx]
