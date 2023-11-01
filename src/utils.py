#import os
#import random
#import numpy as np
#import pandas as pd

def get_small_fastqs(wildcards):
    '''
    Get fastq files of given sample-unit.
    '''
    fastqs = units_small.loc[
        (wildcards.sample, wildcards.tissue), 
        ['layout', 'fastq_1', 'fastq_2']
    ]
    # check if sample has paired layout 
    if 'paired' in fastqs.layout:
        return {
            'r1': fastqs.fastq_1, 
            'r2': fastqs.fastq_2
        }
    else: # NOTE setting r2 to the same r1 which is checked in process_small
        return {
            'r1': fastqs.fastq_1,
            'r2': fastqs.fastq_1
        }

def exp_locus(df, loci):
    '''
    for a given loci (e.g. multiple overlapping candidate precursors), determine
    the locus with the largest summed read count. break ties randomly.
    '''
    random.seed(9030) # set random seed
    d = {i:df.loc[i].sum() for i in loci.split('|')}
    max_locus = [k for k, v in d.items() if v == max(d.values())]
    # break ties random if necessary
    if len(max_locus) > 1:
        return random.choice(max_locus)
    else:
        return max_locus[0]

def get_substrings(sequence, size):
    '''
    generate all substrings of a specified size from a starting sequence
    '''
    return [str(sequence[i:i+size]) for i in range(len(sequence) - size + 1)]

def query_alignment(s):
    '''
    from an alignment row in a sam file, extract chrom,
    flag, start, and length
    '''
    # enforced perfect match with bowtie so all cigar strings
    # in column should be NNM
    row = s.split()
    flag = row[1]
    strand = '+'
    seq = Seq(row[9])
    if flag == '16':
        seq = seq.reverse_complement()
        strand = '-'
    seq = str(seq)
    
    return {
        'chrom'  : row[2],
        'strand' : strand,
        'flag'   : flag,
        'start'  : int(row[3]),
        'length' : int(row[5][:-1]),
        'hpins'  : row[0].split('|'),
        'seq'    : seq
    }

def hairpin_contained(hpin_locs, query):
    '''
    checks if the query sequence alignment is contained within at least
    one of the hairpin sequences from which the query sequence substring
    originated
    '''
    # if query chrom is not in any of the hairpin id chroms
    # it is thus not exclusive
    query_chrom = query['chrom']
    hpin_chroms = [hpin_locs[i]['chrom'] for i in query['hpins']]
    if not any(i in query_chrom for i in hpin_chroms):
        return False
    
    # set query range based on alignment flag
    query_range = range(query['start'], query['start'] + query['length'] - 1)
    
    # subset hpin_locs to only include loci with chrom of query
    query_hpins = {i:hpin_locs[i] for i in query['hpins'] if i in hpin_locs}
    
    # next the query must be contained within at least one
    # hairpin loci or else it is not exclusive
    exclusive = any(
        set(query_range).issubset(range(v['start'], v['end'])) 
        for i,v in query_hpins.items()
        if v['strand'] == query['strand']
    )
    
    return exclusive

def replace_snp(s, pos, snp):
    s_list = list(s)
    s_list[pos] = snp
    mod_s = ''.join(s_list)
    return mod_s
    
def snp_contained_kmers(s, s2, k, pos):
    '''
    s is the original hairpin sequence, s2 is the version containing the snp
    '''
    #kmers = []
    d_kmers = {}
    for i in range(len(s) - k + 1):
        if pos >= i and pos < i + k:
            #d_kmers[s[i:i+k]] = s2[i:i+k]
            d_kmers[s2[i:i+k]] = s[i:i+k]
            #kmers.append(kmer)
    return d_kmers

def name_and_hash(s):
    '''
    from a relative file path, return the filename and md5 hash
    '''
    return [os.path.basename(s), hashlib.md5(open(s,'rb').read()).hexdigest()]

def prep_candidate_mirna(novel_mirna):
    '''
        parse mirdeep2 novel mirna output to return filtered dataframe
    '''
    # read mirdeep2 novel mirna file, convert to dna
    novels = []
    with open(novel_mirna,"r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("chr"):
                keep = [line.split("\t")[i] for i in [0,1,4,8,13,15,16]]
                keep[-2] = keep[-2].upper().replace("U","T")
                keep[-3] = keep[-3].upper().replace("U","T")
                novels.append(keep)
            elif line.startswith("mature"):
                break
    # read into dataframe            
    cols = ["prov_id","score","read_count","sig_fold_p","mat_seq","precur_seq","precur_coord"]        
    df = pd.DataFrame(novels,columns = cols)
    # add column with sample name
    df['sample'] = os.path.basename(novel_mirna).split("_small")[0]
    # filter by score, randfold p-value and read count
    df = df.loc[
        (pd.to_numeric(df['score']) >= 1) \
        & (df['sig_fold_p'] == "yes") \
        & (pd.to_numeric(df['read_count']) >= 10)
    ]
    # group by mature and precursor sequences and keep one with max mirdeep2 
    # score
    idx = df.groupby(['mat_seq','precur_coord'])['score'].transform(max) == df['score']
        
    return df[idx]


def generate_isomir_matrix(d,counts=False,seed=20200930):
    '''
    Modified from XICRA (github - )

    Parameters
    ----------
    d : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    all_data = pd.DataFrame()
    seq_all_data = pd.DataFrame()
    iso_filt = []
    iso_all = []
    
    # get keys and shuffle
    keys = list(d.keys())
    random.seed(seed)
    random.shuffle(keys)
  
    # possible isomirs
    isomirs = [
        'iso_3p', 'iso_add3p', 'iso_5p', 'iso_add5p', 'iso_snv',
        'iso_snv_central', 'iso_snv_central_offset', 'iso_snv_central_supp',
        'iso_snv_seed', 'NA'
    ]

    # function extract key from each element
    def sort_key(element):
        return isomirs.index(element.split(':')[0])
    # this was added to deal with situations where the isomir license plates
    # are identical but the actual modifications are just swapped order (e.g.
    # iso_3p:+1,iso_snv vs iso_snv,iso_3p:+1  

    drop = False 
    for k in keys:
        new_data = pd.DataFrame()
        data = pd.read_csv(d[k],sep='\t')
        # set NAs and generate unique id
        data['Variant'].fillna('NA',inplace=True)
       #data['unique_id'] = data.apply(lambda data: data['miRNA'] + '&' + data['Variant'] + '&' + data['UID'], axis=1)
        data['unique_id'] = data.apply(lambda data:
            data['miRNA']
            + '&'
            + ','.join(sorted(data['Variant'].split(','), key=sort_key))
            + '&'
            + data['UID'],
            axis=1
        )
        # miraligner mirtop creates a column containing sample name and other tags (trim, joined, fastq...)
        new_data = data.filter(['unique_id',k],axis=1).set_index('unique_id')
        # sequence information
        seq_data = data.filter(['UID','Read'],axis=1).set_index('UID')
        # seq_data = seq_data.set_index('UID')
        seq_all_data = pd.concat([seq_all_data,seq_data],sort=True).drop_duplicates('Read')
        # all data
        all_data = pd.concat([all_data, new_data], axis=1, sort=True)
        # add number of successive de-duplicated isomirs 
        if not drop:
            iso_filt.append(len(all_data.index))
            drop = True
        else:
            all_data_filt,all_data_dup = discard_UID_duplicated(all_data)
            iso_filt.append(len(all_data_filt.index))
        iso_all.append(len(all_data.index))
    # return matrix and seqs or list of counts
    if not counts:
        return (all_data,seq_all_data)
    else:
        return (iso_filt,iso_all)

def discard_UID_duplicated(df_data, type_res="miRNA"):
    '''
    Copied from XICRA (github - )

    Parameters
    ----------
    d : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    ## get data index
    df_data['ID'] = df_data.index
    new_data = df_data.filter(['ID'], axis=1)   

    # split ID (hsa-let-7a-2-3p&NA&qNkjr6Ov2) into miRNA, variant and UID
    tmp = new_data['ID'].str.split('&', expand = True)
    
    new_data[type_res]  = tmp[0]
    new_data['variant']  = tmp[1]
    new_data['UID']  = tmp[2]
    
    ## count 
    count_groups = new_data.groupby('UID').count()
    ## print to file?
    
    ## get duplicated
    bigger1count = count_groups[ count_groups['ID'] > 1 ]
    ## print counts to file?
    
    ## get list of UIDs duplicate
    bigger1count_list = bigger1count.index.to_list()
    duplicates = new_data[new_data['UID'].isin(bigger1count_list)]
    ## print duplicate to file?

    ## get duplicated data
    duplicates_indes_list = duplicates.index.to_list()
    duplicates_expression = df_data[df_data.index.isin(duplicates_indes_list)]
    #duplicates_expression['UID'] = duplicates_expression['ID'].str.split('&', expand = True)[2]
    duplicates_expression = duplicates_expression.drop(['ID'], axis=1)
    duplicates_expression.index.name = "ID"
    
    ## get clean data
    clean_data_expression = df_data[~df_data.index.isin(duplicates_indes_list)]
    #clean_data_expression['UID'] = clean_data_expression['ID'].str.split('&', expand = True)[2]
    clean_data_expression = clean_data_expression.drop(['ID'], axis=1)
    clean_data_expression.index.name = "ID"
    
    return (clean_data_expression, duplicates_expression)
