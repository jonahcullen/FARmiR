
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


