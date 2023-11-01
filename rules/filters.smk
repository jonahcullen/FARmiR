import os
import json
import pandas as pd
from scipy import stats

rule filter_samples_reads:
    input:
        pre_jsons = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.pre_trim.fastp.json',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            ))
        ),
        post_jsons = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mm0{mm}.fastp.json',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
                mm=['0']
            ))
        ),
    output:
        samples = S3.remote('{bucket}/private/small/quant/preproc_read_counts.all_samples.tsv', keep_local=True),
       #inc_samples = S3.remote('{bucket}/private/small/quant/preproc_passed_reads.include_samples.tsv'),
       #exc_samples = S3.remote('{bucket}/private/small/quant/preproc_passed_reads.exclude_samples.tsv'),
    run:
        # zip pre and post fastp jsons together to generate df
        l = []
        for i,j in zip(input.pre_jsons, input.post_jsons):
            # ensure tissue_sample identical
            pre_ts = os.path.basename(i).split('.pre_trim')[0]
            post_ts = os.path.basename(j).split('.mm00')[0]
            if pre_ts != post_ts:
                print('ERROR: pre and post jsons no matched!')
                break
            tissue, sample = i.split('/')[4:6]
            # specific to this dataset drop REP sample
            if 'REP' in sample:
                continue
            tissue_sample = f'{tissue}_{sample}'
            with open(i, 'r') as f1_in, open(j, 'r') as f2_in:
                # pre-process read counts
                pre = json.load(f1_in)
                pre_reads = pre['filtering_result']['passed_filter_reads']
                # final post-process read counts
                post = json.load(f2_in)
                post_reads = post['filtering_result']['passed_filter_reads']
            l.append([tissue, sample, tissue_sample, pre_reads, post_reads])

        # create dataframe
        df = pd.DataFrame(l, columns=['tissue', 'sample', 'tissue_sample', 'pre_counts', 'post_counts'])
        
        # calculate the log2 median and median absolute deviation
        df['log2_passed'] = np.log2(df['post_counts'])
        med = df['log2_passed'].median(axis=0)
        # get appropriate constant for mad
        a_n = -0.804168866*len(df)**(-1.008922)
        adj_scale = 1/(stats.norm.ppf(0.75)*(1+a_n))
        print('a_n: ', a_n)
        print('adj_scale: ', adj_scale)
        mad = stats.median_abs_deviation(df['log2_passed'], scale=1/adj_scale)
        print('mad :', mad)
        # determine samples with passed read counts < 3 mads below all med
        df['mad_passed'] = np.where(df['log2_passed'] >= (med - 3*mad), 'PASS', 'FAIL')
       
        # write to outputs
        df.to_csv(output.samples, sep='\t', index=False)
       #inc_samples.to_csv(output.inc_samples, sep='\t', index=False)
       #exc_samples.to_csv(output.exc_samples, sep='\t', index=False)
        
rule filter_samples_biotypes:
    input:
        types = S3.remote(
            expand(
                '{bucket}/private/small/align/bowtie/{release}/all_samples.featcts.csv',
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
            )
        )
    output:
        samples = S3.remote('{bucket}/private/small/quant/preproc_biotype_props.all_samples.tsv'),
    run:
        # read in biotypes, calculate the proportion by biotype for each
        # sample-tissue-pre/post
        tmp = pd.read_csv(input.types[0])
        # specific to this dataset drop REP sample
        df = tmp[~tmp['sample'].str.contains('REP')]
        df['biotype_prop'] = df.groupby(['sample', 'tissue', 'seq_set'])['count'].transform(lambda x: x / x.sum())

        def mirna_fails(grp):
            grp['mirna_passed'] = 'FAIL' if any(grp['biotype'].eq('miRNA') & grp['biotype_prop'].lt(0.1)) else 'PASS'
            return grp

        # label the sample as fail if proportion of mirna < 0.1 and write file
        df2 = df.groupby(['sample', 'tissue', 'seq_set'], group_keys=True).apply(mirna_fails)
        df2.to_csv(output.samples, sep='\t', index=False)
        
rule combine_filters:
    input:
        filt_counts = S3.remote('{bucket}/private/small/quant/preproc_read_counts.all_samples.tsv', keep_local=True),
        filt_types = S3.remote('{bucket}/private/small/quant/preproc_biotype_props.all_samples.tsv', keep_local=True),
    output:
        all_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.all_samples.tsv', keep_local=True),
        inc_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv', keep_local=True)
    run:
        reads_filt = pd.read_csv(input.filt_counts, sep="\t")
        types = pd.read_csv(input.filt_types, sep="\t")
        # add tissue_sample column to match reads_filt and remove REP
        types['tissue_sample'] = types['tissue'] + '_' + types['sample']

        # select need columns
        select_cols = ['tissue_sample', 'mirna_passed']
        filt_types = types[(types['seq_set'] == 'post_filter') & (types['biotype'] == 'miRNA')][select_cols]

        # combine, select PASS, and count samples per tissue
        filters = pd.merge(reads_filt, filt_types)
        
        inc_samples = filters[
            (filters['mad_passed'] == 'PASS') &
            (filters['mirna_passed'] == 'PASS')
        ].copy()

        inc_samples['tissue_count'] = inc_samples.groupby('tissue')['tissue'].transform('count')

        # new column to identify tissues with less than 2 samples remaining
        inc_samples['tissue_passed'] = np.where(inc_samples['tissue_count'] < 2, 'FAIL', 'PASS')

        # select passing samples
        inc_out = inc_samples[
            (inc_samples['mad_passed'] == 'PASS') &
            (inc_samples['mirna_passed'] == 'PASS') &
            (inc_samples['tissue_passed'] == 'PASS')
        ]

        # write to files
        filters.to_csv(output.all_filt, sep='\t', index=False)
        inc_out.to_csv(output.inc_filt, sep='\t', index=False)

