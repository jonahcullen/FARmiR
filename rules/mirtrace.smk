
rule mirtrace_pre_process_csv:
    input:
        fq_join = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_join.corrected.clean.fastq.gz',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            ))
        ),
    output:
        samples_csv = S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/pre_process.csv')
    run:
        # write each fastq and sample name to file
        with open(output.samples_csv, 'w') as f_out:
            for i in input.fq_join:
                print(
                    i, 
                    os.path.basename(i).split('_join')[0],
                    sep=',',file=f_out
                )

rule mirtrace_pre_process:
    input:
        fq_join = S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_join.corrected.clean.fastq.gz',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            )
        ),
        samples_csv = S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/pre_process.csv')
    output:
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-results.json'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-report.html'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-phred.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-qcstatus.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-length.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-rnatype.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-mirna-complexity.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-contamination_basic.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/pre_process/mirtrace-stats-contamination_detailed.tsv'),
    params:
        conda_env  = config['conda_envs']['trace'],
        base       = lambda wildcards, input: os.path.dirname(input.samples_csv)
    threads: 24
    resources:
        time   = 60,
        mem_mb = 240000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e

            mirtrace qc \
                --species eca \
                -c {input.samples_csv} \
                -o MIRTRACE_PRE \
                --force \
                -t {threads}

            cp MIRTRACE_PRE/* {params.base} && rm -rf MIRTRACE_PRE
        '''

rule mirtrace_post_process_csv:
    input:
        fq_join = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mirna.unaligned_mm00.fastq',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
            ))
        ),
    output:
        samples_csv = S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/post_process.csv')
    run:
        # write each fastq and sample name to file
        with open(output.samples_csv, 'w') as f_out:
            for i in input.fq_join:
                print(
                    i, 
                    os.path.basename(i).split('.')[0],
                    sep=',',file=f_out
                )

rule mirtrace_post_process:
    input:
        fq_join = S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mirna.unaligned_mm00.fastq',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
            )
        ),
        samples_csv = S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/post_process.csv')
    output:
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-results.json'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-report.html'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-phred.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-qcstatus.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-length.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-rnatype.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-mirna-complexity.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-contamination_basic.tsv'),
        S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-contamination_detailed.tsv'),
    params:
        conda_env  = config['conda_envs']['trace'],
        base       = lambda wildcards, input: os.path.dirname(input.samples_csv)
    threads: 24
    resources:
        time   = 60,
        mem_mb = 60000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e

            mirtrace qc \
                --species eca \
                -c {input.samples_csv} \
                -o MIRTRACE_POST \
                --force \
                -t {threads}
            
            cp MIRTRACE_POST/* {params.base} && rm -rf MIRTRACE_POST
        '''

# filter samples where <10% of reads were classified as mirnas based on 
# mirtrace above. NOTE that mirtrace uses mirbase 21 by default however
# for the horse there is no difference in mature/hairpin numbers between
# v21 and v22. If a new version is released or a different species is used
# this will need to be reconsidered by generating a custom database set for
# mirtrace
rule distinguish_exp_mirna_samples:
    input:
       #fq_join = S3.remote(
       #    expand(
       #        '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mirna.unaligned_mm00.fastq',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
        cts = S3.remote('{bucket}/private/qc/{release}/small/mirtrace/post_process/mirtrace-stats-rnatype.tsv'),
    output:
        samp_exclude = S3.remote('{bucket}/private/fastq/trimmed/mirna_exp_samples/{release}/rnatype_excluded.tsv'),
        samp_include = S3.remote('{bucket}/private/fastq/trimmed/mirna_exp_samples/{release}/rnatype_included.tsv')
    run:
        # read mirtrace rnatype count table
        df = pd.read_csv(input.cts,sep="\t",index_col='RNA_TYPE')
        # and calculate proportion by rnatype
        props = df.apply(lambda x: x /x.sum(), axis=0)

        # select samples with lt 10% of mirna reads
        exclude = props.loc[:, props.loc['miRNA'] < 0.1]
        exclude.reset_index().melt(
            id_vars='RNA_TYPE', 
            var_name='sample', 
            value_name='prop') \
            .to_csv(output.samp_exclude, sep='\t', index=False)
        # select samples with ge 10% of mirna reads
        included = props.loc[:, props.loc['miRNA'] >= 0.1]
        included.reset_index().melt(
            id_vars='RNA_TYPE', 
            var_name='sample', 
            value_name='prop') \
            .to_csv(output.samp_include, sep='\t', index=False)











