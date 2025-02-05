
rule multiqc:
    input: 
        S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.{step}.fastp.json',
                u=units_small.itertuples(),
                bucket=config['bucket'],
                step=['pre_trim','trim','clean','post_join']
            ))
        ),
        S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/all_small_join.post_join.tsv',
                bucket=config['bucket'],
            )
        ),
        S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mm0{mm}.fastp.json',
               #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/noncode_filt/{release}/{u.tissue}_{u.sample}.max_sized.mm0{mm}.fastp.json',
                u=units_small.itertuples(),
                bucket=config['bucket'],
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
                mm='0',
            ))
        ),
        S3.remote(
            sorted(expand(
                '{bucket}/private/small/align/bowtie/filters/{release}/{u.tissue}/{u.sample}/mismatch_00/all/{u.tissue}_{u.sample}.mm00_all.{ref}.bwt.log',
                u=units_small.itertuples(),
                bucket=config['bucket'],
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
                ref=['ncref','ncref_rfam']
            ))
        ),
        S3.remote('{bucket}/private/fastq/trimmed/all_small_join.sized.clean.lengths.csv'),
        S3.remote('{bucket}/private/fastq/trimmed/all_small_join.post_join.tsv'),
    output: 
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_report.html',keep_local=True),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc.log'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_sources.txt'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_data.json'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_fastp.txt'),
    params:
        star_one  = '*STARpass1*',
        conda_env = config['conda_envs']['qc'],
        outdir    = '{bucket}/private/qc/{release}/small'
    threads: 4
    resources:
        time   = 360,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            multiqc {wildcards.bucket}/private/ \
                --interactive \
                --ignore {params.star_one} \
                --force \
                -o {params.outdir}
        '''
