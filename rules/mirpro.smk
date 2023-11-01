
rule novoindex_genome:
    input:
        S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna'),
    output:
        S3.remote('{bucket}/public/refgen/{release}/NOVOINDEX/{release}_genomic.nice.idx')
    params:
        conda_env = config['conda_envs']['novel_small']
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            novoindex {output} {input}            
        '''

# NOTE that miFam.dat is from v21 of miRbase as there does not appear to be
# a family file for v22
rule mirpro_quant:
    input:
        sized = S3.remote(
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
        hairpin  = S3.remote('{bucket}/public/mirbase/v22/hairpin.fa',keep_local=True),
        mature   = S3.remote('{bucket}/public/mirbase/v22/mature.fa',keep_local=True),
        family   = S3.remote('{bucket}/public/mirbase/v21/miFam.dat',keep_local=True),
        nice_fna = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna',keep_local=True),
        nice_idx = S3.remote('{bucket}/public/refgen/{release}/NOVOINDEX/{release}_genomic.nice.idx',keep_local=True),
        inc_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv', keep_local=True)
    output:
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/3_other_form.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/5_other_form.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/arm.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/mapping_report.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_mature.fa'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_precursor.fa'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_precursor.str'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_quantifier_report.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/precursor_mature_relation.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/quantifier_report.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_family.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_mature.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_precursor.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_novel_mature.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_novel_precursor.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/statistics.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/statistics_collapsed.csv'),
        S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/variation_report.csv'),
    params:
        conda_env  = config['conda_envs']['novel_small'],
        outdir     = '{bucket}/private/small/quant/mirpro/{release}/multi_sample',
        filt_sized = lambda wildcards, input: ' -i '.join(
            [
                i for i in input.sized 
                if os.path.basename(i).split('.mirna')[0] in pd.read_csv(input.inc_filt, sep='\t')['tissue_sample'].tolist()
            ]
        ),
    threads: 64
    resources:
        time   = 1440,
        mem_mb = 450000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
       
            mirpro \
                -t {threads} \
                -i {params.filt_sized} \
                -m {input.mature} \
                -p {input.hairpin} \
                -f {input.family} \
                -s eca \
                -a null \
                --strand yes \
                --novel 1 \
                -g {input.nice_fna} \
                --index {input.nice_idx} \
                --other hsa --other mmu --other cfa --other bta --other ssc \
                -d {params.outdir}
        '''

