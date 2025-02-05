
rule bowtie_index_ref:
    input:
        S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.fna", keep_local=True)
    output:
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.1.ebwt",keep_local=True),
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.2.ebwt",keep_local=True),
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.3.ebwt",keep_local=True),
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.4.ebwt",keep_local=True),
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.rev.1.ebwt",keep_local=True),
        S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna.rev.2.ebwt",keep_local=True),
        bwt_fa  = S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna",keep_local=True),
        idx_dir = directory("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC")
    params:
        conda_env = config['conda_envs']['small'],
    threads: 1
    resources:
        time   = 360,
        mem_mb = 24000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            cp {input} {output.idx_dir}
            bowtie-build {output.bwt_fa} {output.bwt_fa}
        '''

rule convert_fastq_to_fasta:
    input:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq.gz'),
    output:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.pre_filter.clean.fasta'),
    params:
        conda_env = config['conda_envs']['quant']
    threads: 1
    resources:
        time   = 60,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            seqtk seq -A {input} > {output}
            # remove whitespace from headers
            sed 's, ,-,g' -i {output}
        '''

rule convert_nc_filt_to_fasta:
    input:
        S3.remote(f"{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config['release'][0]}/{{tissue}}_{{sample}}.mirna.unaligned_mm00.fastq"),
    output:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.post_filter.clean.fasta'),
    params:
        conda_env = config['conda_envs']['quant']
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            seqtk seq -A {input} > {output}
            # remove whitespace from headers
            sed 's, ,-,g' -i {output}
        '''

rule small_align_ref:
    input:
        bt_idx   = rules.bowtie_index_ref.output.idx_dir,
        fna_idx  = S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna',keep_local=True),
        small_fa = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.{seq_set}.clean.fasta'),
    output:
        sam    = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.sam'),
        unmap  = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.unaligned.fasta'),
        bt_log = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.bwt.log'),
    params:
        conda_env = config['conda_envs']['small']
    threads: 12
    resources:
        time   = 1440,
        mem_mb = 36000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -f \
                -n {wildcards.mm} \
                -e 70 \
                -l 8 \
                -a \
                -m 6 \
                --best \
                --strata \
                -p {threads} \
                {input.fna_idx} \
                {input.small_fa} \
                -S {output.sam} \
                --un {output.unmap} \
                2> {output.bt_log}
        '''

# individual sample feature counts
rule feature_counts_ref:
    input:
        sam     = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.sam'),
        ref_gtf = S3.remote('{bucket}/public/refgen/{release}/Equus_caballus.EquCab3.0.103_genomic.nice.gtf', keep_local=True),
    output:
        ct_mat = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.featcts.txt'),
        fc_sum = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.featcts.txt.summary'),
        fc_log = S3.remote('{bucket}/private/small/align/bowtie/{release}/{tissue}/{sample}/mismatch_0{mm}/{seq_set}/{tissue}_{sample}.mm0{mm}_{seq_set}.ref.featcts.log'),
    params:
        conda_env = config['conda_envs']['small'],
        prefix    = lambda wildcards, output: os.path.splitext(output.ct_mat)[0],
    threads: 2
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            featureCounts \
                -T {threads} \
                -s 1 \
                -M \
                -t exon \
                -g transcript_biotype \
                -a {input.ref_gtf} \
                -o {output.ct_mat} \
                {input.sam} \
                2> {output.fc_log}
        '''

rule gather_featcts_biotypes:
    input:
        ct_mat = S3.remote(
            expand(
                '{bucket}/private/small/align/bowtie/{release}/{u.tissue}/{u.sample}/mismatch_0{mm}/{seq_set}/{u.tissue}_{u.sample}.mm0{mm}_{seq_set}.ref.featcts.txt',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                    'Equus_caballus.EquCab3.0.103'
                ],
                mm='0',
                seq_set=['pre_filter','post_filter']
            )
        )
    output:
        types = S3.remote('{bucket}/private/small/align/bowtie/{release}/all_samples.featcts.csv', keep_local=True)
    threads: 1
    resources:
        time   = 30,
        mem_mb = 6000
    run:
        # extract biotype counts per mismatch/set/tissue and write to file
        with open(output.types, 'w') as f_out:
            print('sample,tissue,mismatch,seq_set,biotype,count',file=f_out)
            # use list(set) due to the way the input expand works
            for i in list(set(input.ct_mat)):
                # get mismatch/set/tissue from file name
                tissue = i.split('/')[6]
                sample = i.split('/')[7]
                mis,seq_set = os.path.basename(i).split('.')[1].split('_',1)
                with open(i, 'r') as f_in:
                    for line in f_in:
                        if line.startswith(tuple(['#', 'Geneid'])):
                            continue
                        else:
                            biotype = line.split('\t')[0]
                            count = int(line.split('\t')[-1])
                            # write to file
                            print(
                                sample, tissue, mis, seq_set, biotype, count, 
                                sep=',', file=f_out
                            )

#################################################
# combine joined fastas, collapse, map, and
# count features by tissue
#################################################
rule combine_pre_filter_x_tissue:
    input:
        fastas = S3.remote(
            expand(
                "{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/{u.tissue}_{u.sample}_join.pre_filter.clean.fasta",
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            )
        ),
    output:
        all_fa  = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.pre_filter.clean.fasta"),
    params:
        tissue_fa = lambda wildcards, input: " ".join(
            set([i for i in input.fastas if wildcards.tissue in i])
        ),
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            # combine fastas
            cat {params.tissue_fa} > {output.all_fa}
        '''

rule combine_post_filter_x_tissue:
    input:
        fastas = S3.remote(
            expand(
                "{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/{u.tissue}_{u.sample}_join.post_filter.clean.fasta",
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            )
        ),
    output:
        post_fa  = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.post_filter.clean.fasta"),
    params:
        tissue_fa = lambda wildcards, input: " ".join(
            set([i for i in input.fastas if wildcards.tissue in i])
        ),
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            # combine fastas
            cat {params.tissue_fa} > {output.post_fa}
        '''

rule small_align_tissue_ref:
    input:
        bt_idx    = rules.bowtie_index_ref.output.idx_dir,
        fna_idx   = S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna",keep_local=True),
        tissue_fa = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.{seq_set}.clean.fasta"),
    output:
        sam    = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.sam"),
        unmap  = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.unaligned.fasta"),
        bt_log = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.bwt.log"),
    params:
        conda_env = config['conda_envs']['small']
    threads: 12
    resources:
        time   = 240,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -f \
                -n {wildcards.mm} \
                -e 70 \
                -l 8 \
                -a \
                -m 10 \
                --best \
                --strata \
                -p {threads} \
                {input.fna_idx} \
                {input.tissue_fa} \
                -S {output.sam} \
                --un {output.unmap} \
                2> {output.bt_log}
        '''

rule feature_counts_tissue_ref:
    input:
        sam     = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.sam"),
        ref_gtf = S3.remote('{bucket}/public/refgen/{release}/Equus_caballus.EquCab3.0.103_genomic.nice.gtf', keep_local=True),
    output:
        ct_mat = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.featcts.txt"),
        fc_sum = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.featcts.txt.summary"),
        fc_log = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_{seq_set}.ref.featcts.log"),
    params:
        conda_env = config['conda_envs']['small'],
        prefix    = lambda wildcards, output: os.path.splitext(output.ct_mat)[0],
    threads: 2
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            featureCounts \
                -T {threads} \
                -s 1 \
                -M \
                -t exon \
                -g transcript_biotype \
                -a {input.ref_gtf} \
                -o {output.ct_mat} \
                {input.sam} \
                2> {output.fc_log}
        '''

