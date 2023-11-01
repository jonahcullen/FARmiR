
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

#rule seqkit_dedup_fastq:
#    input:
#        sized = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/small/{sample}/{tissue}_{sample}_{lane}_join.sized.clean.fastq.gz"),
#    output:
#        dedup = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/small/{sample}/{tissue}_{sample}_{lane}_join.dedup.sized.clean.fastq.gz"),
#        dupes = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/small/{sample}/{tissue}_{sample}_{lane}_join.dupes.sized.clean.fastq.gz"),
#        deets = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/small/{sample}/{tissue}_{sample}_{lane}_join.deets.sized.clean.txt"),
#    params:
#        conda_env = config['conda_envs']['coding'],
#    threads: 1
#    resources:
#        time   = 60,
#        mem_mb = 24000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            zcat {input.sized} \
#                | seqkit rmdup -s -i \
#                -o {output.dedup} \
#                -d {output.dupes} \
#                -D {output.deets}
#        '''

# USE MIREC CORRECTED OR ONLY ADAPTER REMOVED AND JOINED???
# USING THE LATER
rule convert_fastq_to_fasta:
    input:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq.gz'),
       #S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.clean.fastq.gz")
    output:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.pre_filter.clean.fasta'),
       #S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.pre_filter.clean.fasta"),
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
       #S3.remote(f"{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/noncode_filt/{config['release'][0]}/{{tissue}}_{{sample}}_join.sized.clean.mm01_ncref_rfam.unaligned.fastq"),
    output:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.post_filter.clean.fasta'),
       #S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.post_filter.clean.fasta"),
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

#rule convert_fastq_to_fasta_ge30:
#    input:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_{lane}_join.clean.fastq.gz"),
#    output:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_{lane}_join.ge30.sized.clean.fasta"),
#    params:
#        conda_env = config['conda_envs']['quant']
#    threads: 1
#    resources:
#        time   = 10,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            seqtk seq -A -L 30 {input} > {output}
#            # remove whitespace from headers
#            sed 's, ,-,g' -i {output}
#        '''
#rule convert_fastq_to_fasta_gt30:
#    input:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.clean.fastq.gz"),
#    output:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.gt30.sized.clean.fasta"),
#    params:
#        conda_env = config['conda_envs']['quant']
#    threads: 1
#    resources:
#        time   = 10,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            # drop seqs with length shorter than 31 (aka keep > 30)
#            seqtk seq -A -L 31 {input} > {output}
#            # remove whitespace from headers
#            sed 's, ,-,g' -i {output}
#        '''

# SKIP COLLAPSING STUFF - NO VALUE? DETERMINED IN small_hoof.smk ALL RULE
#rule collapse_sample_fasta:
#    input:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.sized.clean.fasta"),
#    output:
#        coll_fa = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.collapsed.sized.clean.fasta"),
#        all_fa  = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.all.sized.clean.fasta"),
#    params:
#        conda_env = config['conda_envs']['small']
#    threads: 1
#    resources:
#        time   = 30,
#        mem_mb = 12000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            mapper.pl \
#                {input} \
#                -c \
#                -m \
#                -s {output.coll_fa}
#            
#            # duplicate the input with modified name for ease of input to align
#            cp {input} {output.all_fa}
#
#        '''

# bowtie2 -x genomeindex -U out.fq|samtools view -bS - > out.bam
#samtools view -b {input.aligned} \ # SAM TO BAM
#    | samtools sort -O BAM \
#    | tee {output.sort_bam} \
#    | samtools index - {output.sort_bai}
rule small_align_ref:
    input:
        bt_idx   = rules.bowtie_index_ref.output.idx_dir,
        fna_idx  = S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna',keep_local=True),
        small_fa = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.{seq_set}.clean.fasta'),
       #small_fa = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.{seq_set}.clean.fasta"),
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
# ORIGINAL STRATEGY - PERHAPS TOO STRONG FOR UNIQUE...
#bowtie \
#    -f \
#    -n {wildcards.mm} \
#    -e 70 \
#    -l 18 \
#    -m 1 \
#    --norc \
#    -p {threads} \
#    {input.fna_idx} \
#    {input.small_fa} \
#    -S {output.sam} \
#    --un {output.unmap} \
#    2> {output.bt_log}

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
               #"{bucket}/private/small/align/bowtie/{release}/{u.tissue}/combine/mismatch_0{mm}/{seq_set}/{u.tissue}.mm0{mm}_{seq_set}.ref.featcts.txt",
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
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
                #tmp = os.path.basename(i).split(".")
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
        #oll_fa = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.collapsed.clean.fasta"),
        all_fa  = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.pre_filter.clean.fasta"),
    params:
        tissue_fa = lambda wildcards, input: " ".join(
            set([i for i in input.fastas if wildcards.tissue in i])
        ),
       #conda_env = config['conda_envs']['small']
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            # combine fastas
            cat {params.tissue_fa} > {output.all_fa}
        '''
# previously collapsed per tissue
#mapper.pl \
#    {output.all_fa} \
#    -c \
#    -m \
#    -s {output.coll_fa}

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
       #conda_env = config['conda_envs']['small']
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

# SHOULD THIS BE -s 2 SINCE I BELIEVE THIS IS REVERSE STRANDED (FR-FIRSTRAND)??
# BASED ON INFER_EXPERIMENT THE SMALL RNA IS FR-SECONDSTRAND OR FORWARD SO IT SHOULD BE 1
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

#rule gather_featcts_biotypes:
#    input:
#        ct_mat = S3.remote(
#            expand(
#                "{bucket}/private/small/align/bowtie/{release}/{u.tissue}/combine/mismatch_0{mm}/{seq_set}/{u.tissue}.mm0{mm}_{seq_set}.ref.featcts.txt",
#                u=units_small.itertuples(), 
#                bucket=config['bucket'], 
#                release=[
#                   #'GCF_002863925.1_EquCab3.0',
#                    'Equus_caballus.EquCab3.0.103'
#                ],
#                mm=list(range(3)),
#               #seq_set=['all','gt25']
#                seq_set=['pre_filter','post_filter']
#            )
#        )
#    output:
#        types = S3.remote("{bucket}/private/small/align/bowtie/{release}/all_tissues.featcts.csv", keep_local=True)
#    threads: 1
#    resources:
#        time   = 30,
#        mem_mb = 6000
#    run:
#        # extract biotype counts per mismatch/set/tissue and write to file
#        with open(output.types,"w") as out:
#            print("tissue,mismatch,seq_set,biotype,count",file=out)
#            # use list(set) due to the way the input expand works
#            for i in list(set(input.ct_mat)):
#                # get mismatch/set/tissue from file name
#                tmp = os.path.basename(i).split(".")
#                tissue = tmp[0]
#                mis,seq_set = tmp[1].rsplit("_",1)
#                with open(i,"r") as f:
#                    for line in f:
#                        if line.startswith(tuple(["#","Geneid"])):
#                            continue
#                        else:
#                            biotype = line.split("\t")[0]
#                            count = int(line.split('\t')[-1])
#                            # write to file
#                            print(tissue,mis,seq_set,biotype,count,sep=",",file=out)

########################################
## rerun above but with fasta input 
## >30 length due to peak around 31bp
########################################
#
#rule convert_fastq_to_fasta_ge30:
#    input:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_{lane}_join.clean.fastq.gz"),
#    output:
#        S3.remote("{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_{lane}_join.ge30.clean.fasta"),
#    params:
#        conda_env = config['conda_envs']['quant']
#    threads: 1
#    resources:
#        time   = 10,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            seqtk seq -A -L 30 {input} > {output}
#            # remove whitespace from headers
#            sed 's, ,-,g' -i {output}
#        '''
#
#rule combine_x_tissue_collapse_ge30:
#    input:
#        fastas = S3.remote(
#            expand(
#                "{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_{u.lane}_join.ge30.clean.fasta",
#                u=units_small.itertuples(), 
#                bucket=config['bucket'], 
#            )
#        ),
#    output:
#        all_fa  = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.all.ge30.clean.fasta"),
#    params:
#        tissue_fa = lambda wildcards, input: " ".join(
#            set([i for i in input.fastas if wildcards.tissue in i])
#        ),
#        conda_env = config['conda_envs']['small']
#    threads: 1
#    resources:
#        time   = 10,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            # combine fastas
#            cat {params.tissue_fa} > {output.all_fa}
#        '''
#
#rule small_align_tissue_ref_ge30:
#    input:
#        bt_idx    = rules.bowtie_index_ref.output.idx_dir,
#        fna_idx   = S3.remote("{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna",keep_local=True),
#        tissue_fa = S3.remote("{bucket}/private/fastq/trimmed/{tissue}/combine/small/{tissue}_join.all.ge30.clean.fasta"),
#    output:
#        sam    = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.sam"),
#        unmap  = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.unaligned.fasta"),
#        bt_log = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.bwt.log"),
#    params:
#        conda_env = config['conda_envs']['small']
#    threads: 12
#    resources:
#        time   = 240,
#        mem_mb = 12000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            bowtie \
#                -f \
#                -n {wildcards.mm} \
#                -e 70 \
#                -l 18 \
#                -a \
#                -m 10 \
#                --best \
#                --strata \
#                -p {threads} \
#                {input.fna_idx} \
#                {input.tissue_fa} \
#                -S {output.sam} \
#                --un {output.unmap} \
#                2> {output.bt_log}
#        '''
#
#rule feature_counts_tissue_ref_ge30:
#    input:
#        sam     = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.sam"),
#        ref_gtf = S3.remote('{bucket}/public/refgen/{release}/Equus_caballus.EquCab3.0.103_genomic.nice.gtf', keep_local=True),
#    output:
#        ct_mat = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.featcts.txt"),
#        fc_sum = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.featcts.txt.summary"),
#        fc_log = S3.remote("{bucket}/private/small/align/bowtie/{release}/{tissue}/combine/mismatch_0{mm}/{seq_set}/{tissue}.mm0{mm}_all_ge30.ref.featcts.log"),
#    params:
#        conda_env = config['conda_envs']['small'],
#        prefix    = lambda wildcards, output: os.path.splitext(output.ct_mat)[0],
#    threads: 2
#    resources:
#        time   = 30,
#        mem_mb = 12000
#    shell:
#        '''
#            set +eu
#            source activate {params.conda_env}
#
#            featureCounts \
#                -T {threads} \
#                -s 1 \
#                -M \
#                -t exon \
#                -g transcript_biotype \
#                -a {input.ref_gtf} \
#                -o {output.ct_mat} \
#                {input.sam} \
#                2> {output.fc_log}
#        '''
