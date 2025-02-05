
rule fastp_small_pre_trim:
    input:
        unpack(get_small_fastqs)
    output:
        html = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.pre_trim.fastp.html'),
        json = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.pre_trim.fastp.json')
    params:
        conda_env = config['conda_envs']['trim'],
        title     = '{tissue}_{sample}_with_adapters',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            if [ {input.r1} == {input.r2} ]; then
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as SINGLE"
                fastp \
                    -i {input.r1} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --disable_quality_filtering \
                    --disable_length_filtering \
                    --disable_adapter_trimming \
                    --disable_trim_poly_g \
                    --thread {threads}
            else
                echo "Processing {wildcards.sample}i ({wildcards.tissue}) as PAIRED"
                fastp \
                    -i {input.r1} -I {input.r2} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --disable_quality_filtering \
                    --disable_length_filtering \
                    --disable_adapter_trimming \
                    --disable_trim_poly_g \
                    --thread {threads}
            fi
        '''

# due to observed issues with fastp automatically finding adapters from se
# samples (could be older samples perhaps?), added a step to find adapters
# and the sequence of which can be passed directly to fastp
rule remove_adapters:
    input:
        unpack(get_small_fastqs)
    output:
        master_adapt = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/master_adapters.csv'),
        se_trimmed   = temp('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/good-mapping/{tissue}_{sample}_trimmed.fastq'),
        r1_trim      = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R1.trim.fastq.gz'),
        r2_trim      = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R2.trim.fastq.gz'),
        html         = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.trim.fastp.html'),
        json         = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.trim.fastp.json'),
    params:
        r1_base    = lambda wildcards, output: os.path.dirname(output.master_adapt),
        r1_fastq   = lambda wildcards, input: os.path.basename(str(input.r1)),
        trim_fastq = lambda wildcards, output: os.path.basename(output.se_trimmed),
        cut_env    = config['conda_envs']['cutadapt'],
        fastp_env  = config['conda_envs']['trim'],
        title      = '{tissue}_{sample}_remove_adapters',
    threads: 4
    resources:
        time   = 120,
        mem_mb = 60000
    shell:
        '''
            set +eu

            if [ {input.r1} == {input.r2} ]; then
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as SINGLE"
                source activate {params.cut_env}

                # due to the way adapt_find works, need to cd to each sample
                # directory prior to running otherwise master_adapters.csv will
                # be overwritten by each sample running in parallel
                cp {input.r1} {params.r1_base} && cd {params.r1_base}
                # from here run adapt_find, the adapter removed fastq will be
                # available in the good-mapping directory, and return
                adapt_find.py ILLUMINA --files {params.r1_fastq} --output_path ./ &&
                ADAPT_SEQ=$(tail -n +2 master_adapters.csv | cut -d ',' -f 4)
                if [[ "$ADAPT_SEQ" != "na" && -n "$ADAPT_SEQ" ]]; then
                    ADAPT=$(realpath good-mapping/*.fastq)
                    # rename the trimmed output
                    mv $ADAPT good-mapping/{params.trim_fastq}
                else
                    echo "No adapter found (adapter == $ADAPT_SEQ)"
                    mkdir -p ./good-mapping ./aux_files/{wildcards.sample}
                    echo "Running dnapi.py in iterative mode followed by cutadapt"
                    ADAPT_SEQ=$(python /opt/hoof/src/DNApi/dnapi.py -k 9:11:2 -r 1.2:1.4:0.1 {params.r1_fastq}) &&
                    cutadapt \
                        -q 20 \
                        -m 15 \
                        -M 50 \
                        -a "$ADAPT_SEQ" \
                        -o ./good-mapping/{params.trim_fastq} \
                        {params.r1_fastq} \
                        > ./aux_files/{wildcards.sample}/{wildcards.tissue}_{wildcards.sample}_cutadapt.txt
                    echo "Found and removed adapter $ADAPT_SEQ"
                    # for at least one sample (mus. sac. 11435), adapt_find failed
                    # leaving ADAPT_SEQ blank and failing to make master_adat
                    if [[ ! -e {output.master_adapt} ]]; then
                        echo "Adapt seq empty, creating empty master adapters"
                        touch master_adapters.csv
                    fi
                fi
                cd -
                
                conda deactivate
                source activate {params.fastp_env}
           
                fastp \
                    -i {output.se_trimmed} \
                    -o {output.r1_trim} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --disable_quality_filtering \
                    --disable_length_filtering \
                    --disable_adapter_trimming \
                    --disable_trim_poly_g \
                    --thread {threads}
                
                echo "Create empty R2 read"
                # touch for empty r2
                touch {output.r2_trim}
            else
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as PAIRED"
                source activate {params.fastp_env}

                fastp \
                    -i {input.r1} -I {input.r2} \
                    -o {output.r1_trim} -O {output.r2_trim} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --detect_adapter_for_pe \
                    --disable_quality_filtering \
                    --disable_length_filtering \
                    --disable_trim_poly_g \
                    --thread {threads}

                # touch for empty barcode and good-mapping files
                touch {output.master_adapt} {output.se_trimmed}
                echo "PE data, adapter detected automatically"
            fi
        '''

rule fastp_small_post_trim:
    input:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.pre_trim.fastp.json'),
        r1_trim = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R1.trim.fastq.gz'),
        r2_trim = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R2.trim.fastq.gz'),
    output:
        r1_clean = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R1.clean.fastq.gz'),
        r2_clean = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R2.clean.fastq.gz'),
        html     = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.clean.fastp.html'),
        json     = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.clean.fastp.json')
    params:
        title     = '{tissue}_{sample}_clean',
        conda_env = config['conda_envs']['trim']
    threads: 4
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            if [ -s {input.r2_trim} ]; then
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as PAIRED"
                fastp \
                    -i {input.r1_trim} -I {input.r2_trim} \
                    -o {output.r1_clean} -O {output.r2_clean} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --disable_length_filtering \
                    --disable_adapter_trimming \
                    --thread {threads}
            else
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as SINGLE"
                fastp \
                    -i {input.r1_trim} \
                    -o {output.r1_clean} \
                    --dont_overwrite \
                    --html {output.html} \
                    --json {output.json} \
                    --report_title {params.title} \
                    --disable_length_filtering \
                    --disable_adapter_trimming \
                    --thread {threads}
                
                # touch for empty r2
                touch {output.r2_clean}
            fi
        '''

rule fastq_join_small_clean:
    input:
        r1_clean = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R1.clean.fastq.gz'),
        r2_clean = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_R2.clean.fastq.gz'),
    output:
        fq_un1  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_un1.clean.fastq.gz'),
        fq_un2  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_un2.clean.fastq.gz'),
        fq_join = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq.gz'),
    params:
        conda_env = config['conda_envs']['cutadapt'],
        fq_join   = lambda wildcards, output: output.fq_un1.rsplit('un1',1)[0]+'%.clean.fastq.gz'
    threads: 1
    resources:
        time   = 120,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            if [ -s {input.r2_clean} ]; then
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as PAIRED"
                fastq-join \
                    -p 0 \
                    -m 6 \
                    {input.r1_clean} {input.r2_clean} \
                    -o {params.fq_join}
            else
                echo "Processing {wildcards.sample} ({wildcards.tissue}) as SINGLE"
                cp {input.r1_clean} {output.fq_join}
                touch {output.fq_un1} {output.fq_un2}
            fi
        '''

rule mirec_cleanup_fastq:
    input:
        fq_join = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq.gz'),
    output:
        fq_join   = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq'),
        corrected = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.corrected.clean.fastq.gz'),
    params:
        base      = lambda wildcards, input: os.path.dirname(input.fq_join),
        in_fastq  = lambda wildcards, output: os.path.basename(output.fq_join),
        out_fastq = lambda wildcards, output: os.path.basename(output.corrected).rsplit('.',1)[0]
    threads: 20
    resources:
        time   = 720,
        mem_mb = 160000
    shell:
        '''
            set -e

            # gunzip input
            gunzip {input.fq_join}

            # due to the way mirec works, need to process each sample
            # separately within a different directory
            mkdir -p {params.base}
            cd {params.base}

            /opt/hoof/src/miREC/miREC.sh \
                -f {params.in_fastq} \
                -s 8 \
                -e 15 \
                -t {threads} \
                -o {params.out_fastq} \
            
            gzip {params.out_fastq}
        '''

rule fastp_small_post_join:
    input:
        corrected = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.corrected.clean.fastq.gz'),
    output:
        sized = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.min_sized.corrected.clean.fastq.gz'),
        html  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.post_join.fastp.html'),
        json  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}.post_join.fastp.json')
    params:
        conda_env  = config['conda_envs']['trim'],
        title      = '{tissue}_{sample}_post_join',
        min_length = 17,
    threads: 4
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            fastp \
                -i {input.corrected} \
                -o {output.sized} \
                --length_required {params.min_length} \
                --dont_overwrite \
                --html {output.html} \
                --json {output.json} \
                --report_title {params.title} \
                --disable_quality_filtering \
                --disable_adapter_trimming \
                --disable_trim_poly_g \
                --thread {threads}
        '''

localrules: gather_post_join_logs
rule gather_post_join_logs:
    input:
        jsons = S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.post_join.fastp.json',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            )
        ),
    output:
        posts = S3.remote('{bucket}/private/fastq/trimmed/all_small_join.post_join.tsv')
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    run:
        # read post join jsons and extract pre- and post-filtered counts
        d = {}
        for i in input.jsons:
            res = re.match('(^[^0-9]+)_(.*)$', os.path.basename(i).rsplit('.', 3)[0])
            sample = res.group(2)
            tissue = res.group(1)
            key = f'{sample}_{tissue}'
            # load json
            with open(i, 'r') as json_f:
                    data=json_f.read()
            obj = json.loads(data)
            d[key] = [
                sample,
                tissue,
                obj['summary']['before_filtering']['total_reads'], 
                obj['summary']['after_filtering']['total_reads']
            ]
        # write d to file
        with open(output.posts, 'w') as out:
            print('sample\ttissue\tpre_filt\tpost_filt', file=out)
            for v in d.values():
                print(*v,sep='\t',file=out)

localrules: filter_mirna_nc_ref
rule filter_mirna_nc_ref:
    input:
       ref_ncrna = S3.remote('{bucket}/public/refgen/{release}/Equus_caballus.EquCab3.0.ncrna.fa.gz')
    output:
        ref_sans_mirna = S3.remote('{bucket}/public/refgen/{release}/{release}.ncrna_sans_mirna.fa')
    run:
        # filter out mirna sequences from the ref noncoding fasta
        with open(output.ref_sans_mirna,'w') as out, gzip.open(input.ref_ncrna,'rt') as fa:
            for rec in SeqIO.parse(fa,'fasta'):
                if 'transcript_biotype:miRNA' not in rec.description:
                    print(rec.format('fasta'),end='',file=out)

localrules: bowtie_index_ncrna_ref
rule bowtie_index_ncrna_ref:
    input:
        ref_sans_mirna  = S3.remote('{bucket}/public/refgen/{release}/{release}.ncrna_sans_mirna.fa'),
    output:
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.2.ebwt',keep_local=True),
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.3.ebwt',keep_local=True),
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.4.ebwt',keep_local=True),
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.rev.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa.rev.2.ebwt',keep_local=True),
        ref_fa      = S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa',keep_local=True),
        ref_idx_dir = directory('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF'),
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

            # build bowtie index for ref ncrna sans mirna
            cp {input.ref_sans_mirna} {output.ref_idx_dir}
            bowtie-build {output.ref_fa} {output.ref_fa}
        '''

# get fastq of reads not aligned against noncoding ref sans mirna
# to be further filtered against rfam sans mirna
rule small_align_ncref_sans_mirna:
    input:
        bt_idx   = rules.bowtie_index_ncrna_ref.output.ref_idx_dir,
        ref_fa   = S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/NCREF/{release}.ncrna_sans_mirna.fa',keep_local=True),
        sized_fq = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.min_sized.corrected.clean.fastq.gz'),
    output:
        sam    = S3.remote('{bucket}/private/small/align/bowtie/filters/{release}/{tissue}/{sample}/mismatch_0{mm}/all/{tissue}_{sample}.mm0{mm}_all.ncref.sam'),
        bt_log = S3.remote('{bucket}/private/small/align/bowtie/filters/{release}/{tissue}/{sample}/mismatch_0{mm}/all/{tissue}_{sample}.mm0{mm}_all.ncref.bwt.log'),
        unmap  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.corrected.clean.mm0{mm}_ncref.unaligned.fastq'),
    params:
        conda_env = config['conda_envs']['small']
    threads: 12
    resources:
        time   = 120,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -q \
                -n {wildcards.mm} \
                -m 1 \
                -e 70 \
                --norc \
                --no-unal \
                -p {threads} \
                {input.ref_fa} \
                {input.sized_fq} \
                -S {output.sam} \
                --un {output.unmap} \
                2> {output.bt_log}
        '''

localrules: bowtie_index_ncrna_rfam
rule bowtie_index_ncrna_rfam:
    input:
        rfam_sans_mirna = S3.remote('{bucket}/public/rfam/v14/rfam148_ncrna.sans_mirna.fa')
    output:
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.2.ebwt',keep_local=True),
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.3.ebwt',keep_local=True),
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.4.ebwt',keep_local=True),
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.rev.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa.rev.2.ebwt',keep_local=True),
        rfam_fa      = S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa',keep_local=True),
        rfam_idx_dir = directory('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM'),
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

            # replace U with T prior to indexing
            sed '/^[^>]/s/u/t/g' {input.rfam_sans_mirna} > {output.rfam_fa}

            # build bowtie index for rfam ncrna sans mirna
            bowtie-build {output.rfam_fa} {output.rfam_fa}
        '''

rule small_align_rfam_sans_mirna:
    input:
        bt_idx = rules.bowtie_index_ncrna_rfam.output.rfam_idx_dir,
        ref_fa = S3.remote('{bucket}/public/rfam/v14/BOWTIE_INDICES/RFAM/rfam148_ncrna.sans_mirna.dna.fa',keep_local=True),
        unmap  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.corrected.clean.mm0{mm}_ncref.unaligned.fastq'),
    output:
        sam    = S3.remote('{bucket}/private/small/align/bowtie/filters/{release}/{tissue}/{sample}/mismatch_0{mm}/all/{tissue}_{sample}.mm0{mm}_all.ncref_rfam.sam'),
        bt_log = S3.remote('{bucket}/private/small/align/bowtie/filters/{release}/{tissue}/{sample}/mismatch_0{mm}/all/{tissue}_{sample}.mm0{mm}_all.ncref_rfam.bwt.log'),
        unmap  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.corrected.clean.mm0{mm}_ncref_rfam.unaligned.fastq'),
    params:
        conda_env = config['conda_envs']['small']
    threads: 12
    resources:
        time   = 120,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -q \
                -n {wildcards.mm} \
                -m 1 \
                -e 70 \
                --norc \
                --no-unal \
                -p {threads} \
                {input.ref_fa} \
                {input.unmap} \
                -S {output.sam} \
                --un {output.unmap} \
                2> {output.bt_log}
        '''

rule fastp_size_select_unmap:
    input:
        unmap  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.corrected.clean.mm0{mm}_ncref_rfam.unaligned.fastq'),
    output:
        sized = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna.unaligned_mm0{mm}.fastq'),
        html  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mm0{mm}.fastp.html'),
        json  = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mm0{mm}.fastp.json')
    params:
        conda_env  = config['conda_envs']['trim'],
        title      = '{tissue}_{sample}_mm0{mm}_max_sized',
        max_length = 25
    threads: 4
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            fastp \
                -i {input.unmap} \
                -o {output.sized} \
                --length_limit {params.max_length} \
                --dont_overwrite \
                --html {output.html} \
                --json {output.json} \
                --report_title {params.title} \
                --disable_quality_filtering \
                --disable_adapter_trimming \
                --disable_trim_poly_g \
                --thread {threads}
        '''

localrules: gather_fastp_jsons
rule gather_fastp_jsons:
    input: 
        process_steps = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.{step}.fastp.json',
                u=units_small.itertuples(),
                bucket=config['bucket'],
                step=['pre_trim', 'trim', 'post_join']
            ))
        ),
        final_step = S3.remote(
            sorted(expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mm0{mm}.fastp.json',
                u=units_small.itertuples(),
                bucket=config['bucket'],
                release=[
                    'Equus_caballus.EquCab3.0.103'
                ],
                mm='0',
            ))
        )
    output:
        step_log = S3.remote('{bucket}/private/fastq/process_fastqs_counts.csv')
    run:
        def samp_tiss(s):
            return os.path.basename(s).split('.',1)[0]
        # json process order
        order = ['pre_trim', 'trim', 'post_join', 'mm00']
        def json_sorting(s):
            match = re.search(r'\.(pre_trim|trim|post_join|mm00)\.', s)
            if match:
                return order.index(match.group(1))
            return len(order)
        d = {}
        # combine and sort process and final steps
        all_steps = sorted(input.process_steps + input.final_step)
        for i in range(0, len(all_steps), 4):
            # ensure tissue_sample identical for each set of 4
            sample_set = all_steps[i:i+4]
            if not all(samp_tiss(i) == samp_tiss(sample_set[0]) for i in sample_set):
                print('ERROR: pre and post jsons not matched!')
                break
            tissue, sample = sample_set[0].split('/')[4:6]
            # specific to this dataset drop REP sample
            if 'REP' in sample:
                continue
            tissue_sample = f'{tissue}_{sample}'
            d[tissue_sample] = [0]*6
            for ind, val in enumerate(sorted(sample_set, key=json_sorting)):
                print(val)
                with open(val, 'r') as f_in:
                    step = json.load(f_in)
                    if ind < 2: # only need passed reads from first two steps (pre_trim, trim)
                        count = step['filtering_result']['passed_filter_reads']
                        d[tissue_sample][ind] = count
                    elif ind == 2:
                        pre_count = step['read1_before_filtering']['total_reads']
                        post_count = step['filtering_result']['passed_filter_reads']
                        d[tissue_sample][2] = pre_count
                        d[tissue_sample][3] = post_count
                    else:
                        pre_count = step['read1_before_filtering']['total_reads']
                        post_count = step['filtering_result']['passed_filter_reads']
                        d[tissue_sample][4] = pre_count
                        d[tissue_sample][5] = post_count

            print(d)
        # write to file
        with open(output.step_log, 'w') as f_out:
            print('tissue_sample,seq_depth,trimmed_reads,clean_min_length,post_join_mirec,pre_mirnas,max_length', file=f_out)
            for k,v in d.items():
                print(k, ','.join(map(str,v)), sep=',', file=f_out)

rule convert_corrected_to_fasta:
    input:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna.unaligned_mm00.fastq'),
    output:
        S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna.unaligned_mm00.fasta'),
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
        '''

