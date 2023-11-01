#rule junc_modules_x_sample:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        sorted_bai    = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#        junc_sat_pdf  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junctionSaturation_plot.pdf"),
#        junc_sat_r    = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junctionSaturation_plot.r"),
#        junc_ann_xls  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junction.xls"),
#        junc_ann_r    = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junction_plot.r"),
#        junc_ann_splc = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.splice_events.pdf"),
#        junc_ann_junc = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.splice_junction.pdf"),
#        junc_ann_bed  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junction.bed"),
#        junc_ann_iner = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.junction.Interact.bed")
#    threads: 1
#    params:
#        base = lambda wildcards, output: str(Path.joinpath(
#            Path(output.junc_sat_pdf).parent,
#            Path(output.junc_sat_pdf).stem.split(".")[0]
#            )
#        ),
#        conda_env = config['conda_envs']['qc']
#    resources:
#        time   = 40,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#
#            samtools index {input.sorted_bam}
#
#            junction_saturation.py \
#                -i {input.sorted_bam} \
#                -r {input.bed} \
#                -o {params.base}
#
#            junction_annotation.py \
#                -i {input.sorted_bam} \
#                -r {input.bed} \
#                -o {params.base}
#        '''
#
#rule infer_experiment:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        infer_out = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.infer_experiment.txt"),
#    params:
#        conda_env = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 5,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            infer_experiment.py \
#                -r {input.bed} \
#                -i {input.sorted_bam} \
#                > {output.infer_out}
#        '''
#
#rule genebody_coverage:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        genbod_cov_txt = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.geneBodyCoverage.txt"),
#        genbod_cov_pdf = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.geneBodyCoverage.curves.pdf"),
#        genbod_cov_r   = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.geneBodyCoverage.r"),
#    params:
#        out_prefix = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}"),
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 360,
#        mem_mb = 24000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            geneBody_coverage.py \
#                -r {input.bed} \
#                -i {input.sorted_bam} \
#                -o {params.out_prefix}
#        '''
#
#rule inner_distance:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        inner_dis_txt  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.inner_distance.txt"),
#        inner_dis_freq = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.inner_distance_freq.txt"),
#        inner_dis_pdf  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.inner_distance_plot.pdf"),
#        inner_dis_r    = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.inner_distance_plot.r"),
#    params:
#        out_prefix = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}"),
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 20,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            inner_distance.py \
#                -r {input.bed} \
#                -i {input.sorted_bam} \
#                -o {params.out_prefix}
##        '''
#
#rule read_gc:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#    output:
#        read_gc_xls = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.GC.xls"),
#        read_gc_pdf = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.GC_plot.pdf"),
#        read_gc_r   = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.GC_plot.r"),
#    params:
#        out_prefix = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}"),
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 40,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            read_GC.py \
#                -i {input.sorted_bam} \
#                -o {params.out_prefix}
#        '''
#
#rule read_distribution:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        read_dis = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.read_dist.txt"),
#    params:
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 120,
#        mem_mb = 6000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            read_distribution.py \
#                -r {input.bed} \
#                -i {input.sorted_bam} \
#                > {output.read_dis}
#        '''
#
#rule read_duplication:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#        bed        = S3.remote("{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf.bed")
#    output:
#        pos_dup_xls  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.pos.DupRate.xls"),
#        seq_dup_xls  = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.seq.DupRate.xls"),
#        read_dup_pdf = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.DupRate_plot.pdf"),
#        read_dup_r   = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.DupRate_plot.r"),
#    params:
#        out_prefix = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}"),
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 480,
#        mem_mb = 24000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            read_duplication.py \
#                -i {input.sorted_bam} \
#                -o {params.out_prefix}
#        '''
#
#rule bam_stat:
#    input:
#        sorted_bam = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam"),
#        sorted_bai = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/{tissue}_{sample}_Aligned.sortedByCoord.out.bam.bai"),
#    output:
#        stat_out = S3.remote("{bucket}/private/align/star2pass_x_tissue/{release}/{tissue}/{sample}/rseqc/{tissue}_{sample}.bam_stat.txt"),
#    params:
#        conda_env  = config['conda_envs']['qc']
#    threads: 1
#    resources:
#        time   = 30,
#        mem_mb = 12000
#    shell:
#        '''
#            set +eu
#            
#            source activate {params.conda_env}
#            
#            bam_stat.py \
#                -i {input.sorted_bam} \
#                > {output.stat_out}
#        '''

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
       #S3.remote(
       #    sorted(expand(
       #        "{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}.trim.fastp.json",
       #        u=units_small.itertuples(),
       #        bucket=config['bucket'],
       #    ))
       #),
       #S3.remote(
       #    sorted(expand(
       #        "{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}.clean.fastp.json",
       #        u=units_small.itertuples(),
       #        bucket=config['bucket'],
       #    ))
       #),
       #S3.remote(
       #    sorted(expand(
       #        "{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}.post_join.fastp.json",
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #    ))
       #),
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
       #S3.remote(
       #    sorted(expand(
       #        '{bucket}/private/small/align/bowtie/filters/{release}/{u.tissue}/{u.sample}/mismatch_00/all/{u.tissue}_{u.sample}.mm00_all.ncref_rfam.bwt.log',
       #        u=units_small.itertuples(),
       #        bucket=config['bucket'],
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    ))
       #),
       #S3.remote(
       #    sorted(expand(
       #        '{bucket}/private/small/align/bowtie/{release}/{u.tissue}/{u.sample}/mismatch_0{mm}/{seq_set}/{u.tissue}_{u.sample}.mm0{mm}_{seq_set}.ref.featcts.txt.summary',
       #        u=units_small.itertuples(),
       #        bucket=config['bucket'],
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        mm=list(range(3)),
       #        seq_set=['pre_filter','post_filter']
       #    ))
       #),
       #S3.remote(
       #    sorted(expand(
       #        '{bucket}/private/small/align/bowtie/{release}/{u.tissue}/combine/mismatch_0{mm}/{seq_set}/{u.tissue}.mm0{mm}_{seq_set}.ref.featcts.txt.summary',
       #        u=units_small.itertuples(),
       #        bucket=config['bucket'],
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        mm=list(range(3)),
       #        seq_set=['pre_filter','post_filter']
       #    ))
       #),
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/align/bowtie/{release}/all_tissues.featcts.csv',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
        S3.remote('{bucket}/private/fastq/trimmed/all_small_join.sized.clean.lengths.csv'),
        S3.remote('{bucket}/private/fastq/trimmed/all_small_join.post_join.tsv'),
    output: 
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_report.html',keep_local=True),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc.log'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_sources.txt'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_data.json'),
        S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_fastp.txt'),
       #S3.remote('{bucket}/private/qc/{release}/small/multiqc_data/multiqc_featureCounts.txt'),
       #S3.remote('{bucket}/private/qc/{release}/multiqc_data/multiqc_general_stats.txt'),
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
