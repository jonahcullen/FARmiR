
rule convert_fastq_tabular:
    input:
        sized = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}.mirna.unaligned_mm00.fastq'),
       #sized = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}_join.max_sized.clean.mm00_ncref_rfam.unaligned.fastq'),
    output:
        sized_tab = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}.mirna.unaligned_mm00.freq.tsv'),
    threads: 4
    resources:
        time   = 60,
        mem_mb = 60000
    run:
        from Bio import SeqIO
        from collections import Counter

        seqs = {}
        # extract the sequnece to seqs dictionary and count
        with open(input.sized, 'r') as in_f:
            for record in SeqIO.parse(in_f, 'fastq'):
                seq = str(record.seq)
                if seq in seqs:
                    seqs[seq] += 1
                else:
                    seqs[seq] = 1 
        
        ## print in file
        with open(output.sized_tab, 'w') as out_f:
            for seq, count in sorted(seqs.items(), key=lambda x: x[1], reverse=True):
                print(f'{seq}\t{str(count)}', file=out_f)

#rule convert_fastq_tabular:
#    input:
#        unmap = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}_join.clean.mm00_ncref_rfam.unaligned.fastq'),
#    output:
#        sized_tab = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}_join.clean.mm00_ncref_rfam.unaligned.freq.tsv'),
#    threads: 4
#    resources:
#        time   = 30,
#        mem_mb = 12000
#    run:
#        from collections import defaultdict
#
#        # modified from https://github.com/HCGB-IGTP/HCGB_python_functions
#        def process_fastq(lines):
#            ks = ['name','sequence','optional','quality']
#            return {k: v for k, v in zip(ks, lines)}
#
#        # generate frequency file that contains the count of each join sized
#        # clean fastq sequence
#        freq_fasta = defaultdict(int)
#        
#        # read fastq    
#        n = 4
#        with open(input.unmap,'r') as fh:
#            lines = []
#            for line in fh:
#                lines.append(line.rstrip())
#                if len(lines) == n:
#                    record = process_fastq(lines)
#                    lines = []
#                    ## add sequences & count
#                    freq_fasta[record['sequence']] += 1
#
#        ## print in file
#        with open(output.sized_tab,'w') as out:
#            for k in sorted(freq_fasta.keys()):
#                print(k,freq_fasta[k],sep='\t',file=out)

# NOTE - the default is -sub 1 but based on Schmauch et al 2021 maybe that is not the best idea,
# perhaps should run at 1, 2, and 3 or 1-7 for the length of the seed sequence? NOPE, does not work.
rule miralign_fastq_tabular:        
    input:
       #sized_tab = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}_join.max_sized.clean.mm00_ncref_rfam.unaligned.freq.tsv'),
        sized_tab = S3.remote(f'{{bucket}}/private/fastq/trimmed/{{tissue}}/{{sample}}/small/noncode_filt/{config["release"][0]}/{{tissue}}_{{sample}}.mirna.unaligned_mm00.freq.tsv'),
        eca_pre   = S3.remote('{bucket}/public/mirbase/v22/eca_db/hairpin.fa',keep_local=True),
        eca_mat   = S3.remote('{bucket}/public/mirbase/v22/eca_db/mature.fa',keep_local=True),
        eca_str   = S3.remote('{bucket}/public/mirbase/v22/eca_db/miRNA.str',keep_local=True),
        eca3_gff  = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3.ens.gff',keep_local=True),
    output:
        mirna   = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/{tissue}_{sample}.mirna',keep_local=True),
        opts    = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/{tissue}_{sample}.mirna.opt'),
        nomap   = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/{tissue}_{sample}.mirna.nomap'),
        logfile = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/{tissue}_{sample}.mirna.log'),
    params:
        conda_env = config['conda_envs']['small'],
        mirna_db  = lambda wildcards, input: os.path.dirname(input.eca_pre),
        prefix    = lambda wildcards, output: output.mirna.rsplit(".",1)[0]
    threads: 4
    resources:
        time   = 60,
        mem_mb = 32000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            miraligner -Xms750m -Xmx32g \
                -db {params.mirna_db} \
                -sub 1 \
                -add 3 \
                -trim 3 \
                -s eca \
                -i {input.sized_tab} \
                -o {params.prefix} \
                2> {output.logfile}   
        '''

rule mirgff3_miraligner:
    input:
        mirna           = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/{tissue}_{sample}.mirna'),
        eca_pre         = S3.remote('{bucket}/public/mirbase/v22/eca_db/hairpin.fa',keep_local=True),
        eca3_mirtop_gff = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3.ens_mirtop.gff',keep_local=True),
    output:
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/mirtop.gff'),
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/{tissue}_{sample}.gff'),
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/log/trace.log'),
    params:
        conda_env = config['conda_envs']['small'],
        outdir    = lambda wildcards, input: os.path.join(os.path.dirname(input.mirna),'mirtop')
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            mirtop gff \
                --format seqbuster \
                --sps eca \
                --hairpin {input.eca_pre} \
                --gtf {input.eca3_mirtop_gff} \
                -o {params.outdir} \
                {input.mirna}
        '''

rule mirgff3_stats_counts_export:
    input:
        sample_gff      = S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/mirtop.gff'),
        eca_pre         = S3.remote('{bucket}/public/mirbase/v22/eca_db/hairpin.fa',keep_local=True),
        eca3_mirtop_gff = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3.ens_mirtop.gff',keep_local=True),
    output:
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/stats/mirtop_stats.txt'),
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/counts/mirtop.tsv'),
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/export/mirtop_rawData.tsv'),
        S3.remote('{bucket}/private/small/quant/miraligner/{tissue}/{sample}/mirtop/export/log/trace.log'),
    params:
        conda_env     = config['conda_envs']['small'],
        outdir_stats  = lambda wildcards, input: os.path.join(os.path.dirname(input.sample_gff),'stats'),
        outdir_counts = lambda wildcards, input: os.path.join(os.path.dirname(input.sample_gff),'counts'),
        outdir_export = lambda wildcards, input: os.path.join(os.path.dirname(input.sample_gff),'export')
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e

            # stats
            mirtop stats -o {params.outdir_stats} {input.sample_gff}

            # counts
            mirtop counts \
                -o {params.outdir_counts} \
                --gff {input.sample_gff} \
                --hairpin {input.eca_pre} \
                --sps eca \
                --gtf {input.eca3_mirtop_gff}

            # export
            mirtop export \
                -o {params.outdir_export} \
                --hairpin {input.eca_pre} \
                --sps eca \
                --gtf {input.eca3_mirtop_gff} \
                --format isomir \
                {input.sample_gff}
        '''

#rule tissue_miraligner_matrix:
#    input:
#        mirtops = S3.remote(
#            expand(
#                '{bucket}/private/small/quant/miraligner/{u.tissue}/{u.sample}/mirtop/counts/mirtop.tsv',
#                u=units_small.itertuples(), 
#                bucket=config['bucket'],
#            )
#        ),
#        samp_exclude = S3.remote(f'{{bucket}}/private/fastq/trimmed/mirna_exp_samples/{config["release"][0]}/rnatype_excluded.tsv'),
#    output:
#        filt  = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_counts.csv',keep_local=True),
#        dupes = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_counts.duplicates.csv',keep_local=True),
#        seqs  = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_seqs.csv',keep_local=True),
#    params:
#        mirtops_sub = lambda wildcards, input: ",".join(
#            [
#                i for i in input.mirtops
#                if wildcards.tissue )
#            ]
#    ),
#    threads: 1
#    resources:
#        time   = 120,
#        mem_mb = 60000
#    run:
#        # get list of samples to exclude
#        df = pd.read_csv(input.samp_exclude,sep='\t',index_col='RNA_TYPE')
#        drop = df['sample'].unique().tolist()
#        
#        # get all mirtops into a dictionary
#       #d = {list(pd.read_csv(i,sep='\t').columns)[-1]:i for i in input.mirtops}
#        d = {}
#        for i in input.mirtops:
#            tmp = i.split('/')
#            sample = f'{tmp[-5]}_{tmp[-4]}'
#            if sample in drop:
#                continue
#            d[list(pd.read_csv(i,sep='\t').columns)[-1]] = i
#        # generate complete matrix with duplicates
#        all_isomirs,all_seqs = generate_isomir_matrix(d)
#        # extract duplicate uids
#        all_isomirs_filt,all_isomirs_dup = discard_UID_duplicated(all_isomirs)
#        # write to files
#        all_isomirs_filt.to_csv(output.filt,quoting=csv.QUOTE_NONNUMERIC)
#        all_isomirs_dup.to_csv(output.dupes,quoting=csv.QUOTE_NONNUMERIC)
#        all_seqs.to_csv(output.seqs,quoting=csv.QUOTE_NONNUMERIC)

# exclude samples with lt 10% mirna reads as identified following mirtrace
# with rule distinguish_exp_mirna_samples
rule all_tissues_miraligner_matrix:
    input:
        mirtops = S3.remote(
            expand(
                '{bucket}/private/small/quant/miraligner/{u.tissue}/{u.sample}/mirtop/counts/mirtop.tsv',
                u=units_small.itertuples(), 
                bucket=config['bucket'],
            )
        ),
        samp_exclude = S3.remote(f'{{bucket}}/private/fastq/trimmed/mirna_exp_samples/{config["release"][0]}/rnatype_excluded.tsv'),
    output:
        filt  = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_counts.DUPE_TEST.csv',keep_local=True),
        dupes = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_counts.duplicates.DUPE_TEST.csv',keep_local=True),
        seqs  = S3.remote('{bucket}/private/small/quant/miraligner/combine/isomir_seqs.DUPE_TEST.csv',keep_local=True),
    threads: 1
    resources:
        time   = 120,
        mem_mb = 60000
    run:
        # get list of samples to exclude
        df = pd.read_csv(input.samp_exclude,sep='\t',index_col='RNA_TYPE')
        drop = df['sample'].unique().tolist()
        
        # get all mirtops into a dictionary
       #d = {list(pd.read_csv(i,sep='\t').columns)[-1]:i for i in input.mirtops}
        d = {}
        for i in input.mirtops:
            tmp = i.split('/')
            sample = f'{tmp[-5]}_{tmp[-4]}'
            if sample in drop:
                continue
            d[list(pd.read_csv(i,sep='\t').columns)[-1]] = i
        # generate complete matrix with duplicates
        all_isomirs,all_seqs = generate_isomir_matrix(d)
        # extract duplicate uids
        all_isomirs_filt,all_isomirs_dup = discard_UID_duplicated(all_isomirs)
        # write to files
        all_isomirs_filt.to_csv(output.filt,quoting=csv.QUOTE_NONNUMERIC)
        all_isomirs_dup.to_csv(output.dupes,quoting=csv.QUOTE_NONNUMERIC)
        all_seqs.to_csv(output.seqs,quoting=csv.QUOTE_NONNUMERIC)

#rule isomir_merge_count:
#    input:
#        mirtops = S3.remote(
#            expand(
#                "{bucket}/private/small/quant/miraligner/{release}/{u.tissue}/{u.sample}/mirtop/counts/mirtop.tsv",
#                u=units_small.itertuples(), 
#                bucket=config['bucket'], 
#                release=[
#                   #'GCF_002863925.1_EquCab3.0',
#                    'Equus_caballus.EquCab3.0.103'
#                ],
#            )
#        )
#    output:
#        shuffle = temp("{bucket}/private/small/quant/miraligner/{release}/combine/shuffle/isomir_merge.shuffle_{perm}.csv"),
#    threads: 1
#    resources:
#        time   = 120,
#        mem_mb = 60000
#    run:
#        # get all mirtops into a dictionary
#        d = {list(pd.read_csv(i,sep='\t').columns)[-1]:i for i in input.mirtops}
#        # generate complete matrix with duplicates
#        filt_cts,all_cts = generate_isomir_matrix(d,True,2+int(wildcards.perm))
#        # write to output
#        with open(output.shuffle,'w') as out:
#            print("step,perm,filtered,all",file=out)
#            for i,(j,k) in enumerate(zip(filt_cts,all_cts)):
#                print(
#                    f"step_{str(i+1).zfill(4)}",
#                    f"perm_{str(wildcards.perm).zfill(4)}",
#                    j,k,
#                    sep=",",
#                    file=out
#                )
#
#
#rule isomir_all_perms:
#    input:
#        mirtops = expand(
#            "{bucket}/private/small/quant/miraligner/{release}/combine/shuffle/isomir_merge.shuffle_{perm}.csv",
#            bucket=config['bucket'], 
#            release=[
#               #'GCF_002863925.1_EquCab3.0',
#                'Equus_caballus.EquCab3.0.103'
#            ],
#            perm=[str(i).zfill(4) for i in range(1,101)]
#        )
#    output:
#        all_perms = S3.remote("{bucket}/private/small/quant/miraligner/{release}/combine/shuffle/shuffle_all.csv",keep_local=True),
#    threads: 1
#    resources:
#        time   = 120,
#        mem_mb = 60000
#    run:
#        # combine perms into a list and concat all
#        dfs = []
#        for f in list(input.mirtops):
#            df = pd.read_csv(f)
#            dfs.append(df)
#        
#        # write each resample round tpm to file
#        df = pd.concat(dfs)
#        df.to_csv(output.all_perms,index=False)
#
