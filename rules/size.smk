
rule small_fastq_pre_filtered_lengths:
    input:
       #fq_join = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/{tissue}_{sample}_join.clean.fastq.gz'),
        fq_join = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.clean.fastq.gz'),
    output:
        dist = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/{tissue}_{sample}_join.pre_filtered.clean.length.txt'),
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            zcat {input.fq_join} \
                | awk '{{if(NR%4==2) print length($1)}}' \
                | sort -n \
                | uniq -c \
                > {output.dist}
        '''

rule small_fastq_filtered_lengths:
    input:
       #unmap = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.clean.mm0{mm}_ncref_rfam.unaligned.fastq'),
        unmap = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.min_sized.corrected.clean.mm0{mm}_ncref_rfam.unaligned.fastq'),
    output:
        dist = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}_join.mm0{mm}_filtered.clean.length.txt'),
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            cat {input.unmap} \
                | awk '{{if(NR%4==2) print length($1)}}' \
                | sort -n \
                | uniq -c \
                > {output.dist}
        '''

rule gather_small_lengths:
    input:
        pre_filt = S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_join.pre_filtered.clean.length.txt',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
            )
        ),
        post_filt = S3.remote(
            expand(
                '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}_join.mm0{mm}_filtered.clean.length.txt',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
               #mm=list(range(3)),
                mm=['0'] # NOTE - only considering perfect matches for non-coding filtering
            )
        ),
    output:
        lengths = S3.remote('{bucket}/private/fastq/trimmed/all_small_join.sized.clean.lengths.csv')
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    run:
        # process each length file and combine
        dist = input.pre_filt + input.post_filt
        tpm = None
        for f in sorted(list(set(dist))):
            df = pd.read_table(
                f, 
                sep=' ',
                names=['count', 'length'],
                skipinitialspace=True
            )
            # get sample, tissue, and mismatch
            tmp = f.split('/')
            tissue,sample = tmp[4],tmp[5]
            df['sample'] = sample
            df['tissue'] = tissue
            # mismatch
            df['mismatch'] = 'pre_filtered'
            if 'mm0' in f:
                df['mismatch'] = os.path.basename(f).split('_')[-2].split('.')[-1]
            # combine
            if tpm is None:
                tpm = df
            else:
                tpm = tpm.append(df)

        # write to file
        tpm.to_csv(output.lengths,sep=',',index=False)

