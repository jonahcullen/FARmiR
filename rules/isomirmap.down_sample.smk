
# get dataframe of minimum post count sample per tissue
rule minimum_tissue_count:
    input:
        inc_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv'),
        all_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.all_samples.tsv')
    output:
        tiss_ct = S3.remote('{bucket}/private/small/quant/minimum_tissue_count.inc_samples.tsv')
    run:
        df = pd.read_csv(input.inc_filt, sep='\t', usecols=['tissue','post_counts']) 
        # get minimum post count sample per tissue
        df_min = df.loc[df.groupby('tissue').post_counts.idxmin()]
        # get tissues that were already removed and not in included samples
        # (e.g. uterine fluid), minimum post count per tissue
        # note this is required as the downsampling rule is applied to all
        # samples regardless of filter status
        df_all = pd.read_csv(input.all_filt, sep='\t', usecols=['tissue','post_counts'])
        filt_df = df_all[~df_all['tissue'].isin(df_min['tissue'].tolist())]
        filt_min = filt_df.loc[filt_df.groupby('tissue').post_counts.idxmin()]
        # add to included minimum
        df_min = pd.concat([df_min, filt_min], ignore_index=True)     
        # if minimum post count less 1m, down sample to 1m in others
        df_min['sample_to'] = np.where(
            df_min['post_counts'] < 1000000, 
            1000000, 
            df_min['post_counts']
        )  
        # and write to file
        df_min.to_csv(output.tiss_ct, sep='\t', index=False)

rule downsample_sized_fastqs:
    input:
        sized   = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna.unaligned_mm00.fastq'),
        tiss_ct = S3.remote('{bucket}/private/small/quant/minimum_tissue_count.inc_samples.tsv')
    output:
        ds_sized = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna_sampled.unaligned_mm00.fastq'),
    threads: 4
    params:
       #seed      = lambda wildcards: wildcards.sample_n,
        conda_env = config['conda_envs']['quant'],
        seed_size = lambda wildcards, input: \
            pd.read_csv(input.tiss_ct,sep='\t')['sample_to'] \
                .loc[pd.read_csv(input.tiss_ct,sep='\t')['tissue'] == input.sized.split('/')[4]] \
                .item() 
    resources:
        time   = 60,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            seqtk sample -s{params.seed_size} {input.sized} {params.seed_size} > {output.ds_sized}
        '''

# NOTE - running all samples and will filter out those to be excluded
# during an aggregate rule
rule isomirmap_downsample_quant:
    input:
        ds_sized   = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{sample}/small/noncode_filt/{release}/{tissue}_{sample}.mirna_sampled.unaligned_mm00.fastq'),
        tables_cfg = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/tables.cfg'),
        space_fa   = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
        lookups    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_lookups.txt'),
        repeats    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_repeats.txt'),
        coords     = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_coords.txt'),
        snp_bundle = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_snps.txt'),
        inc_filt   = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv')
    output:
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-exclusive-isomiRs.expression.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-exclusive-isomiRs.expression.html'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-exclusive-isomiRs.countsmeta.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-ambiguous-isomiRs.expression.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-ambiguous-isomiRs.expression.html'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-ambiguous-isomiRs.countsmeta.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-snps-isomiRs.expression.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-snps-isomiRs.expression.html'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-snps-isomiRs.countsmeta.txt'),
        S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-isomiRs.miRCarta.gff3'),
        mirbase = S3.remote('{bucket}/private/small/quant/downsample/{release}/{tissue}/{sample}/{tissue}_{sample}-IsoMiRmap_v5-isomiRs.miRBase.gff3'),
    params:
        conda_env = config['conda_envs']['novel_small'],
        table_dir = lambda wildcards, input: os.path.dirname(input.tables_cfg),
        prefix    = lambda wildcards, output: os.path.join(
            os.path.dirname(output.mirbase), 
            f'{wildcards.tissue}_{wildcards.sample}'
        )
    threads: 4
    resources:
        time   = 60,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            python /opt/hoof/src/isoMiRmap/IsoMiRmap.py \
                --m {params.table_dir} \
                --p {params.prefix} \
                {input.ds_sized}
        '''

rule isomirmap_downsample_combine:
    input:
        mirbase = S3.remote(
            sorted(expand(
                '{bucket}/private/small/quant/downsample/{release}/{u.tissue}/{u.sample}/{u.tissue}_{u.sample}-IsoMiRmap_v5-isomiRs.miRBase.gff3',
                u=units_small.itertuples(), 
                bucket=config['bucket'], 
                release=[
                   #'GCF_002863925.1_EquCab3.0',
                    'Equus_caballus.EquCab3.0.103'
                ],
            ))
        ),
        inc_filt = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv')
    output:
        ct_mat  = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_cts.ds_isomirmap.tsv'),
        rpm_mat = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_rpm.ds_isomirmap.tsv'),
    params:
        filt_gffs = lambda wildcards, input:
            [
                i for i in input.mirbase 
                if os.path.basename(i).split('-IsoMiRmap')[0] in pd.read_csv(input.inc_filt, sep='\t')['tissue_sample'].tolist()
            ]
    threads: 4
    resources:
        time   = 360,
        mem_mb = 32000
    run:
        # ignore the append/concat future warning - should be changed eventually
        import warnings
        warnings.simplefilter(action='ignore', category=FutureWarning)
        # combine parsed gffs per included sample
        keeps = ['Read','id','Name','Variant','Parent','Expression','Filter']
        gff = None
        
        # loop through each and extract keeps info plus normalized RPM expression
        for i in params.filt_gffs:
            tmp = pd.DataFrame()
            sample = os.path.basename(i).rsplit('-IsoMiRmap')[0] # -IsoMiRmap
            with open(i, 'r') as f_in:
                for line in f_in:
                    if line.startswith('#') or 'pre_miRNA' in line: # drop pre miRNA as not mature
                        continue
                    line = line.strip().split('\t')
                    attr = line[-1].split('; ')
                    d = {}
                    for i in attr:
                        k,v = i.split('=')
                        if any([j in k for j in keeps]):
                            if 'Filter' in k:
                                d[k] = v[:-1] # drop the semicolon
                            else:
                                d[k] = v
                        if 'normalized' in k:
                            for j in v.split(','):
                                k2,v2 = j.split(':')
                                if 'RPM_inputFile' in k2:
                                    d[k2] = v2
                    # add to df
                    tmp = tmp.append(d, ignore_index=True)
            # rename count/rpm columns
            tmp = tmp.rename(
                columns={
                    'Expression': f'{sample}.Expression',
                    'RPM_inputFile': f'{sample}.RPM_inputFile'
                }
            )
            # set index to everything but expression/rpm cols    
            tmp = tmp.set_index(tmp.columns.difference(
                [f'{sample}.Expression', f'{sample}.RPM_inputFile']).tolist()
            )
            # add to growing matrix
            if gff is None:
                gff = tmp
            else:
                gff = gff.join(tmp, how='outer')

        # fill missing with 0 to indicate sample did not have mature seq
        gff.fillna(0, inplace=True)
        # subset for count matrix
        gff_ct = gff.filter(like='Expression')
        gff_ct.columns = gff_ct.columns.str.replace('.Expression', '')
        gff_ct.to_csv(output.ct_mat, sep='\t')
        # subset for count matrix
        gff_rpm = gff.filter(like='RPM_inputFile')
        gff_rpm.columns = gff_rpm.columns.str.replace('.RPM_inputFile', '')
        gff_rpm.to_csv(output.rpm_mat, sep='\t')

rule isomirmap_downsample_filter:
    input:
        ct_mat  = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_cts.ds_isomirmap.tsv'),
        rpm_mat = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_rpm.ds_isomirmap.tsv'),
    output:
        rpm_exp_pass = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_rpm.ds_isomirmap_exp_filt.tsv'),
        rpm_exp_fail = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_rpm.ds_isomirmap_exp_filt_sample_fail.tsv'),
        ct_exp_pass  = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_cts.ds_isomirmap_exp_filt.tsv'),
    run:
        # read in rpm and filter to include only mature mirans (rows) with rpm
        # values > 20 in at least 2 samples AND at least >= 10% of non-zero
        # counts
        rpm_df = pd.read_csv(
            input.rpm_mat,
            index_col=['Filter', 'Name', 'Parent', 'Read', 'Variant', 'id'],
            sep='\t'
        )
        exp_rpm_df = rpm_df[(rpm_df > 20).sum(axis=1) >= 2]
        # drop duplicates
        dedup_df = exp_rpm_df.groupby(exp_rpm_df.index.get_level_values('id')).first()
        # mean proportion of zeros per sample and identify those >= 0.90
        mean_props = (dedup_df == 0).mean()
        fail_samples = mean_props[mean_props >= 0.90].index.tolist()
        # subset to dfs with fail and pass samples
        fail_df = exp_rpm_df[fail_samples]
        pass_exp_rpm_df = exp_rpm_df[exp_rpm_df.columns.difference(fail_samples)]

        # read in cts and filter to include the same samples and mirnas as rpm 
        # exp filt
        ct_df = pd.read_csv(
            input.ct_mat,
            index_col=['Filter', 'Name', 'Parent', 'Read', 'Variant', 'id'],
            sep='\t'
        )
        pass_ct_df = ct_df[ct_df.columns.difference(fail_samples)]
        pass_exp_ct_df = pass_ct_df[pass_ct_df.index.isin(pass_exp_rpm_df.index)] 
        #######################################################################
        # NEED TO EITHER CREATE A NEW RULE OR FIGURE OUT HOW TO ID SAMPLES THAT
        # NOW HAVE < 2 SAMPLES/TISSUE
        #######################################################################
        # write to out
        pass_exp_rpm_df.to_csv(output.rpm_exp_pass, sep='\t', na_rep='NA')
        fail_df.to_csv(output.rpm_exp_fail, sep='\t', na_rep='NA')
        pass_exp_ct_df.to_csv(output.ct_exp_pass, sep='\t', na_rep='NA')

# for each mature sequence (uid), get the numbers of different precursors, ref
# matures found, isomirs, isomirs with seed shifts (5p), and filter type
rule isomirmap_downsample_descript:
    input:
        rpm_exp_pass = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_rpm.ds_isomirmap_exp_filt.tsv'),
    output:
        descript = S3.remote('{bucket}/private/small/quant/downsample/{release}/merge/isomir_descript.ds.tsv'),
    run:
        d = defaultdict(dict)

        with open(input.rpm_exp_pass, 'r') as f_in:
            next(f_in)
            for line in f_in:
                fltr,name,parent,read,variant,uid = line.strip().split()[:6]
                if name not in d:
                    d[name]['isomirs'] = 0
                    d[name]['seed_shift'] = 0
                    d[name]['parent_loci'] = []
                    d[name]['ref'] = 0
                    d[name]['ref_uid'] = 'isomir_only'
                d[name]['filter'] = fltr
                # number of isomirs
                if 'NA' not in variant:
                    d[name]['isomirs'] += 1
                else: # variant is NA == ref mir
                    d[name]['ref'] += 1
                    d[name]['ref_uid'] = uid
                # check if seed shifted
                if '_5p' in variant:
                    d[name]['seed_shift'] += 1
                # parent loci
                d[name]['parent_loci'].append(parent) if parent not in d[name]['parent_loci'] else None

        with open(output.descript, 'w') as f_out:
            print(
                'name', 'ref_uid', 'precur_num', 'ref_num', 
                'isomirs', 'seed_shift', 'filter', 
                sep='\t', file=f_out
            )
            for k,v in d.items():
                print(
                    k, v['ref_uid'], len(v['parent_loci']), v['ref'], 
                    v['isomirs'], v['seed_shift'], v['filter'], 
                    sep='\t', file=f_out
                )

