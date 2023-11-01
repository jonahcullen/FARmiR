
# NOTE - this is not generic yet...
rule blast_db_rnacentral:
    input:
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa')
    output:
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nsq'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.ntf'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nto'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nhr'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nin'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.njs'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.not'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.ndb'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            makeblastdb -in {input} -dbtype nucl
        '''

rule blast_novel_precurs_rnacent:
    input:
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nsq'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.ntf'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nto'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nhr'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.nin'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.njs'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.not'),
        S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa.ndb'),
        rnacent = S3.remote('{bucket}/public/rnacentral/v22/Equus_caballus_rnacentral.fa'),
        precurs = S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_precursor.fa'),
    output:
        S3.remote('{bucket}/private/small/novel/{release}/mirpro_novel_blastn.out'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            blastn -db {input.rnacent} \
                -query {input.precurs} \
                -evalue 1e-7 \
                -outfmt 6 \
                > {output}
        '''

rule filter_blast_rnacent:
    input:
        precurs = S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_precursor.fa'),
        blastn  = S3.remote('{bucket}/private/small/novel/{release}/mirpro_novel_blastn.out'),
    output:
        filter_fa  = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.fa'),
        filter_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.bed'),
    run:
        # remove sequences with hits from novel_precursor.fa
        df = pd.read_csv(input.blastn, sep='\t', header = None)
        seq_hits = list(set(df.iloc[:,0]))

        with open(input.precurs, 'r') as f_in, \
            open(output.filter_fa, 'w') as f_out:
            for seq in SeqIO.parse(f_in, 'fasta'):
                if seq.id not in seq_hits:
                    SeqIO.write(seq, f_out, 'fasta')

        # convert novel_precursor.nc_filter.fa to bed
        with open(output.filter_fa, 'r') as f_in, \
            open(output.filter_bed, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        line = line.strip().split('\t')
                        name = line[0][1:]
                        chrom,strand,start,end = line[-1].split(':')
                        print(
                            chrom, max(int(start)-1,0), # bed 0-based index 
                            end, name, '.', strand, 
                            file=f_out, sep='\t'
                        )

# merged intervals with a maximum overlap of 17 (ie <17 will not be merged or 
# treated separately) along the same chrom and strand
rule merge_precurs:
    input:
        ref_dict   = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.dict'),
        filter_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.bed'),
    output:
        genome    = S3.remote('{bucket}/private/small/novel/{release}/genome.txt'),
        sort_bed  = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.bed'),
        temp_bed  = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.TMP.bed'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            # genome file for sorting
            tail -n +2 {input.ref_dict} | cut -f 2,3 > {output.genome}
            sed -i 's/SN:\|LN://g' {output.genome}

            bedtools sort \
                -i {input.filter_bed} \
                -g {output.genome} \
                > {output.sort_bed}

            bedtools merge -s \
                -d -17 \
                -c 4,6 \
                -o distinct \
                -delim "|" \
                -i {output.sort_bed} \
                > {output.temp_bed}
        '''

rule reformat_merged_bed:
    input:
        temp_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.TMP.bed'),
    output:
        merge_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.merged.bed'),
    run:
        # the merged output needs to be reformatted
        with open(input.temp_bed, 'r') as f_in, \
            open(output.merge_bed, 'w') as f_out:
            for line in f_in:
                line = line.strip().split()
                start = int(line[1])+1
                mod_name = f'{line[-2]}&{line[0]}|{line[-1]}|{start}|{line[2]}'
                print(
                    line[0], line[1], line[2], 
                    mod_name, '.', line[-1], 
                    sep='\t', file=f_out
                )

rule blast_db_ref:
    input:
        S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna'),
    output:
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nsq'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.ntf'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nto'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nhr'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nin'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.njs'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.not'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.ndb'),
        blast_ref = S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            # modify fasta to remove white space after first ' '
            sed 's/\s.*$//' {input} > {output.blast_ref}

            makeblastdb -in {output.blast_ref} -dbtype nucl
        '''

rule blast_merged_ref:
    input:
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nsq'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.ntf'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nto'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nhr'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.nin'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.njs'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.not'),
        S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna.ndb'),
        ref_fa    = S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna'),
        merge_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.merged.bed'),
    output:
        merge_fa  = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.merged.fa'),
        merge_fai = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.nc_filter.sorted.merged.fa.fai'),
        blastn    = S3.remote('{bucket}/private/small/novel/{release}/filtered_novel_blastn.out'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            # get fasta 
            bedtools getfasta \
                -nameOnly \
                -fi {input.ref_fa} \
                -bed {input.merge_bed} \
                -fo {output.merge_fa}

            # index and blast against reference
            samtools faidx {output.merge_fa}

            blastn -db {input.ref_fa} \
                -query {output.merge_fa} \
                -evalue 1e-7 \
                -outfmt 6 \
                > {output.blastn}
        '''

# count the number of hits for each one as "interspersed repeats" to be removed
rule get_novel_candidates:
    input:
        blastn = S3.remote('{bucket}/private/small/novel/{release}/filtered_novel_blastn.out'),
    output:
        candidates = S3.remote('{bucket}/private/small/novel/{release}/novel_candidates.list'),
        many_hits  = S3.remote('{bucket}/private/small/novel/{release}/novel_candidates.many_hits.list'),
    run:
        all_hpins = []
        d = defaultdict(int)
        with open(output.candidates, 'w') as f_out, \
            open(output.many_hits, 'w') as many_out:
            with open(input.blastn, 'r') as f_in:
                for line in f_in:
                    line = line.strip().split('\t')
                    d[line[0]] += 1
                    if line[0] not in all_hpins:
                        all_hpins.append(line[0])
            # write to outs
            for k,v in d.items():
                if v <= 15:
                    print(k, file=f_out)
                else:
                    print(k, file=many_out)

# for each candidate from novel_candidates.list, if there is a '|' that means 
# the precursors were merged. to determine which one to keep: 1. keep the one 
# with the highest total expression by summing the counts. 2. if equal counts, 
# randomly pick
rule filter_by_norm_exp:
    input:
        precur_cts = S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/result_novel_precursor.csv'),
        candidates = S3.remote('{bucket}/private/small/novel/{release}/novel_candidates.list'),
        inc_filt   = S3.remote('{bucket}/private/small/quant/preproc_all_filters.inc_samples.tsv'),
        precurs    = S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_precursor.fa'),
        matures    = S3.remote('{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_mature.fa'),
    output:
        exp_precurs = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.fa'),
        exp_matures = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.exp_cands.fa'),
    run:
        # precursor count matrix
        pre_cts = pd.read_csv(
            input.precur_cts, 
            skiprows=2, header=0, index_col='id'
        )

        # modify sample names
        pre_cts.columns = pre_cts.columns.str.replace(
            '.mirna.unaligned_mm00_novel count', ''
        )

        keep_cands = []
        with open(input.candidates, 'r') as f_in:
            for line in f_in:
                line = line.strip().split('&')[0]
                if '|' in line:
                    keep_cands.append(exp_locus(pre_cts, line))
                else:
                    keep_cands.append(line)

        # convert to rpm using the final output from filters.smk for post process cts
        inc_samp = pd.read_csv(
            input.inc_filt, 
            sep='\t', usecols=['tissue_sample', 'post_counts']
        ).set_index('tissue_sample')['post_counts'].to_dict()

        samp_sums = pd.Series(inc_samp, index=pre_cts.columns)

        # divide each count by the column sum and multiply by one million
        rpm_precurs = pre_cts.div(samp_sums) * 1e6
        # select keep_cands rows
        rpm_keep_cands = rpm_precurs[rpm_precurs.index.isin(keep_cands)]
        # select precurs with rpm >= 20 in at least two samples
        exp_precurs = rpm_keep_cands[(rpm_keep_cands >= 20).sum(axis=1) >= 2]
        # write to list
        keep_exp_precurs = exp_precurs.index.tolist()

        # use keep_exp_precurs to filter both novel_precursor.fa and novel_mature.fa
        with open(input.precurs, 'r') as f_in, \
            open(output.exp_precurs, 'w') as f_out:
            for seq in SeqIO.parse(f_in, 'fasta'):
                if seq.id in keep_exp_precurs:
                    SeqIO.write(seq, f_out, 'fasta')

        with open(input.matures, 'r') as f_in, \
            open(output.exp_matures, 'w') as f_out:
                for seq in SeqIO.parse(f_in, 'fasta'):
                    mod_seq = seq.id.lower()
                    if seq.id.count('-') > 3:
                        mod_seq = seq.id.lower().rsplit('-',1)[0]
                    # check if in keep list
                    if mod_seq in keep_exp_precurs:
                        SeqIO.write(seq, f_out, 'fasta')

# create a bed file from novel_precursor.exp_cands.fa and check for overlaps 
# with eca3 mirna gff generated as part of make_nice.smk plus the manual 
# addition of ucsc lift off (eca3_primary.ens_manual.sorted.bed)
rule get_exp_precur_bed:
    input:
        exp_fa = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.fa'),
    output:
        exp_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.bed'),
    run:
        with open(input.exp_fa, 'r') as f_in, \
            open(output.exp_bed, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    line = line.strip().split('\t')
                    name = line[0][1:]
                    chrom,strand,start,end = line[-1].split(':')
                    mod_name = f'{name}&{chrom}|{strand}|{start}|{end}'
                    print(chrom, max(int(start)-1,0), end, mod_name, '.', strand, file=f_out, sep='\t')

rule overlap_exp_precur_and_ref:
    input:
        ref_mir = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3_primary.ens_manual.sorted.bed'),
        exp_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.bed'),
    output:
        exp_sort      = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.sorted.bed'),
        exp_ref_inter = S3.remote('{bucket}/private/small/novel/{release}/eca3_novel_intersect.bed'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            # sort bed and intersect
            sort -k 1,1 -k2,2n {input.exp_bed} > {output.exp_sort}

            bedtools intersect \
                -a {output.exp_sort} \
                -b {input.ref_mir} \
                -s -wo \
                > {output.exp_ref_inter}
        '''
        
# find novel ids with overlap greater than 16 (>=17)
rule remove_precurs_overlap_ref:
    input: 
        exp_sort      = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.exp_cands.sorted.bed'),
        exp_ref_inter = S3.remote('{bucket}/private/small/novel/{release}/eca3_novel_intersect.bed'),
        exp_matures   = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.exp_cands.fa'),
    output:
        precur_temp = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.tmp.bed'),
        final_mats   = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.final.fa'),
    run:
        inter_remove = []
        with open(input.exp_ref_inter, 'r') as f_in:
            for line in f_in:
                line = line.strip().split('\t')
                if int(line[-1]) > 16:
                    inter_remove.append(line[3])

        # generate pre-FINAL novel precursor bed and mature fa
        with open(input.exp_sort, 'r') as f_in, \
            open(output.precur_temp, 'w') as f_out:
            for line in f_in:
                line = line.strip().split()
                if line[3] not in inter_remove:
                    print('\t'.join(line), file=f_out)
            
        with open(input.exp_matures, 'r') as f_in, \
            open(output.final_mats, 'w') as f_out:
            for seq in SeqIO.parse(f_in, 'fasta'):
                mod_seq = seq.id.lower()
                if seq.id.count('-') > 3:
                    mod_seq = seq.id.lower().rsplit('-',1)[0]
                # check if in keep list
                if mod_seq not in [i.split('&')[0] for i in inter_remove]:
                    dna_seq = seq.seq.back_transcribe()
                    seq.seq = dna_seq
                    SeqIO.write(seq, f_out, 'fasta')

# get final precursor fasta
rule final_precur_fa:
    input:
        ref_fa      = S3.remote('{bucket}/public/refgen/{release}/BLAST_DB/{release}_genomic.nice.fna'),
        precur_temp = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.tmp.bed'),
    output:
        precur_fa = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.fa'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            bedtools getfasta \
                -nameOnly \
                -s \
                -fi {input.ref_fa} \
                -bed {input.precur_temp} \
                -fo {output.precur_fa}
        '''

rule final_precur_mature_beds:
    input:
        precur_fa   = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.fa'),
        precur_temp = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.tmp.bed'),
        final_mats  = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.final.fa'),
    output:
        final_mats_bed   = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.final.bed'),
        final_precur_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.bed'),
    run:
        # get into d and strip strand info from name [eg ...(-)]
        d_tmp = SeqIO.to_dict(SeqIO.parse(input.precur_fa, 'fasta'))
        d_prefa = {k[:-3]:v for k,v in d_tmp.items()} 

        # dictionary of precursors based on final fasta
        d_npins = defaultdict(dict)
        with open(input.precur_temp, 'r') as f_in:
            for line in f_in:
                line = line.strip().split()
                seq = str(d_prefa[line[3]].seq)
                hpin = line[3].split('&')[0]
                for i in range(0, 5):
                    mat = f'{hpin.replace("mir", "miR")}-{str(i)}'
                    if i == 4:
                        mat = hpin.replace('mir', 'miR')
                    d_npins[mat]['name'] = d_prefa[line[3]].id # start in name is 1-based
                    d_npins[mat]['chrom'] = line[0]
                    d_npins[mat]['start'] = line[1] # bed 0-based
                    d_npins[mat]['end'] = line[2]
                    d_npins[mat]['strand'] = line[-1]
                    d_npins[mat]['seq'] = seq
        
        # the final mature mirnas
        d_mat_fa = SeqIO.to_dict(SeqIO.parse(input.final_mats, 'fasta'))

        d_nmats = defaultdict(dict)
        for k,v in d_mat_fa.items():
            if k in d_npins:
                parent = k.replace('miR', 'mir')
                d_nmats[k]['parent'] = parent
                if k.count('-') > 3:
                    d_nmats[k]['parent'] = parent[:-2]
                d_nmats[k]['chrom'] = d_npins[k]['chrom']
                d_nmats[k]['strand'] = d_npins[k]['strand']
                if '+' in d_npins[k]['strand']:
                    start = int(d_npins[k]['start']) + d_npins[k]['seq'].find(str(v.seq))
                    end = start + len(str(v.seq))
                else: # negative strand
                    end = int(d_npins[k]['end']) - d_npins[k]['seq'].find(str(v.seq))
                    start = end - len(mat)
                # add to dict
                d_nmats[k]['mirna_start'] = start
                d_nmats[k]['mirna_end'] = end

            # write final mature to bed
            with open(output.final_mats_bed, 'w') as f_out:
                source = 'mirpro'
                type = 'miRNA'
                for k,v in d_nmats.items():
                    parent = f'MINO00{v["parent"].split("-")[-1]}'
                    alias = f'MIMATNO00{v["parent"].split("-")[-1]}'
                    if k.count('-') > 3:
                        mat_id = f'{alias}_{int(k[-1])+1}'
                    else:
                        mat_id = alias
                    name = f'ID={mat_id};Alias={alias};Name={k};Parent={parent}'
                    print(
                        v['chrom'], v['mirna_start'], v['mirna_end'],
                        '.', '.', v['strand'],
                        source, type, '.',
                        name,
                        sep='\t', file=f_out
                    )
        
            # write final haripin to bed
            with open(input.precur_temp, 'r') as f_in, \
                open(output.final_precur_bed, 'w') as f_out:
                source = 'mirpro'
                type = 'miRNA_primary_transcript'
                for line in f_in:
                    chrom,start,end,tmp,_,strand = line.strip().split()
                    mir_name = tmp.split('&')[0]
                    #if 'eca-novel-mir-7147' in mir_name:
                    #       print(chrom,start,end,tmp,strand)
                    alias = f'MINO00{mir_name.split("-")[-1]}'
                    name = f'ID={alias};Alias={alias};Name={mir_name}'
                    print(
                        chrom, start, end,
                        '.', '.', strand,
                        source, type, '.',
                        name,
                        sep='\t', file=f_out
                    )

rule combine_final_beds:
    input:
        final_mats_bed   = S3.remote('{bucket}/private/small/novel/{release}/novel_mature.final.bed'),
        final_precur_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_precursor.final.bed'),
    output:
        final_all_bed = S3.remote('{bucket}/private/small/novel/{release}/novel_mirna.sorted.bed', keep_local=True),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            cat {input.final_mats_bed} {input.final_precur_bed} \
                | sort -k 1,1 -k2,2n \
                > {output.final_all_bed}
        '''

