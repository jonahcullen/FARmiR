# generation of the species mapping bundle for isomirmap

import os
from textwrap import dedent
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# final output was saved to secondary manually at {bucket}/public/mirbase/v22/eca_db/eca3.ens_manual_novel.sorted.bed
# NOTE - very equine specific
localrules: refseq_ensembl_table
rule refseq_ensembl_table:
    input:
        ens = S3.remote('{bucket}/public/refgen/Equus_caballus.EquCab3.0.103/Equus_caballus.EquCab3.0.103_genomic.nice.dict'),
        ref = S3.remote('{bucket}/public/refgen/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.dict')
    output:
        ref_ens_table = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/refseq_ens_table.tsv')
    run:
        # convert refseq chrUn to ensembl's PJAA
        refseq = pd.read_csv(
            input.ref, 
            sep='\t', skiprows=1, usecols = [1,3],
            names = ['refseq_contig', 'm5']
        )

        ensembl = pd.read_csv(
            input.ens, 
            sep='\t', skiprows=1, usecols = [1,3],
            names = ['ensembl_contig', 'm5']
        )

        merged_df = pd.merge(refseq, ensembl, on='m5')
        merged_df['refseq_mod'] = merged_df['refseq_contig'].str.replace('SN:', 'chrUn_').str.replace('\.1', 'v1')
        merged_df['ensembl_mod'] = merged_df['ensembl_contig'].str[3:]
        merged_df.to_csv(output.ref_ens_table, index=False, sep='\t')

#####################
# mirna space
#####################
# add 6 bp up- and down-stream of hairpins
# NOTE - input.bed was generated manually above so input.novel is 
# included as an input here to ensure those rules are completed prior to
# generating the mapping bundle
localrules: hairpin_flanks
rule hairpin_flanks:
    input:
        novel = S3.remote('{bucket}/private/small/novel/{release}/novel_mirna.sorted.bed'),
        bed   = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3.ens_manual_novel.sorted.bed')
    output:
        flanks = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/eca3.ens_manual_novel.sorted_flanks.bed')
    run:
        with open(input.bed, 'r') as f_in, open(output.flanks, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                if 'miRNA_primary_transcript' in line:
                    line = line.strip().split('\t')
                    start = max(int(line[1]) - 6,0) # note bedtools is 0-based
                    end = int(line[2]) + 6 # and inclusive
                    strand = line[5]
                    chrom = line[0]
                    # get id and hairpin name
                    mi_id, hp_id = [
                        x.split('=')[1] 
                        for x in line[-1].split(';') 
                        if x.startswith('ID=') or x.startswith('Name=')
                    ]
                    name = f'{mi_id}|{hp_id}&WithFlank&{chrom}|{strand}|{start+1}|{end}'
                    print(chrom, start, end, name, '0', strand, sep='\t', file=f_out)

rule flank_fa:
    input:
        ref_fa = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna'),
        flanks = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/eca3.ens_manual_novel.sorted_flanks.bed')
    output:
        flanks_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/eca3.ens_manual_novel.sorted_flanks.fa'),
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            bedtools getfasta \
                -nameOnly \
                -s \
                -fi {input.ref_fa} \
                -bed {input.flanks} \
                -fo {output.flanks_fa}
        '''

# clean up flanks fa by replacing .1 with v1 in PJAA contigs (eca3 specific?)
localrules: mirnaspace_bundle
rule mirnaspace_bundle:
    input:
        flanks_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/eca3.ens_manual_novel.sorted_flanks.fa'),
    output:
        space_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
    run:
        with open(input.flanks_fa, 'r') as f_in, open(output.space_fa, 'w') as f_out:
            for line in f_in:
                line = line.strip()
                if line.startswith('>'):
                    line = line[:-3]
                    if 'PJAA' in line:
                        line = line.replace('.1','v1')
                print(line, file=f_out)

#####################
# lookup table
#####################
# for each hairpin, generate all possible substrings of length 17 to 25
localrules: hairpin_substrings
rule hairpin_substrings:
    input:
        bundle_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
    output:
        substrings_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_substrings.fa'),
    run:
        hpin_d = defaultdict(list)

        with open(input.bundle_fa, 'r') as f_in:
            for seq in SeqIO.parse(f_in, 'fasta'):
                for i in range(17,26):
                    hpin_subs = get_substrings(seq.seq, i)
                    for j in hpin_subs:
                        hpin_d[j].append(seq.id.split('|')[0])

        # write the substrings to file and align
        with open(output.substrings_fa, 'w') as f_out:
            for k,v in hpin_d.items():
                print('>'+'|'.join(v), k, sep='\n', file=f_out)


rule substrings_align_ref:
    input:
        bt_idx  = rules.bowtie_index_ref.output.idx_dir,
        fna_idx = S3.remote('{bucket}/public/refgen/{release}/BOWTIE_INDICES/GENOMIC/{release}_genomic.nice.fna'),
        subs_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_substrings.fa'),
    output:
        sam    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_ref.sam'),
       #unmap  = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_ref.unmap'),
        bt_log = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_ref.log'),
    params:
        conda_env = config['conda_envs']['novel_small']
    threads: 6
    resources:
        time   = 30,
        mem_mb = 60000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -f \
                -a \
                -n 0 \
                -p 6 \
                --sam-nosq \
                {input.fna_idx} \
                {input.subs_fa} \
                -S {output.sam} \
                2> {output.bt_log}
        '''

localrules: lookups_bundle
rule lookups_bundle:
    input:
        space_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
        hpin_sam = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_ref.sam'),
    output: 
        lookups = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_lookups.txt'),
    run:
        # dict of mi hairpin_ids
        hpin_locs = defaultdict(dict)

        with open(input.space_fa, 'r') as f_in:
            for seq in SeqIO.parse(f_in, 'fasta'):
                tmp = seq.id.split('|') 
                hpin,_,chrom = tmp[1].split('&')
                # add hpin info to dict
                hpin_locs[tmp[0]]['name'] = hpin
                hpin_locs[tmp[0]]['chrom'] = chrom
                hpin_locs[tmp[0]]['strand'] = tmp[2]
                hpin_locs[tmp[0]]['start'] = int(tmp[-2])
                hpin_locs[tmp[0]]['end'] = int(tmp[-1])

        # dict of yes/no exclusive
        lookups = defaultdict(str)

        with open(input.hpin_sam, 'r') as f_in:
            for line in f_in:
                if line.startswith('@'):
                    continue
                # get each alignment into a dict
                query = query_alignment(line.strip())
                seq = Seq(line.strip().split()[9])
                if query['flag'] == '16':
                    seq = seq.reverse_complement()
                seq = str(seq)
                # add sequence to targets with default 'Y' and only
                # subsequently to 'N' if the sequence fails later
                
                # skip sequences that have already been called 'N'
                if 'N' in lookups[seq]:
                    continue
                # check if contained
                if not hairpin_contained(hpin_locs, query):
                    lookups[seq] = 'N'
                else:
                    lookups[seq] = 'Y'

        with open(output.lookups, 'w') as f_out:
            for k,v in lookups.items():
                print(k,v, sep='\t', file=f_out)

#####################
# other annotations
#####################
# get repeat classes and convert to bed
rule repeats_to_bed:
    output:
        reps     = S3.remote('{bucket}/public/repeats/{release}/equCab3.fa.out'),
        reps_bed = S3.remote('{bucket}/public/repeats/{release}/equCab3.fa.bed'),
    params:
        work_dir     = lambda wildcards, output: os.path.dirname(output.reps_bed),
        rep_bed_name = lambda wildcards, output: os.path.basename(output.reps_bed),
        conda_env    = config['conda_envs']['novel_small']
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            cd {params.work_dir}
            wget http://hgdownload.soe.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.out.gz
            gunzip equCab3.fa.out.gz            

            # convert repeatmasker from ucsc to bed
            rmsk2bed < equCab3.fa.out > {params.rep_bed_name}
        '''

checkpoint separate_repeat_types:
    input:
        reps_bed = S3.remote('{bucket}/public/repeats/{release}/equCab3.fa.bed'),
    output:
        reps_dir = directory('{bucket}/public/repeats/{release}/beds'),
    run:
        # loop through bed file to get each type into list
        d_types = defaultdict(list)
        with open(input.reps_bed, 'r') as f_in:
            for line in f_in:
                line = line.strip().split()
                rm_type = line[-5]
                if '/' in rm_type:
                    rm_type = rm_type.split('/')[0]
                d_types[rm_type].append('\t'.join(line))                

        # write each type to bed
        os.makedirs(output.reps_dir, exist_ok=True)
        for k,v in d_types.items():
            if '?' in k:
                k = k.replace('?','_q')
            bed = os.path.join(output.reps_dir, f'{k}.eca3.bed')
            with open(bed, 'w') as f_out:  
                print('\n'.join(v), file=f_out)

rule create_class_islands:
    input:
        reps = '{bucket}/public/repeats/{release}/beds/{rep}.eca3.bed',
    output:
        merge_reps = S3.remote('{bucket}/public/repeats/{release}/merged/merged_{rep}.eca3.bed'),
    params:
        conda_env = config['conda_envs']['novel_small']
    threads: 1
    resources:
        time   = 10,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            bedtools merge -s \
                -d -1 \
                -c 4,6,11 \
                -o distinct \
                -delim "|" \
                -i {input.reps} \
                > {output.merge_reps}
        '''

def get_bed_reps(wildcards):
    # get repeat class bed based on wildcards
    checkpt_out = checkpoints.separate_repeat_types.get(**wildcards).output[0]
    # glob repeat class names
    REPS, = glob_wildcards(os.path.join(checkpt_out, '{rep}.eca3.bed'))
    # return class islands beds
    return expand(
        '{bucket}/public/repeats/{release}/merged/merged_{rep}.eca3.bed',
        bucket=config['bucket'], 
        release=[
           #'GCF_002863925.1_EquCab3.0',
            'Equus_caballus.EquCab3.0.103'
        ],
        rep=REPS
    )

localrules: combine_class_islands
rule combine_class_islands:
    input:
        get_bed_reps
    output:
        S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.bed'), 
    shell:
        '''
            cat {input} | sort -k 1,1 -k2,2n > {output}
        '''

# clean up the bed file names and include the score column ('.') in appropriate position
localrules: cleanup_combined_islands
rule cleanup_combined_islands:
    input:
        ref_ens = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/refseq_ens_table.tsv'),
        merged  = S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.bed'), 
    output:
        cleaned = S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.cleaned.bed')
    run:
        # get refseq ensembl conversion table into dict
        d_convert = pd.read_csv(input.ref_ens,sep='\t').set_index('refseq_mod')['ensembl_mod'].to_dict()
        
        with open(input.merged, 'r') as f_in, \
            open(output.cleaned, 'w') as f_out:
            for line in f_in:
                line = line.strip().split()
                if '/' in line[-1]:
                    line[-1] = line[-1].split('/')[0]
                # replace chrUn (refseq) with PJAA (ensembl)
                if line[0] in d_convert:
                    line[0] = d_convert[line[0]]
                # and chrM
                if 'chrM' in line[0]:
                    line[0] = 'chrMt'
                # set bed name column to repeat types
                line[3] = line[-1]
                # add score to position
                line.insert(4,'.')
                print('\t'.join(line), file=f_out) 

rule repeats_fa:
    input:
        ref_fa  = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna'),
        cleaned = S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.cleaned.bed')
    output:
        cleaned_fa = S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.cleaned.fa')
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 30,
        mem_mb = 6000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            
            bedtools getfasta \
                -fi {input.ref_fa} \
                -s \
                -name \
                -bed {input.cleaned} \
                -fo {output.cleaned_fa}
        '''

rule bowtie_index_repeats:
    input:
        S3.remote('{bucket}/public/repeats/{release}/merged_all_repeats.eca3.sorted.cleaned.fa')
    output:
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.2.ebwt',keep_local=True),
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.3.ebwt',keep_local=True),
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.4.ebwt',keep_local=True),
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.rev.1.ebwt',keep_local=True),
        S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa.rev.2.ebwt',keep_local=True),
        ebwt_base = S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa'),
        idx_dir   = directory('{bucket}/public/repeats/{release}/BOWTIE_INDICES')
    params:
        conda_env = config['conda_envs']['novel_small'],
    threads: 4
    resources:
        time   = 60,
        mem_mb = 16000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            # copy fa to index dir
            cp {input} {output.ebwt_base}
            # build bowtie index for repeat class islands
            bowtie-build {output.ebwt_base} {output.ebwt_base}
        '''

rule substrings_align_repeats:
    input:
        bt_idx    = rules.bowtie_index_repeats.output.idx_dir,
        ebwt_base = S3.remote('{bucket}/public/repeats/{release}/BOWTIE_INDICES/merged_all_repeats.eca3.sorted.cleaned.fa'),
        subs_fa   = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_substrings.fa'),
    output:
        sam    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_repeats.sam'),
        unmap  = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_repeats.unmap'),
        bt_log = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_repeats.log'),
    params:
        conda_env = config['conda_envs']['novel_small']
    threads: 6
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bowtie \
                -f \
                -a \
                -n 0 \
                -p 6 \
                --sam-nosq \
                --norc \
                {input.ebwt_base} \
                {input.subs_fa} \
                -S {output.sam} \
                --un {output.unmap} \
                2> {output.bt_log}
        '''

localrules: repeats_bundle
rule repeats_bundle:
    input:
        sam = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/hairpin_sub_repeats.sam'),
    output:
        repeats = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_repeats.txt'),
    run:
        # get substrings into dict with values being a list of unique repeat types
        repeats = defaultdict(list)
        with open(input.sam, 'r') as f_in:
            for line in f_in:
                if line.startswith('@'):
                    continue
                # get each alignment row
                row_tmp = line.strip().split()
                if row_tmp[1] == '4':
                    seq = row_tmp[9]
                    rep_class = ''
                else:
                    query = query_alignment(line.strip())
                    seq = Seq(line.strip().split()[9])
                    if query['flag'] == '16':
                        seq = seq.reverse_complement()
                    seq = str(seq)
                    rep_class = line.strip().split()[2].split(':')[0]
                # add to dict
                if rep_class not in repeats[seq]:
                    repeats[seq].append(rep_class)

        # write to file
        with open(output.repeats, 'w') as f_out:
            for k,v in repeats.items():
                print(k,','.join(v),sep='\t',file=f_out)

#####################
# meta coordinates
#####################
localrules: meta_coords_bundle
rule meta_coords_bundle:
    input:
        bed   = S3.remote('{bucket}/public/mirbase/v22/eca_db/eca3.ens_manual_novel.sorted.bed')
    output:
        coords = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_coords.txt'),
    run:
        # extract mature seqs (MIMAT and MIMATNO - novel mature) from bed
        # and prepare bundle file
        with open(input.bed, 'r') as f_in, \
            open(output.coords, 'w') as f_out:
            for line in f_in:
                if any([i in line for i in ['MIMAT','MIMATNO']]):
                    chrom,start,end,_,_,strand,_,_,_,name = line.strip().split()
                    name = name.split(';')
                    mimat = name[0].split('ID=')[1]
                    mir = name[2].split('Name=')[1]
                    print(
                        chrom, strand, 
                        int(start)+1, end, 
                        f'{mimat}&{mir}', 
                        sep='\t', file=f_out
                    )
        
#####################
# snp table
#####################
# NOTE - need to think about this one. the snp_vcf should be set in the config
# or required elsewhere for more flexibility
localrules: convert_snps_to_bed
rule convert_snps_to_bed:
    input:
        snp_vcf = S3.remote('Hf2hrt/public/refgen/Equus_caballus.EquCab3.0.103/MNEc2M.EquCab3.09182018.recode.rmGT.refseqnames.chrIDmod.vcf')
    output:
        snp_bed = S3.remote('Hf2hrt/public/refgen/Equus_caballus.EquCab3.0.103/MNEc2M.EquCab3.bed')
    run:
        # python to convert the standard snp vcf format to bed suitable for overlap with hsa.bed
        with open(input.snp_vcf, 'r') as f_in, \
            open(output.snp_bed, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                chrom = line[0]
                pos = int(line[1])
                snp_id = line[2]
                ref = line[3]
                alt = line[4]
                # if either multi or indel, skip
                if any(len(s) > 1 for s in [ref, alt]):
                    continue
                print(chrom, pos - 1, pos, ref, alt, snp_id, sep='\t', file=f_out)

rule intersect_snps_flanks:
    input:
        flanks  = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/eca3.ens_manual_novel.sorted_flanks.bed'),
        snp_bed = S3.remote('Hf2hrt/public/refgen/Equus_caballus.EquCab3.0.103/MNEc2M.EquCab3.bed')
    output:
        overlap = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/MNEc2M_hairpin.overlap')
    params:
        conda_env = config['conda_envs']['novel_small']
    threads: 6
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}

            bedtools intersect \
                -a {input.snp_bed} \
                -b {input.flanks} \
                -wa -wb > \
                {output.overlap}
        '''

# for each row in the intersection, check if snp is contained, and generate 
# dictionary with the key being the miRNA sequence containing the snp
localrules: snp_contained
rule snp_contained:
    input:
        space_fa = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
        overlap  = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/temp/MNEc2M_hairpin.overlap')
    output:
        snp_bundle = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_snps.txt'),
    run:
        tmp = SeqIO.to_dict(SeqIO.parse(input.space_fa, 'fasta'))
        # dict of hairpin ids (MIs) : seq
        hpin_seq = {}
        for k,v in tmp.items():
            k2 = k.split('|')[0]
            hpin_seq[k2] = str(v.seq)
  
        comps = {'A':'T', 'T':'A',
                 'G':'C', 'C':'G'}
        # for each overlap between snp location and hairpin plus flanks
        d_rows = {}
        with open(input.overlap, 'r') as f_in:
            for line in f_in:
                line = line.strip().split()
                if ',' in line[4] or '.' in line[4]: # skip multiallelic...for now
                    continue
                hpin_loci = line[-3]
                hpin_name = hpin_loci.split('|')[0]
                snp_name = line[5]
                pos = int(line[1]) - int(line[7])
                snp = line[4]
                if '-' in line[-1]:
                    pos = (int(line[8]) - 1) - int(line[1])
                    snp = comps[snp]
                hpin_mod = replace_snp(hpin_seq[hpin_name], pos, snp)
                for i in range(17,26):
                    tmp = snp_contained_kmers(hpin_seq[hpin_name], hpin_mod, i, pos)
                    for k,v in tmp.items():
                        # tuple of snp_name, hpin_loci, and ref seq in bundle format
                        row_item = (f'[{snp_name}]',f'[{hpin_loci}]',f'[{v}]')
                        if k not in d_rows:
                            d_rows[k] = [row_item]
                        else:
                            d_rows[k].append(row_item)

        # write to file
        with open(output.snp_bundle, 'w') as f_out:
            for k,v in d_rows.items():
                if len(v) > 1:
                    snp_ids, hpin_ids, ref_seqs = zip(*v)
                    print(
                        k,
                        ', '.join(snp_ids),
                        ', '.join(hpin_ids),
                        ', '.join(ref_seqs),
                        sep='\t', file=f_out 
                    )           
                else:
                    print(k, '\t'.join(str(i) for i in v[0]), sep='\t', file=f_out)

#####################
# all bundle
#####################
localrules: bundle_config
rule bundle_config:
    input:
        space_fa   = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_mirnaspace.fa'),
        lookups    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_lookups.txt'),
        repeats    = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_repeats.txt'),
        coords     = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_coords.txt'),
        snp_bundle = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/eca3_bundle_snps.txt'),
    output:
        tables_cfg = S3.remote('{bucket}/private/small/quant/isomirmap/{release}/bundle/tables.cfg')
    run:
        # generate the table config file
        with open(output.tables_cfg, 'w') as f_out:
            print(dedent('''\
                MAPPINGTABLES_NAME=miRBase_v22
                MAPPINGTABLES_REV=BundleRev_001
                MAPPINGTABLES_DESCRIPTION=miR-space: miRBase v22, assembly: EquCab3'''
            ), file=f_out)
            print(f'FRAGMENTS={name_and_hash(input.lookups)[0]} MD5SUM={name_and_hash(input.lookups)[1]}', file=f_out)
            print(f'HAIRPINSEQUENCES={name_and_hash(input.space_fa)[0]} MD5SUM={name_and_hash(input.space_fa)[1]}', file=f_out)
            print(f'OTHERANNOTATIONS={name_and_hash(input.repeats)[0]} MD5SUM={name_and_hash(input.repeats)[1]}', file=f_out)
            print(f'METACOORDINATES={name_and_hash(input.coords)[0]} MD5SUM={name_and_hash(input.coords)[1]}', file=f_out)
            print(f'SNPS={name_and_hash(input.snp_bundle)[0]} MD5SUM={name_and_hash(input.snp_bundle)[1]}', file=f_out)

 

