import os
import re
import csv
import gzip
import json
import shutil
import random
import hashlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict

localrules: gather_se_adapter_logs,
            gather_post_join_logs,
            filter_mirna_nc_ref,
            fastq_lengths, 
            modify_lifted_gff,
            bowtie_index_ncrna_ref,
            bowtie_index_ncrna_rfam,
            bowtie_index_eca_hairpin,
            gather_fastp_jsons,
            filter_samples_reads,
            filter_samples_biotypes,
            combine_filters,
            convert_unaligned_to_fasta,
            get_mirbase_fastas,
            save_mirdeep2_counts,
            combine_other_mature,
            save_mirna_novel_mirdeep2,
            upset_novel_mirna,
            decomp_clean_fqs,
            mirtrace_pre_process_csv,
            mirtrace_post_process_csv,
            distinguish_exp_mirna_samples,
            filter_blast_rnacent,
            reformat_merged_bed,
            get_novel_candidates,
            filter_by_norm_exp,
            get_exp_precur_bed,
            remove_precurs_overlap_ref,
            final_precur_mature_beds,
            refseq_ensembl_table,            
            hairpin_flanks,
            mirnaspace_bundle,
            hairpin_substrings,
            lookups_bundle,
            separate_repeat_types,
            combine_class_islands,
            cleanup_combined_islands,
            repeats_bundle,
            meta_coords_bundle,
            convert_snps_to_bed,
            snp_contained,
            bundle_config,
            isomirmap_filter,
            isomirmap_descript,
            minimum_tissue_count,
            isomirmap_downsample_filter,
            isomirmap_downsample_descript,


singularity: config['sif']
include: "src/utils.py"

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id,
    secret_access_key=s3_access_key
)

# small
units_small = pd.read_table(config['samples'], dtype=str) \
    .set_index(['sample', 'tissue'], drop=False) \
    .sort_index()

num_tissues = len(units_small['tissue'].unique().tolist())

rule all:
    input:
       #######################
       ## MAIN PIPELINE
       #######################
       #S3.remote(
       #    expand(
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.pre_trim.fastp.html',
       #        '{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.trim.fastp.html',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}.clean.fastp.html',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_join.clean.fastq.gz',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mirna.unaligned_mm0{mm}.fastq',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/{u.tissue}_{u.sample}_{u.lane}_join.sized.clean.freq.tsv',
       #       #'{bucket}/public/refgen/{release}/BBMAP_INDICES/GENOMIC/ref/genome/1/info.txt',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.corrected.mm00.fastq',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}_join.clean.mm00_ncref_rfam.unaligned.freq.tsv',
       #       #'{bucket}/private/small/quant/miraligner/{u.tissue}/{u.sample}/{u.tissue}_{u.sample}.mirna',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #       #mm=list(range(3)),
       #        mm=['0'] # NOTE - only include mismatch of 0 for now
       #    )
       #),
       # PRE-PROCESS FILTERING PRIOR TO QUANT
       #expand(
       #    '{bucket}/private/small/quant/preproc_all_filters.all_samples.tsv',
       #    bucket=config['bucket'], 
       #   #fltr=['preproc_read_counts', 'preproc_biotype_props']
       #),
       ## mirdeep2 canon counts 
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/quant/mirdeep2/canon/combine/canonical_counts.tsv',
       #        bucket=config['bucket'], 
       #    )
       #),
       ## miraligner counts
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/quant/miraligner/{u.tissue}/{u.sample}/mirtop/counts/mirtop.tsv',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
       ## final miralign matrix
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/quant/miraligner/combine/isomir_counts.DUPE_TEST.csv',
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
       ## srnabench - PAUSE ON THIS TILL I CAN FIND v1.5
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/quant/srnabench/{u.tissue}/{u.sample}/mirtop/mirtop.gff',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
       ## mirpro quant and novel
       #S3.remote(
       #    expand(
       #       #'{bucket}/private/small/quant/mirpro/{release}/multi_sample/result/novel_quantifier_report.csv',
       #        '{bucket}/private/small/novel/{release}/novel_mirna.sorted.bed',
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
       ## isomirmap bundle and quant
       #S3.remote(
       #    expand(
       #       #'{bucket}/private/small/quant/isomirmap/{release}/temp/refseq_ens_table.tsv',
       #       #'{bucket}/private/small/quant/isomirmap/{release}/bundle/tables.cfg',
       #       #'{bucket}/private/small/quant/isomirmap/{release}/{u.tissue}/{u.sample}/{u.tissue}_{u.sample}-IsoMiRmap_v5-isomiRs.miRBase.gff3',
       #       #'{bucket}/private/small/quant/isomirmap/{release}/merge/isomir_{exp}.isomirmap_exp_filt.tsv',
       #        '{bucket}/private/small/quant/isomirmap/{release}/merge/isomir_descript.tsv',
       #       #u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        exp=['cts','rpm']
       #    )
       #),
       ## FINAL STEP - DON'T DELETE
       ## isomirmap with or without down sampling
       #S3.remote(
       #    expand(
       #       #'{bucket}/private/small/quant/minimum_tissue_count.inc_samples.tsv',
       #       #'{bucket}/private/fastq/trimmed/{u.tissue}/{u.sample}/small/noncode_filt/{release}/{u.tissue}_{u.sample}.mirna_sampled.unaligned_mm00.fastq',
       #       #'{bucket}/private/small/quant/downsample/{release}/merge/isomir_{exp}.isomirmap.tsv',
       #        '{bucket}/private/small/quant/isomirmap/{release}/merge/isomir_descript.tsv',
       #       #'{bucket}/private/small/quant/downsample/{release}/merge/isomir_descript.ds.tsv',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        exp=['cts','rpm']
       #    )
       #),
        S3.remote(
            expand(
                '{bucket}/private/fastq/process_fastqs_counts.csv',
                bucket=config['bucket']
            )
        ),
       ## PER SAMPLE FEATURE COUNTS - LEAVE OPEN UNTIL MULTIQC RULE WORKS
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/align/bowtie/{release}/{u.tissue}/{u.sample}/mismatch_0{mm}/{seq_set}/{u.tissue}_{u.sample}.mm0{mm}_{seq_set}.ref.featcts.txt',
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        mm='0',
       #       #mm=list(range(3)),
       #        seq_set=['pre_filter','post_filter']
       #    )
       #),
       #S3.remote(
       #    expand(
       #        '{bucket}/private/small/align/bowtie/{release}/all_samples.featcts.csv',
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #),
       ## reshuffle
       #S3.remote(
       #    expand(
       #       #"{bucket}/private/small/quant/miraligner/{release}/combine/shuffle/isomir_merge.shuffle_{perm}.csv",
       #        "{bucket}/private/small/quant/miraligner/{release}/combine/shuffle/shuffle_all.csv",
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #       #perm=[str(i).zfill(4) for i in range(1,11)]
       #    )
       #),
       #######################
       ## QC
       #######################
       # SOMETHING WEIRD WITH MULTIQC AND MISSING OUTPUT FILES THAT DO IN FACT
       # EXIST...NOTE - MIRTRACE OUTPUT NEEDS TO BE MODED TO QC.SMK
       #S3.remote(
       #    expand(
       #        '{bucket}/private/qc/{release}/small/mirtrace/{step}_process/mirtrace-report.html',                
       #        u=units_small.itertuples(), 
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #        step=['pre','post']
       #    )
       #),
       #S3.remote(
       #    expand(
       #        '{bucket}/private/qc/{release}/small/multiqc_report.html',
       #        bucket=config['bucket'], 
       #        release=[
       #           #'GCF_002863925.1_EquCab3.0',
       #            'Equus_caballus.EquCab3.0.103'
       #        ],
       #    )
       #)

# ----------------------------------------------------------
#       various plots
# ----------------------------------------------------------

#include: 'rules/viz.smk'

# ----------------------------------------------------------
#       small rna processing and quant
# ----------------------------------------------------------

include: 'rules/make_nice.smk'
include: 'rules/process_small.smk'
#include: 'rules/size.smk'
##include: 'rules/mirtrace.smk'
#include: 'rules/featcts_small.smk'
#include: 'rules/filters.smk'
##include: 'rules/miraligner.smk'
#include: 'rules/mirpro.smk'
#include: 'rules/novel.smk'
#include: 'rules/mapping_bundle.smk'
#include: 'rules/isomirmap.smk'
#include: 'rules/isomirmap.down_sample.smk'

#include: 'rules/prost_NEW.smk'
#include: 'rules/mirdeep2.smk'
#include: 'rules/srnabench.smk'

#include: 'rules/known_mirna.smk'
#include: 'rules/novel_mirna.smk'
#include: 'rules/quant.known_mirna.smk'
#include: 'rules/quant.novel_mirna.smk'

# ----------------------------------------------------------
#       QC
# ----------------------------------------------------------

include: 'rules/qc.smk'

