
rule ec2_nice_fasta:
    input:
       #fna = FTP.remote("http://ftp.ensembl.org/pub/release-94/fasta/equus_caballus/dna/Equus_caballus.EquCab2.dna.toplevel.fa.gz")
        fna = S3.remote("{bucket}/public/refgen/Equus_caballus.EquCab2.94/Equus_caballus.EquCab2.dna.toplevel.fa")
    output:
        nice_fna = S3.remote("{bucket}/public/refgen/Equus_caballus.EquCab2.94/Equus_caballus.EquCab2.94_genomic.nice.fna")
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    run:
        # open depending on extension
        class RawFile(object):
            def __init__(self,filename):
                self.filename = filename
                if filename.endswith('.gz'):
                    self.handle = gzip.open(filename,'rt')
                elif filename.endswith('bz2'):
                    self.handle = bz2.open(filename,'rt')
                elif filename.endswith('xz'):
                    self.handle = lzma.open(filenaem,'rt')
                else:
                    self.handle = open(filename,'r')
            def __enter__(self):
                return self.handle
            def __exit__(self,dtype,value,traceback):
                self.handle.close()

        # dictionary of ensembl ids to nice ids
        ensem_map = {}
        
        def inc_range(start, end):
            return range(start, end+1)
                                
        for i in inc_range(1,31):
            ensem_map[str(i)] = 'chr' + str(i)
            if 'MT' or 'X' not in ensem_map:
                ensem_map['MT'] = 'chrMt'
                ensem_map['X'] = 'chrX'

        # open fasta from ensembl and replace names
        with RawFile(input.fna) as IN, open(output.nice_fna,'w') as OUT:                                                                                                        
            for line in IN:                                        
                if line.startswith('>'):
                    name, *fields = line.lstrip('>').split()
                    if name in ensem_map:
                        new_name = '>' + ensem_map[name]
                        line = ' '.join([new_name, name] + fields + ['\n'])
                print(line,file=OUT,end='')

rule get_mirbase_fastas:
    output:
        mb_pre  = S3.remote("{bucket}/public/mirbase/v22/hairpin.fa",keep_local=True),
        mb_mat  = S3.remote("{bucket}/public/mirbase/v22/mature.fa",keep_local=True),
        eca_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_eca.fa",keep_local=True),
        eca_mat = S3.remote("{bucket}/public/mirbase/v22/mature_eca.fa",keep_local=True),
        bta_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_bta.fa",keep_local=True),
        bta_mat = S3.remote("{bucket}/public/mirbase/v22/mature_bta.fa",keep_local=True),
        oar_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_oar.fa",keep_local=True),
        oar_mat = S3.remote("{bucket}/public/mirbase/v22/mature_oar.fa",keep_local=True),
        ssc_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_ssc.fa",keep_local=True),
        ssc_mat = S3.remote("{bucket}/public/mirbase/v22/mature_ssc.fa",keep_local=True),
        mmu_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_mmu.fa",keep_local=True),
        mmu_mat = S3.remote("{bucket}/public/mirbase/v22/mature_mmu.fa",keep_local=True),
        hsa_pre = S3.remote("{bucket}/public/mirbase/v22/hairpin_hsa.fa",keep_local=True),
        hsa_mat = S3.remote("{bucket}/public/mirbase/v22/mature_hsa.fa",keep_local=True),
    params:
        conda_env = config['conda_envs']['small'],
        base      = "{bucket}/public/mirbase/v22"
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e
            
            # mirbase hairpin.fa
            wget https://mirbase.org/ftp/CURRENT/hairpin.fa.gz
            gunzip hairpin.fa.gz
            # mirbase mature.fa
            wget https://mirbase.org/ftp/CURRENT/mature.fa.gz
            gunzip mature.fa.gz

            # horse
            extract_miRNAs.pl hairpin.fa eca > {output.eca_pre}
            extract_miRNAs.pl mature.fa eca > {output.eca_mat}

            # cow
            extract_miRNAs.pl hairpin.fa bta > {output.bta_pre}
            extract_miRNAs.pl mature.fa bta > {output.bta_mat}
            
            # sheep
            extract_miRNAs.pl hairpin.fa oar > {output.oar_pre}
            extract_miRNAs.pl mature.fa oar > {output.oar_mat}
            
            # pig
            extract_miRNAs.pl hairpin.fa ssc > {output.ssc_pre}
            extract_miRNAs.pl mature.fa ssc > {output.ssc_mat}
            
            # mouse
            extract_miRNAs.pl hairpin.fa mmu > {output.mmu_pre}
            extract_miRNAs.pl mature.fa mmu > {output.mmu_mat}
            
            # human
            extract_miRNAs.pl hairpin.fa hsa > {output.hsa_pre}
            extract_miRNAs.pl mature.fa hsa > {output.hsa_mat}

            # move all to base dir
            mv hairpin.fa {output.mb_pre}
            mv mature.fa {output.mb_mat}

        '''

rule combine_other_mature:
    input:
        bta_mat = S3.remote("{bucket}/public/mirbase/v22/mature_bta.fa",keep_local=True),
        oar_mat = S3.remote("{bucket}/public/mirbase/v22/mature_oar.fa",keep_local=True),
        ssc_mat = S3.remote("{bucket}/public/mirbase/v22/mature_ssc.fa",keep_local=True),
        mmu_mat = S3.remote("{bucket}/public/mirbase/v22/mature_mmu.fa",keep_local=True),
        hsa_mat = S3.remote("{bucket}/public/mirbase/v22/mature_hsa.fa",keep_local=True),
    output:
        other_mat = S3.remote("{bucket}/public/mirbase/v22/mature_other_species.fa",keep_local=True),
    shell:
        '''
            cat {input} > {output}
        '''

rule get_eca_mirbase:
    output:
        eca_pre  = S3.remote("{bucket}/public/mirbase/v22/eca_db/hairpin.fa",keep_local=True),
        eca_mat  = S3.remote("{bucket}/public/mirbase/v22/eca_db/mature.fa",keep_local=True),
        eca_str  = S3.remote("{bucket}/public/mirbase/v22/eca_db/miRNA.str",keep_local=True),
        eca_gff  = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca2.ncbi.gff",keep_local=True)
    params:
        conda_env = config['conda_envs']['small'],
        base      = "{bucket}/public/mirbase/v22"
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e
            
            # mirbase hairpin.fa
            wget https://mirbase.org/ftp/CURRENT/hairpin.fa.gz
            gunzip hairpin.fa.gz
            
            # mirbase mature.fa
            wget https://mirbase.org/ftp/CURRENT/mature.fa.gz
            gunzip mature.fa.gz
            
            # mirbase str
            wget https://mirbase.org/ftp/CURRENT/miRNA.str.gz
            gunzip miRNA.str.gz

            # mirbase gff
            wget https://www.mirbase.org/ftp/CURRENT/genomes/eca.gff3

            # horse
            extract_miRNAs.pl hairpin.fa eca > {output.eca_pre}
            extract_miRNAs.pl mature.fa eca > {output.eca_mat}

            # move all to base dir
            mv miRNA.str {output.eca_str}
            mv eca.gff3 {output.eca_gff}
        '''

# chrom alias originally came from UCSC
#https://hgdownload-test.gi.ucsc.edu/hubs/GCA/000/002/305/GCA_000002305.1/GCA_000002305.1.chromAlias.txt
rule modify_mirbase_gff:
    input:
        eca_gff  = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca2.ncbi.gff"),
        alias    = S3.remote("{bucket}/public/refgen/GCA_000002305.1/GCA_000002305.1.chromAlias.txt"),
        ncbi_ec2 = S3.remote("{bucket}/public/refgen/GCA_000002305.1/GCA_000002305.1_EquCab2.0_genomic.dict"),
        ens_ec2  = S3.remote("{bucket}/public/refgen/Equus_caballus.EquCab2.94/Equus_caballus.EquCab2.94_genomic.nice.dict")
    output:
        ens_eca_gff     = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca2.ens.gff"),
        ens_eca_gff_mod = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca2_mod.ens.gff"),
        feat_types      = S3.remote("{bucket}/public/mirbase/v22/eca_db/feat_types.txt")
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    run:
        # read chrom alias, and ec2 dicts to dataframes
        # alias df
        df_alias = pd.read_csv(
            input.alias,
            sep='\t',
            # index_col="refseq",
            usecols=["refseq","# genbank"]
        )
        # ens ec2 df
        df_ens = pd.read_csv(
            input.ens_ec2,
            sep="\t",
            header=None,
            skiprows=1,
            usecols=[1,3],
            names=['ens_sn','m5']
        )
        df_ens['ens_sn'] = df_ens['ens_sn'].str.replace('SN:','')

        # ncbi ec2 df
        df_ncbi = pd.read_csv(
            input.ncbi_ec2,
            sep="\t",
            header=None,
            skiprows=1,
            usecols=[1,3],
            names=['ncbi_sn','m5']
        )

        df_ncbi['ncbi_sn'] = df_ncbi['ncbi_sn'].str.replace('SN:','')

        # combine ncbi and ens by m5 id
        ncbi_ens = pd.merge(df_ens,df_ncbi,on='m5',how='outer')

        # merge alias using genbank and ncbi_sn
        df_full = df_alias.merge(
            ncbi_ens,
            left_on='# genbank',
            right_on='ncbi_sn',
            # how='outer'
        ).set_index('refseq').drop('# genbank',axis=1)
        # NOTE that NCBI ec2 is missing chrMt

        # get genbank and ensembl sequence names
        refseq = df_full.index.values.tolist()
        ens_sn = df_full.ens_sn.values.tolist()

        # open mirbase gff and rename refseq with ensembl names
        with open(input.eca_gff,'r') as infile, open(output.ens_eca_gff,'w') as outfile:
            for line in infile:
                for ref,ens in zip(refseq,ens_sn):
                    line = line.replace(ref,ens)
                outfile.write(line)

        # modify eca2.ens.gff to work with liftoff due to known issue with 
        # lifting over short sequences - hacky but works and fix has not yet
        # (20221114) been implemented
        body = []
        modify = []

        with open(output.ens_eca_gff,'r') as infile:
            for line in infile:
                if "Derives_from" in line:
                    line = line.replace("Derives_from","Parent")
                body.append(line)
                if (not line.startswith("#")) and ("miRNA_primary_transcript" in line):
                    for r in (("ID","Parent"), ("miRNA_primary_transcript","dummy")):
                        line = line.replace(*r)
                    modify.append(line)
 
        with open(output.ens_eca_gff_mod,'w') as f:
            for line in body+modify:
                f.write(line)

        # write feature types to file for liftoff
        with open(output.feat_types,'w') as f:
            print("miRNA\nmiRNA_primary_transcript",file=f)

rule ens_eca2_liftover_eca3:
    input:
        eca_gff_mod  = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca2_mod.ens.gff"),
        feat_types   = S3.remote("{bucket}/public/mirbase/v22/eca_db/feat_types.txt"),
        nice_ec2_fna = S3.remote("{bucket}/public/refgen/Equus_caballus.EquCab2.94/Equus_caballus.EquCab2.94_genomic.nice.fna"),
        nice_ec3_fna = S3.remote("{bucket}/public/refgen/Equus_caballus.EquCab3.0.103/Equus_caballus.EquCab3.0.103_genomic.nice.fna"),
    output:
        eca3_gff = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca3.ens.gff"),
        unmapped = S3.remote("{bucket}/public/mirbase/v22/eca_db/unmapped_features.txt"),
    params:
        conda_env = config['conda_envs']['small'],
        feats     = S3.remote("{bucket}/public/mirbase/v22/eca_db/feat_types.txt"),
        tmp_gff   = "{bucket}/public/mirbase/v22/eca_db/eca3_tmp.gff",
    threads: 1
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
 
            liftoff \
                -g {input.eca_gff_mod} \
                -f {input.feat_types} \
                -o {params.tmp_gff} \
                {input.nice_ec3_fna} \
                {input.nice_ec2_fna}
        
            # remove dummy lines
            awk '!/dummy/' {params.tmp_gff} > {output.eca3_gff}

            # unmapped features saved to working directory - move to output
            mv unmapped_features.txt {output.unmapped}
        '''

rule modify_lifted_gff:
    input:
        eca3_gff = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca3.ens.gff"),
    output:
        eca3_mirtop_gff = S3.remote("{bucket}/public/mirbase/v22/eca_db/eca3.ens_mirtop.gff"),
    run:
        # in order to use miRTop, need to modify the lifted gff by including
        # mirbase version in header and reverting "Parent" to "Derives_from"
        body = []

        with open(input.eca3_gff,'r') as infile:
            for line in infile:
                if "Parent" in line:
                    line = line.replace("Parent","Derives_from")
                body.append(line)
                # to get mirbase version into gff
                if "##gff-version 3" in line:
                    body.append("# miRBase v22\n")
 
        with open(output.eca3_mirtop_gff,'w') as f:
            for line in body:
                f.write(line)

# there are 13 unmapped features, only two of which (MI0012799 and MI0028353)
# are included in the ensembl 103 (or 109) annotation. these were manually
# added to eca3.ens_mirtop.gff (and sorted), and uploaded as 
# eca3.ens_manual.gff


