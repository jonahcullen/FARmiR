
# upset plots for novel mirnas by tissue
localrules: upset_novel_mirna
rule upset_novel_mirna:
    input:
        full = S3.remote("{bucket}/private/small/quant/mirdeep2/{release}/novel/{tissue}/merge/mismatch_00/all/{tissue}_mirna_all.novel.csv"),
    output:
        upset = S3.remote("{bucket}/private/small/quant/mirdeep2/{release}/novel/{tissue}/merge/mismatch_00/all/{tissue}_mirna_novel.upset.png")
    shell:
        'Rscript src/upset_mirna.R {input.full} {output.upset}'




