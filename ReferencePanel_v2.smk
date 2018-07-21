'''Snakefile for Construction of Reference Panel 0.1'''
# snakemake -s ReferencePanel_v2.smk
# snakemake -s ReferencePanel_v2.smk --dag | dot -Tsvg > dag_ReferencePanel_v2.svg

configfile: "ReferencePanel_config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

rule all:
    input:
        expand("DerivedData/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand("DerivedData/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map']),
        expand("DerivedData/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples']),
        "DerivedData/ReferencePanel/MtMap.txt",
        "DerivedData/ReferencePanel/MtStrand.txt"

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        "scripts/PYTHON/ambiguous2missing.py",
        "data/ReferencePanel_v2/McInerney_Master_Alignment_July18_2018.fasta",
    output:
        temp("DerivedData/ReferencePanel_v2/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"),
    params:
        in_fasta = "data/ReferencePanel_v2/McInerney_Master_Alignment_July18_2018.fasta",
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        out = "data/ReferencePanel_v2/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'

