'''Snakefile for Construction of Reference Panel 0.1'''
# snakemake -s ReferencePanel_v2.smk
# snakemake -s ReferencePanel_v2.smk --dag | dot -Tsvg > dag_ReferencePanel_v2.svg

configfile: "ReferencePanel_config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        "scripts/PYTHON/ambiguous2missing.py",
        "/g/data1a/te53/MitoImpute/data/FASTA/masters/McInerney_Master_Alignment_July18_2018.fasta",
    output:
        temp("/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"),
    params:
        in_fasta = "/g/data1a/te53/MitoImpute/data/FASTA/masters/McInerney_Master_Alignment_July18_2018.fasta",
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        out = "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"
    shell:
        'if [ -f {params.out} ]
        then
        	echo "FOUND... PASSING"
        else
        	python {params.in_script} -i {params.in_fasta} -o {params.out} -v
        fi
        '
        
## 2. Run the fasta2vcf_mtDNA.py script
rule fasta2vcf:
    input:
        "scripts/PYTHON/fasta2vcf_mtDNA.py",
        "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"
    output:
        "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz"
    params:
        in_fasta = "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta",
        in_script = "scripts/PYTHON/fasta2vcf_mtDNA.py",
        out = "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'

#$ <FILL THIS OUT ONCE COMMAND-LINE VERSIONS ARE WRITTEN>