import os
import re
import snakemake.io
from glob import glob

#Run BLASTN for cyanotoxin genes on all GLAMR (project set) samples (assemblies) of interest
rule make_blast_db:
    input: 
        "blast/database/nucl/{blast_name}.fasta"
    output: 
        "blast/database/nucl/{blast_name}.fasta.nin"
    log: 
        "blast/logs/database/{blast_name}.txt"
    conda: "config/blast.yml"
    resources: 
        cpus=1, mem_mb=5000, time_min=10000
    shell: 
        """
        makeblastdb -in {input} -dbtype nucl -logfile {log}
        """

rule blastn_assemblies_UNIVERSAL:
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        assembly = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/StClair_manuscript/blast/query/{run}/{sample}_final.contigs.renamed.fa"
    output:
         "blast/output/{target}/blastn/assemblies/{sample}__{run}_{blast_name}_contigs_blastn.tsv"
    log:
        "blast/logs/output/{target}/blastn/assemblies/{sample}__{run}_{blast_name}_contigs_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.assembly} -out {output} -outfmt '6 std qcovs stitle qlen' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_assemblies_UNIVERSAL_cyanoToxin:
    input: 
        expand("blast/output/cyanoToxin/blastn/assemblies/{sample}__{run}_{blast_name}_contigs_blastn.tsv", sample = glob_wildcards("/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/StClair_manuscript/blast/query/assemblies_vs_cyanoToxin/{sample}_final.contigs.renamed.fa",followlinks=True).sample, run = "assemblies_vs_cyanoToxin", blast_name = "blastdbcyanoToxin")

rule blastn_GToTree:
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        assembly = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/StClair_manuscript/blast/query/GToTree_vs_cyanoToxin/{sample}.fna"
    output: 
         "blast/output/{target}/blastn/mags/{sample}__{run}_{blast_name}_contigs_blastn.tsv"
    log:
        "blast/logs/output/{target}/blastn/mags/{sample}__{run}_{blast_name}_contigs_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.assembly} -out {output} -outfmt '6 std qcovs stitle qlen' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_GToTree_cyanoToxin:
    input: 
        expand("blast/output/cyanoToxin/blastn/mags/{sample}__{run}_{blast_name}_contigs_blastn.tsv", sample = glob_wildcards("/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/StClair_manuscript/blast/query/GToTree_vs_cyanoToxin/{sample}.fna",followlinks=True).sample, run = "GToTree_vs_cyanoToxin", blast_name = "blastdbcyanoToxin")

#Assign taxonomy to bins with cyanotoxin blastn hits


# Rule to run Prokka on cyanotoxin contigs
rule prokka_cyanoToxin_contigs:
    input: 
        fasta = "clinker/contig_compare/contigs_oi/{contig_name}.fa"
    params:
        prefix = "{contig_name}"
    output: 
        gbk = "clinker/contig_compare/contigs_oi/prokka/{contig_name}/{contig_name}.gbk",  # Specify the final .gbk output file
        gb = "clinker/contig_compare/contigs_oi/prokka/gb/{contig_name}.gb",
        folder = directory("clinker/contig_compare/contigs_oi/prokka/{contig_name}")
    conda: 
        "config/prokka_2024.1.yml"
    resources:
        cpus = 8, 
        mem_mb = 16000, 
        time_min = 10000
    shell:
        """
        prokka {input.fasta} \
               --prefix {params.prefix} \
               --outdir {output.folder} \
               --cpus {resources.cpus} \
               --force

        cp {output.gbk} {output.gb}
        """

rule run_prokka_cyanoToxin_contigs:
    input:
        expand("clinker/contig_compare/contigs_oi/prokka/{contig_name}/{contig_name}.gbk", contig_name = glob_wildcards("clinker/contig_compare/contigs_oi/{contig}.fa",followlinks=True).contig)
