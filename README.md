# Lake St. Clair Microseira Wollei - 2023, 2024
'Omics analysis pipeline for M. wollei HABs study in Lake St. Clair

### Approach:
- Analysis utilized the Great Lakes Atlas for Multi-omics Research analysis products (GLAMR - https://greatlakesomics.org)
- R Markdown files direct the workflow 
- Snakemake file can be called through R Markdown code, or command line

### Data Availability:
Raw sequence files associated with this study can be accessed at the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) under BioProject PRJNA1474428. 

### Notebooks:

**SC_notebook1.Rmd**
  - Identify the presence of saxitoxin/cyanotoxin genes in metagenomic data using BLASTN
  - Refine to only include metagenome amplified genomes (>90% completeness, <10% contamination) with saxitoxin genes
  - Examine placement of Lake St. Clair MAG of interest in phylogenomic tree
  - Extract contigs containing saxitoxin for visualization via clinker
  - Use coverm to characterize MAG of interest
    
