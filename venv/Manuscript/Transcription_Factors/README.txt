# Data paths
The human genome hg38 fasta file was downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

JASPAR motifs were downloaded from the following links:

For transcription factors:
http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.zip

for INR and TATA-box:
http://jaspar.genereg.net/api/v1/matrix/POL012.1.meme
http://jaspar.genereg.net/api/v1/matrix/POL002.1.meme

Genes were downloaded from:
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz


# Analysis

The background file was generated using the coomand: fasta-get-markov hg38.fa

Motif extraction was performed with FIMO using the script fimo_motifs.py


# Protein coding genes were filtered for "gene" and gene_type "protein_coding" using the command:
cat gencode.v33.annotation.gtf | grep 'gene_type "protein_coding"'  | awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t"  $7}'




