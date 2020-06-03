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

Step 1: Extraction of transcription factor binding sites for the JASPAR non-redundant motif list
The background file was generated using the coomand: fasta-get-markov hg38.fa

Motif extraction was performed with FIMO using the script fimo_motifs.py
The genome-wide motifs found were formmatted with the script: MEME_motifs_formatter.py


# Step 2: Identification of transcriptional strand asymmetries for each of the transcription factors
Protein coding genes were filtered for "gene" and gene_type "protein_coding" using the command:
cat gencode.v33.annotation.gtf | grep 'gene_type "protein_coding"'  | awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t"  $7}' > gencode.v33.annotation.bed.protein_coding

Promoter upstream and promoter downstream regions were extracted with the script get_promoters.py






