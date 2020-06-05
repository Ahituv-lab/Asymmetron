# Urls for links to download data used in the study
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


Unibind TFBSs for PWMs and peak-caller MACS were downloaded from:
https://unibind.uio.no/static/data/bulk/pwm_tfbs.tar.gz


# Analysis

Step 1: Extraction of transcription factor binding sites for the JASPAR non-redundant motif list
-The background file was generated using the command: fasta-get-markov hg38.fa

-Motif extraction was performed with FIMO using the script fimo_motifs.py which is provided in this directory.

-The genome-wide motif maps were formmatted with the script: MEME_motifs_formatter.py


# Step 2: Identification of transcriptional strand asymmetries for each of the transcription factor genome-wide motif maps using FIMO.
Protein coding genes were filtered for "gene" and gene_type "protein_coding" using the command:
cat gencode.v33.annotation.gtf | grep 'gene_type "protein_coding"'  | awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t"  $7}' > gencode.v33.annotation.bed.protein_coding

The transcriptional strand asymmetry of each transcription factor was calculated with the command:
python contained_asymmetries.py gencode.v33.annotation.bed.protein_coding TF_X.bed


Promoter upstream and promoter downstream regions were extracted with the script get_promoters.py

python contained_asymmetries.py promoters/gencode.v33.annotation.bed.promoters_1kB_upstream TF_X.bed
python contained_asymmetries.py promoters/gencode.v33.annotation.bed.promoters_1kB_downstream TF_X.bed


Step 3: Filtering of the Unibind TFBS files was performed with the script unibind_filtering.py
Command to extract transcriptional strand asymmetry of transcription factor X in a single ChiP-seq experiment across transcribed regions:
python contained_asymmetries.py gencode.v33.annotation.bed.protein_coding Unibind_TF_X.bed


Command to extract transcriptional strand asymmetry of transcription factor X in a single ChiP-seq experiment across transcribed regions:
python contained_asymmetries.py promoters/gencode.v33.annotation.bed.promoters_1kB_upstream Unibind_TF_X.bed
python contained_asymmetries.py promoters/gencode.v33.annotation.bed.promoters_1kB_downstream Unibind_TF_X.bed


Command for pairwise asymmetries between TATA and INR motifs within 1kB around the TSS of protein-coding genes:
cat promoters/gencode.v33.annotation.bed.promoters_1kB_upstream >> promoters/gencode.v33.annotation.bed.promoters_1kB_window
cat promoters/gencode.v33.annotation.bed.promoters_1kB_downstream >> promoters/gencode.v33.annotation.bed.promoters_1kB_window

bedtools intersect -a TATA_box.bed -b promoters/gencode.v33.annotation.bed.promoters_1kB_window -u > TATA_box.bed.at_promoters
bedtools intersect -a INR_box.bed -b promoters/gencode.v33.annotation.bed.promoters_1kB_window -u > INT_box.bed.at_promoters
python pairwise_asymmetries.py TATA_box.bed.at_promoters INR.bed.at_promoters --min_distance=1 --max_distance=100 --bins=10 --plots


