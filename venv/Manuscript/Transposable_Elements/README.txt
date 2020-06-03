# Data paths
The human genome hg38 fasta file was downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

The RepeatMasker dataset was downloaded from:
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz
LINE, SINE and LTR repeats were extracted using the following command:
cat rmsk.txt | grep  -E 'LINE|SINE|LTR' | awk '{print $6 "\t" $7 "\t" $8 "\t" $11 "\t"  $10 "\t" $12 "\t" $13}' > rmsk.txt.formmated

Transposable elements were extracted with the commands:
cat rmsk.txt.formmated | grep LINE > LINEs.bed
cat rmsk.txt.formmated | grep SINE > SINEs.bed
cat rmsk.txt.formmated | grep LTR > LTRs.bed

The transposable sufamilies were derived with the commands:
cat LINEs.bed  | awk '$7=="L1" {print $0}'
cat LINEs.bed  | awk '$7=="L2" {print $0}'
cat SINEs.bed  | awk '$7=="Alu" {print $0}'
cat SINEs.bed  | awk '$7=="MIR" {print $0}'

Structural variants were extracted from: 
https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
Breakpoints were formatted with the command:
cat gnomad_v2.1_sv.sites.bed |  awk '{print "chr"$1 "\t" $2 "\t" $2 "\t" $5}' >> gnomad_v2.1_sv.sites.bed.formmated
cat gnomad_v2.1_sv.sites.bed |  awk '{print "chr"$1 "\t" $3 "\t" $3 "\t" $5}' >> gnomad_v2.1_sv.sites.bed.formmated



