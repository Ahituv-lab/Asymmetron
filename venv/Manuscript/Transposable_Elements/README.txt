# Data paths
The human genome hg38 fasta file was downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

The RepeatMasker dataset was downloaded from:
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz
LINE, SINE and LTR repeats were extracted using the following command:
cat rmsk.txt | grep  -E 'LINE|SINE|LTR' | awk '{print $6 "\t" $7 "\t" $8 "\t" $11 "\t"  $10 "\t" $12 "\t" $13}'


