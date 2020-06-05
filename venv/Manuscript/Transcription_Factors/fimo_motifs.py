import re,os,sys,glob

# The files variable contains the paths to all the transcription factor motifs in meme format
files=glob.glob("MEME_files_non_redundant/*.meme")

# Runs independently each motif as a job in LSF
motif = files[int(sys.argv[1])-1]
os.system("mkdir outs/"+motif.split("/")[-1])

# Performs FIMO motif finding
os.system("~/meme/bin/fimo -oc outs/"+motif.split("/")[-1]+"/ --bfile background_model.hg38 --thresh 1e-4 --max-stored-scores 10000000 "+motif +" hg38.fa")
