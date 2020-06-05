import re,os,sys,glob
import numpy as np

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def size(path):
        DataL=reader(path)
        Sizes=[]
        for v in DataL:
                Sizes+=[int(v[2])-int(v[1])]
        return Sizes

# Genic bins (We use 10 bins here)
Data_plus=reader("gencode.v33.annotation.bed.protein_coding.formatted.plus")
Data_minus=reader("gencode.v33.annotation.bed.protein_coding.formatted.minus")
bins=10
for binned in range(1,bins+1):
        datafile1=open("gencode.v33.annotation.bed.protein_coding.formatted_bin_"+str(binned)+".txt","w")

	# For plus oriented genes
        for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
                datafile1.write(chrom+'\t'+str(start_bin)+'\t'+str(end_bin)+'\t'+str(size)+'\t'+str(v[4])+'\n')
        
	#For minus oriented genes
        for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
                datafile1.write(chrom+'\t'+str(end_bin)+'\t'+str(start_bin)+'\t'+str(size)+'\t'+str(v[4])+'\n')
datafile1.close()

# Upstream 1kB
datafile1=open("gencode.v33.annotation.bed.protein_coding.formatted_upstream_1kB","w")
# For plus oriented genes
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
                if start-1000>0:
                        datafile1.write(chrom+'\t'+str(start-1000)+'\t'+str(start)+'\t'+str(size)+'\t'+str(v[4])+'\n')

# For minus oriented genes
for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
                datafile1.write(chrom+'\t'+str(end)+'\t'+str(end+1000)+'\t'+str(size)+'\t'+str(v[4])+'\n')
datafile1.close()


# Downstream 1kB
datafile1=open("gencode.v33.annotation.bed.protein_coding.formatted_downstream_1kB","w")
# For plus oriented genes
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
                datafile1.write(chrom+'\t'+str(end)+'\t'+str(end+1000)+'\t'+str(size)+'\t'+str(v[4])+'\n')

# For minus oriented genes
for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
                if start-1000>0:
                        datafile1.write(chrom+'\t'+str(start-1000)+'\t'+str(start)+'\t'+str(size)+'\t'+str(v[4])+'\n')
datafile1.close()

