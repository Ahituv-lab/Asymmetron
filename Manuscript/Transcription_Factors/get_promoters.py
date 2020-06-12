import re,os,sys,glob

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data


DataL=reader("gencode.v33.annotation.bed.protein_coding")
datafile=open("promoters/gencode.v33.annotation.bed.promoters_1kB_upstream","w")
datafile2=open("promoters/gencode.v33.annotation.bed.promoters_1kB_downstream","w")
for k in DataL:

	orientation = k[-1]

	# Promoter upstream regions
	# If the gene is oriented in plus oreintation
	if orientation=="+":
		datafile.write(k[0]+'\t'+str(max(0,int(k[1])-1000))+'\t'+str(int(k[1]))+'\t'+k[3]+'\t'+"."+'\t'+orientation+'\n')
	# If the gene is oriented in minus orientation
	elif orientation=="-":
		datafile.write(k[0]+'\t'+str(int(k[2]))+'\t'+str(int(k[2])+1000)+'\t'+k[3]+'\t'+"."+'\t'+orientation+'\n')

	# Promoter downstream regions
	# If the gene is oriented in plus oreintation
	if orientation=="+":
		datafile2.write(k[0]+'\t'+str(int(k[1]))+'\t'+str(int(k[1])+1000)+'\t'+k[3]+'\t'+"."+'\t'+orientation+'\n')
	# If the gene is oriented in minus orientation
	elif orientation=="-":
		datafile2.write(k[0]+'\t'+str(max(0,int(k[2])-1000))+'\t'+str(str(k[2]))+'\t'+k[3]+'\t'+"."+'\t'+orientation+'\n')
datafile.close()
datafile2.close()
