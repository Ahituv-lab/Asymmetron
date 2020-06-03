import re,os,sys,glob

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]	
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

chroms = range(1,23)+["X","Y"]
chromsL=["chr"+str(i) for i in chroms]
files_total=glob.glob("outs/MA*")
os.system("mkdir beds")
for path in files_total:
	name=path+"/fimo.tsv"
	if len(glob.glob(path+"/*tsv"))>0:
		DataL=reader(name)[1:]
		length=len(DataL[0])
		datafile=open("beds/"+path.split("/")[-1].split(".tsv")[0]+".bed","w")
		for i in DataL:
			if len(i)>1:
				datafile.write(i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[6]+'\t'+i[5]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\t'+i[0]+'\t'+i[1]+'\n')
		datafile.close()

