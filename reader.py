import re,os,sys,glob

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

DataL=reader("MCF7_RepliStrand.leading_lagging")
datafile=open("MCF7_RepliStrand.leading_lagging.bed","w")
for i in DataL:
	if i[3]=="Leading":
		datafile.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+"+"+'\t'+i[-2]+'\n')
	elif i[3]=="Lagging":
		datafile.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+"-"+'\t'+i[-2]+'\n')
datafile.close()
