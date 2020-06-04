import re,os,sys,glob

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                if i[0]!="#":
                        Data+=[i.strip().split('\t')]
        return Data


files=glob.glob("*")
for one in files:
	files2=glob.glob(one+"/*")
        for two in files2:
        	DataL=reader(two)
                datafile=open(two+".formatted","w")
                for k in DataL:
                	datafile.write(k[0]+'\t'+k[1]+'\t'+k[2]+'\t'+k[3]+'\t'+k[5]+'\t'+k[4]+'\n')
                datafile.close()
                os.system("mv " + two+".formatted " + two)
