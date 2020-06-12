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

		# Transform the format to be recognizable by Asymmetron as bed file
                for k in DataL:
                	datafile.write(k[0]+'\t'+k[1]+'\t'+k[2]+'\t'+k[3]+'\t'+"."+'\t'+k[5]+'\t'+k[4]+'\n')
                datafile.close()
                os.system("mv " + two+".formatted " + two)
