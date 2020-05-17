import re,os,sys,glob
import functions
import operator
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from operator import itemgetter

def read_BED(path,boot=False):
        """ 
        This function reads bed files.
        If an extra column for a score (e.g. replication timing) is given 
        Then last_Col=True and also returns another list 
        The idea of the last column is to sub-divide the  file  into groups based on the scores and do separate the analysis in each of them
        """
        if boot==False:
                Data = []
                with open(path) as f:
                        for line in f:
                                Data.append(line.strip().split()[:5])
                return Data
    
        elif boot==True:
                Data = []    
                with open(path) as f:
                        for line in f:
                                Data.append(line.strip().split()[:5])
                Data_Boot=[]
                for sim in range(len(Data)):
                        Data_Boot+=[random.choice(Data)]
                return Data_Boot
        else:
                print("ERROR")


cancers = ["BLCA-US","BRCA-US","CMDI-UK","ESAD-UK","KICH-US","LGG-US","LIRI-JP","MELA-AU","PACA-AU","PBCA-DE","READ-US","STAD-US","BOCA-UK","BTCA-SG","COAD-US","GACA-CN","KIRC-US","LICA-FR","LUAD-US","ORCA-IN","PACA-CA","PRAD-CA","RECA-EU","THCA-US","BRCA-EU","CESC-US","DLBC-US","GBM-US","KIRP-US","LIHC-US","LUSC-US","OV-AU","PAEN-AU","PRAD-UK","SARC-US","UCEC-US","BRCA-UK","CLLE-ES","EOPC-DE","HNSC-US","LAML-KR","LINC-JP","MALY-DE","OV-US","PAEN-IT","PRAD-US","SKCM-US"]

tissues={"BRCA-EU":"Breast_Cancer","BRCA-UK":"Breast_Cancer","BRCA-US":"Breast_Cancer","Biliary":"BTCA-SG","Biliary":"LIRI-JP","Bladder":"BLCA-US","Bone/SoftTissue":"BOCA-UK","Bone/SoftTissue":"SARC-US","Bone/SoftTissue":"BOCA-UK","Bone/SoftTissue":"SARC-US","Cervix":"CESC-US","CNS":"PBCA-DE","CNS":"GBM-US","CNS":"LGG-US","COAD-US":"Colon/Rectum","READ-US":"Colon/Rectum","ESAD-UK":"Esophagus","HNSC-US":"Head/Neck","ORCA-IN":"Head/Neck","KIRC-US":"Kidney","RECA-EU":"Kidney","KICH-US":"Kidney","LICA-FR":"Liver","LIHC-US":"Liver","LINC-JP":"Liver","LIRI-JP":"Liver","LUSC-US":"Lung","LUAD-US":"Lung","DLBC-US":"Lymphoid","MALY-DE":"Lymphoid","CLL-ES":"Lymphoid","CLLE-ES":"Lymphoid","CMDI-UK":"Myeloid","LAML-KR":"Myeloid","LAML-US":"Myeloid","OV-AU":"Ovary","OV-US":"Ovary","PACA-AU":"Pancreas","PACA-CA":"Pancreas","PAEN-AU":"Pancreas","EOPC-DE":"Prostate","PRAD-CA":"Prostate","PRAD-UK":"Prostate","PRAD-US":"Prostate","MELA-AU":"Skin","SKCM-US":"Skin","STAD-US":"Stomach","THCA-US":"Thyroid","UCEC-US":"Uterus"}


tissuesD={"Breast_Cancer":["BRCA-EU","BRCA-UK","BRCA-US"],"Biliary":["BTCA-SG","LIRI-JP"],"Bladder":["BLCA-US"],"Bone/SoftTissue":["BOCA-UK","SARC-US","BOCA-UK","SARC-US"],"Cervix":["CESC-US"],"CNS":["PBCA-DE","GBM-US","LGG-US"],"Colon/Rectum":["COAD-US","READ-US"],"Esophagus":["ESAD-UK"],"Head/Neck":["HNSC-US","ORCA-IN"],"Kidney":["KIRC-US","RECA-EU","KICH-US"],"Liver":["LICA-FR","LIHC-US","LINC-JP","LIRI-JP"],"Lung":["LUSC-US","LUAD-US"],"Lymphoid":["DLBC-US","MALY-DE","CLL-ES","CLLE-ES"],"Myeloid":["CMDI-UK","LAML-KR","LAML-US"],"Ovary":["OV-AU","OV-US"],"Pancreas":["PACA-AU","PACA-CA","PAEN-AU"],"Prostate":["PRAD-CA","PRAD-UK","PRAD-US"],"Skin":["MELA-AU","SKCM-US"],"Stomach":["STAD-US"],"Thyroid":["THCA-US"],"Uterus":["UCEC-US"]}


p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/LINE_hg38.txt"),read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))
template_LINEs=((p_p+m_m)/float(p_p+m_m+p_m+m_p))

p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/SINE_hg38.txt"),read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))
template_SINEs=((p_p+m_m)/float(p_p+m_m+p_m+m_p))

p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/LTR_hg38.txt"),read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))
template_LTRs=((p_p+m_m)/float(p_p+m_m+p_m+m_p))

print(template_LINEs,template_SINEs,template_LTRs,"LINEs","SINEs","LTRs")
Scores_LINEs=[];Scores_SINEs=[];Scores_LTRs=[];
Scores_LINEsLL=[];Scores_SINEsLL=[];Scores_LTRsLL=[];
NamesL=[]
for cancer_tissue in tissuesD.keys()[:10]:
	p_pL_LINE=0; m_mL_LINE=0; p_mL_LINE=0; m_pL_LINE=0;
	p_pL_SINE=0; m_mL_SINE=0; p_mL_SINE=0; m_pL_SINE=0;
	p_pL_LTRs=0; m_mL_LTRs=0; p_mL_LTRs=0; m_pL_LTRs=0;
	for cancer in tissuesD[cancer_tissue]:

		files=glob.glob("/nfs/compgen-04/team218/ilias/strand_program/TCGA_ICGC/cancers/"+cancer+"/*formatted.bed")
		perLine_LINE=[]
		perLine_SINE=[]
		perLine_LTRs=[]
		for one in files:
			print(Scores_LINEsLL)
			# Annotate LINEs
			print(one,cancer)
			try:
				Annotation_data_LINEs=functions.strand_annotate_third_BED_overlap(one,read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/LINE_hg38.txt"))
				p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(Annotation_data_LINEs,read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))

				p_pL_LINE+=p_p
				m_mL_LINE+=m_m
				p_mL_LINE+=p_m
				m_pL_LINE+=m_p
        	                if p_p+m_m+p_m+m_p!=0:
                        	        perLine_LINE+=[((p_p+m_m)/float(p_p+m_m+p_m+m_p))/float(template_LINEs)]
				print(p_p+m_m+p_m+m_p)
			except:
				pass


			try:
                	        Annotation_data_SINEs=functions.strand_annotate_third_BED_overlap(one,read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/SINE_hg38.txt"))
                	        p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(Annotation_data_SINEs,read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))
                	        p_pL_SINE+=p_p
                	        m_mL_SINE+=m_m
                	        p_mL_SINE+=p_m
                	        m_pL_SINE+=m_p
                	        if p_p+m_m+p_m+m_p!=0:
                	                perLine_SINE+=[((p_p+m_m)/float(p_p+m_m+p_m+m_p))/float(template_SINEs)]
				print(p_p+m_m+p_m+m_p)
                	except:
                	        pass



                	try:
                	        Annotation_data_LTRs=functions.strand_annotate_third_BED_overlap(one,read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/LTR_hg38.txt"))
                	        p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(Annotation_data_LTRs,read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))
                	        p_pL_LTRs+=p_p
                        	m_mL_LTRs+=m_m
                        	p_mL_LTRs+=p_m
                        	m_pL_LTRs+=m_p
                        	if p_p+m_m+p_m+m_p!=0:
                        	        perLine_LTRs+=[((p_p+m_m)/float(p_p+m_m+p_m+m_p))/float(template_LTRs)]
                	except:	
                	        pass

	if p_pL_LINE+m_mL_LINE+p_mL_LINE+m_pL_LINE!=0:
		Ratio = ((p_pL_LINE+m_mL_LINE)/float(p_pL_LINE+m_mL_LINE+p_mL_LINE+m_pL_LINE))
		Scores_LINEs.append(Ratio/float(template_LINEs))
		Scores_LINEsLL.append(perLine_LINE)
	#else:
	#	Scores_LINEsLL.append(1)
	#	Scores_LINEs.append(1)

	if p_pL_SINE+m_mL_SINE+p_mL_SINE+m_pL_SINE!=0:
        	Ratio = ((p_pL_SINE+m_mL_SINE)/float(p_pL_SINE+m_mL_SINE+p_mL_SINE+m_pL_SINE))
        	Scores_SINEs.append(Ratio/float(template_SINEs))
		Scores_SINEsLL.append(perLine_SINE)
	#else:
	#	Scores_SINEs.append(1)
	#	Scores_SINEsLL.append(1)

	if p_pL_LTRs+m_mL_LTRs+p_mL_LTRs+m_pL_LTRs!=0:
        	Ratio = ((p_pL_LTRs+m_mL_LTRs)/float(p_pL_LTRs+m_mL_LTRs+p_mL_LTRs+m_pL_LTRs))
        	Scores_LTRs.append(Ratio/float(template_LTRs))
		Scores_LTRsLL.append(perLine_LTRs)
	#else:
	#	Scores_LTRs.append(1)
	#	Scores_LTRsLL.append(1)

	NamesL.append(cancer_tissue)
print(NamesL)

import matplotlib.patches as mpatches
labels = []
def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))


Scores_LINEsLL_medians=[np.mean(k) for k in Scores_LINEsLL]
Scores_LINEsLL_medians_,NamesL=[list(x) for x in zip(*sorted(zip(Scores_LINEsLL_medians, NamesL), key=itemgetter(0)))]
Scores_LINEsLL_medians_,Scores_LINEsLL=[list(x) for x in zip(*sorted(zip(Scores_LINEsLL_medians, Scores_LINEsLL), key=itemgetter(0)))]
Scores_LINEsLL_medians_,Scores_SINEsLL=[list(x) for x in zip(*sorted(zip(Scores_LINEsLL_medians, Scores_SINEsLL), key=itemgetter(0)))]
Scores_LINEsLL_medians_,Scores_LTRsLL=[list(x) for x in zip(*sorted(zip(Scores_LINEsLL_medians, Scores_LTRsLL), key=itemgetter(0)))]

ax = plt.subplot(111)
plt.boxplot(Scores_LINEsLL,positions=range(1,len(NamesL)*4+1,4))#,"LINEs")
plt.boxplot(Scores_SINEsLL,positions=range(2,len(NamesL)*4+1,4))#,"SINEs")
plt.boxplot(Scores_LTRsLL,positions=range(3,len(NamesL)*4+1,4))#,"LTRs")
plt.xticks(range(2,len(NamesL)*4+1,4),NamesL,rotation=90)
plt.tight_layout()
plt.xlim(xmin=0)
plt.ylim(ymin=0)
#plt.ylabel("log(Strand Asymmetry)")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.ylabel("Strand Asymmetry")
plt.savefig("Cancers_Strand_Asymmetries_LINE_SINE_LTR.png")
plt.close()
