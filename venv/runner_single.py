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

Annotation_data_LINEs=functions.strand_annotate_third_BED_overlap("/nfs/compgen-04/team218/ilias/strand_program/TCGA_ICGC/cancers/All_cancers.bed",read_BED("/nfs/compgen-04/team218/ilias/strand_program/gnomad/types/LTR_hg38.txt"))
p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(Annotation_data_LINEs,read_BED("/nfs/compgen-04/team218/ilias/strand_program/TFBSs/gencode.v33.annotation.bed"))

print([((p_p+m_m)/float(p_p+m_m+p_m+m_p))/float(template_LTRs)])
