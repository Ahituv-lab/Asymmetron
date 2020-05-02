## ASYMMETRON ###
import itertools,math
import numpy as np
from scipy.stats import binom_test
from pybedtools import BedTool
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import visualizations

# Google drive linke with the draft https://docs.google.com/document/d/1elnyyHShRcY5X406qk9O2-odU5yb7F15vWtwmSzrmrw/edit
def pairs_generator(pathL1,pathL2,NamesL1,NamesL2):
	"""
	Takes as input two sets of lists (comma separated, provided by the user
	Returns all the combinations of two between them
	The idea would be to run something like python asymmetron.py -A path1,path2,path3 -B path4,path5,path6 -Names_A Name1,Name2,Name3 -Names_B Name4,Name5,Name6
	"""
	return list(itertools.product(pathL1, pathL2)), list(itertools.product(NamesL1, NamesL2))

def read_BED(path,last_col=False):
	"""
	This function reads bed files.
	If an extra column for a score (e.g. replication timing) is given 
	Then last_Col=True and also returns another list 
	The idea of the last column is to sub-divide the  file  into groups based on the scores and do separate the analysis in each of them
	"""
	if last_col==False:
		Data = [];
		with open(path) as f:
			for line in f:
				Data.append(line.strip().split()[:5])
		return Data

	elif last_col==True:
		Data = [];Score = [];
		with open(path) as f:
			for line in f:
				Data.append(line.strip().split()[:5])
				Score.append(float(line.strip().split()[-1]))
		return Data,Score
	else:
		print("ERROR")

def readVCFtoBED(path):
	"""
	This function converts a VCF file to BED
	Not important, just for compatibility with multiple input types
	"""
	DataBED=[];
	with open(path)as f:
		for line in f:
			if line[0] == "#":
				var = line.strip().split()
				chrom="chr"+var[0].replace("chr","").replace("Chr","").replace("CHR","")
				pos = int(var[1]) # VCF is 0-based coordinate system, BED is 1-based coordinate system
				ID = var[2]
				DataBED.append([chrom,pos-1,pos,ID])
	return DataBED


def binner(min_size,max_size,bin_no):
        """
        Takes as input the distance range and number of bins and returns the bins in form (min,max) for every bin.
        """
	bin_size = float(max_size-min_size)/float(bin_no)
	Bins =[(min_size+bin_size*k,min_size+bin_size*(k+1)) for k in range(0,bin_no)]
	return Bins

def separate_on_score(path_score,path,number_of_bins,output_plot):
	"""
	Score list is ordered as DataL list of lists and the first is used to bin the second.
	This requires binning this feature and calculating the asymmetry at each bin of this column values. 
	This would be useful e.g. if we want to see if expression levels are associated with mutational strand asymmetry or if replication timing is.
	This function should divide a list of lists representing a file to list of lists (DataL) based on
	the Score bins 
	"""
	#We should check Score to be integer / float in our checks
        DataL,ScoreL=read_BED(path_score,last_col=True)
	DataL2=read_BED(path,last_col=False)
	StepsL = binner(min(ScoreL),max(ScoreL),number_of_bins)
	DataStepsL=[];ScoresStepsL=[];
	for step in StepsL:
		DataStep=[];ScoreStep=[];
		for i in range(len(ScoreL)):
			if ScoreL[i]>=step[0] and ScoreL[i]<step[1]:
				DataStep+=[DataL[i]]
				ScoreStep+=[ScoreL[i]]
		DataStepsL+=[DataStep]
		ScoresStepsL+=[ScoreStep]

	Ratio_Same_Opposite=[];
        for step in range(len(StepsL)):
		p_p_step,m_m_step,p_m_step,m_p_step,same_strand_step,opposite_strand_step,convergent_step,divergent_step=functions.overlap(DataStepsL,DataL2)
                if same_strand_step+opposite_strand_step!=[]:
                	Ratio_Same_Opposite.append(same_strand_step/float(same_strand_step+opposite_strand_step))
		else:
                	Ratio_Same_Opposite.append(0.5)
	visualizations.barplot_single_gen(Ratio_Same_Opposite,Score_names,output_plot)
	return

def asymmetries_single(path,name,window_min,window_max,patterns,bins,plot,threshold,output):
	"""
	This function calculates the strand asymmetry biases in a single file.
	Input file (If multiple files are provided the analysis is done independently in each")
	Minimum and maximum distances between consecutive instances.
	Number of bins to divide the signal in (optional). Default, no binning.
	"""
        DataL = BedTool(path).sort()
	Data_df=DataL.to_dataframe()
        Strand = list(Data_df.iloc[:,-1])
	total_plus = Strand.count("+")
        total_minus = Strand.count("-")
	probability={}
        probability["+"]=total_plus/float(total_plus+total_minus)
        probability["-"] = 1-probability["+"]

	for pattern in patterns:

        	probability_pattern = np.prod([probability[k] for k in list(pattern)])
		
		number_of_tests = len(DataL)-len(pattern)

		# Bonferoni correction
        	consecutive_threshold= next(x for x, val in enumerate(range(len(DataL))) if (probability_pattern**x)*number_of_tests > threshold) 
                Counter_consecutive_real,Counter_consecutive_control,DistancesL=calc_consecutive(list(DataL),window_min,window_max,pattern,consecutive_threshold,output)

                if plot==True:
                    consecutive, times_found = zip(*Counter_consecutive_real.items())
		    consecutive_sorted, times_found_sorted = [list(x) for x in zip(*sorted(zip(consecutive, times_found), key=lambda pair: pair[0]))]
                    indexes = np.arange(len(consecutive_sorted))
                    visualizations.barplot_single_gen(List1,List1_names,output)

		if bins>1:
			Bins=binner(window_min,window_max,bins)
			OccsL=[];
			for min_bin,max_bin in Bins:
				Occs_per_bin=0
				for dist in DistancesL:
					if dist>=min_bin and dist<max_bin:
						Occs_per_bin+=1
				OccsL.append(Occs_per_bin)
			# Plot barplot of occs consecutive in each bin
			visualizations.barplot_single_gen(List1,List1_names,output)
	return 

def extract_pattern(DataL,signS,pattern,threshold,is_real):
	"""
	Takes as input the file data, a string of all consecutive signs or "_" if distance between consecutive not within constraints.
	Takes as input the p-value threshold after which significat consecutive biases are observed. 
	Returns the file with the biased coordinates for that pattern.
	Returns consecutive occurrences list.
	"""
	occs=[];CoordinatesL=[];
	counts=0;Coordinates=[];
	for i in range(0,len(signS)-len(pattern)):
		if signS[i:i+len(pattern)]==pattern:	
			for m in DataL[i:i+len(pattern)]:
				Coordinates.append(m)
			counts+=1
		else:
			if counts>0:
				occs+=[counts]	
				if counts>=threshold:
					CoordinatesL.append(Coordinates)
				counts=0;Coordinates=[];
	if is_real:
		datafile=open("Consecutives_threshold_coordinates","w") # We need to somehow incorporate the naming here
		for line in CoordinatesL:
			datafile.write('\t'.join([str(x) for x in line])+'\n')
		datafile.close()

	return occs

def calc_consecutive(DataL,window_min,window_max,pattern,threshold_consecutiveN,output):
        """
        This function takes as input the list of lists of the fle data, the range of distances to search, the p-value threshold and aa list of strand sign patterns.
        It returns the number of consecutive occurrences in the plus and minus orientation.
        We should also include alternating occurrences
        """
	from random import shuffle #pottentially do Monte carlo instead
	from collections import Counter	
	datafile=open(output,"w")
	Signs='';
	Distances=[];
	for i in range(len(DataL)-1):
		chrom_up,start_up,end_up,name1,strand_up=DataL[i][0:5]
		chrom_down,start_down,end_down,name2,strand_down=DataL[i+1][0:5]
		distance = abs(int(start_down)-int(end_up))
		if distance>=window_min and distance<window_max:
			Distances.append(distance)
			Signs+=strand_up
		else:
			Signs+="_"
			Distances.append("_")
	datafile.close()

	# Random control
	Signs_control = list(Signs)
	shuffle(Signs_control)
	Signs_control = ''.join(Signs_control)

	Occs_pattern=extract_pattern(DataL,Signs,pattern,threshold_consecutiveN,True)
	Occs_pattern_control=extract_pattern(DataL,Signs_control,pattern,threshold_consecutiveN,False)
    	return Counter(Occs_pattern),Counter(Occs_pattern_control),Distances

def strand_annotate_third_BED_overlap(unnotated_path,annotated_path):
	"""
	For a third file that doesn't have its own annotation e.g. mutation files since mutations are on both strands
	This function enables the strand annotation of such a file
	Using an annotated file as mirror, for overlapping instances
	Returns the unnotated file, with strand annotation
	"""
	DataL_unnotated = BedTool(unnotated_path)
	DataL_annotated = BedTool(annotated_path)
	Overlap_strand = DataL_unnotated.intersect(DataL_annotated,wao=True)
	Overlap_strand_df= Overlap_strand.to_dataframe()
	Chromosome = list(Overlap_strand_df.iloc[:,0])
	Start = list(Overlap_strand_df.iloc[:,1])
	End = list(Overlap_strand_df.iloc[:,2])
	ID = list(Overlap_strand_df.iloc[:,3])
	Strand = list(Overlap_strand_df.iloc[:,-2])
	Chromosome,Start,End,ID,Strand = zip(*((chrom, start, end,id_used,strand) for chrom, start, end, id_used, strand in zip(Chromosome, Start, End, ID, Strand) if strand in ["+","-"]))
	DataL=[];
	for i in range(len(Chromosome)):
		DataL.append([Chromosome[i],Start[i],End[i],ID[i],Strand[i]])
	return DataL

def overlap(path1,path2):
	"""
	Uses pybedtools intersect function to find overlapping coordinates
	"""
	DataL1=BedTool(path1).sort()
	DataL2=BedTool(path2).sort()
	overlap=DataL1.intersect(DataL2,wao=True)
	Overlap_df=overlap.to_dataframe()
	Strand1 = list(Overlap_df.iloc[:,4])
	Strand2 = list(Overlap_df.iloc[:,9])
	p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent=orientation(Strand1,Strand2)
	return p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent

def proximal(path1,path2,name1,name2,window_min,window_max,upstream=False,downstream=False,bins=False):
	"""
       This is the main function of pairwise_asymmetries.py
       Uses pybedtools closest function to find proximal coordinates
       Then calculates asymmetry through orientation function for proximal pairs
       # the flags it uses from here https://bedtools.readthedocs.io/en/latest/content/tools/closest.html
       if bins==True then return not the counts but the lists of counts to bin them
       """
	DataL1 = BedTool(path1).sort()
	DataL2 = BedTool(path2).sort()
	if upstream==downstream and upstream==True:
		closest = DataL1.closest(DataL2,D='ref')
	elif upstream==True:
		closest = DataL1.closest(DataL2,D='a',id=False,iu=True)
	elif downstream==True:
		closest = DataL1.closest(DataL2,D='b',iu=False,id=True)
	else:
		closest = DataL1.closest(DataL2,D='ref')

	closest_df = closest.to_dataframe()
	Strand1 = list(closest_df.iloc[:,4])
	Strand2 = list(closest_df.iloc[:,9])
	Distance = [abs(i) for i in list(closest_df.iloc[:,-1])]
	Distance,Strand1,Strand2 = zip(*((dist, strand1, strand2) for dist, strand1, strand2 in zip(Distance, Strand1,Strand2) if dist <window_max and dist>=window_min))

	if bins!=False:
		p_pL_bin=[];m_mL_bin=[];p_mL_bin=[];m_pL_bin=[];same_strandL_bin=[];opposite_strandL_bin=[];convergentL_bin=[];divergentL_bin=[];
		Bins=binner(window_min,window_max,bins)
		for index, bin_i in enumerate(Bins):
			Strand1Bin=[];Strand2Bin=[];
			min_bin,max_bin = bin_i
			for k in range(len(Distance)):
				if Distance[k]>=min_bin and Distance[k]<max_bin:
					Strand1Bin.append(Strand1[k])
					Strand2Bin.append(Strand2[k])

			p_p_bin,m_m_bin,p_m_bin,m_p_bin,same_strand_bin,opposite_strand_bin,convergent_bin,divergent_bin=orientation(Strand1Bin,Strand2Bin)
			p_pL_bin.append(p_p_bin);m_mL_bin.append(m_m_bin);p_mL_bin.append(p_m_bin);m_pL_bin.append(m_p_bin);same_strandL_bin.append(same_strand_bin);opposite_strandL_bin.append(opposite_strand_bin);convergentL_bin.append(convergent_bin);divergentL_bin.append(divergent_bin)
			print p_p_bin,m_m_bin,p_m_bin,m_p_bin,bin_i

		# Same Opposite orientation
		visualizations.barplot_pair_lists_gen(Bins,same_strandL_bin,opposite_strandL_bin,name1,name2,"same_opposite_bins_"+name1+"_"+name2+".png")
		# Convergent Divergent orientation
		visualizations.barplot_pair_lists_gen(Bins,convergentL_bin,divergentL_bin,name1,name2,"convergent_divergent_bins_"+name1+"_"+name2+".png")

	return p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent

def asym_binned(window_min,window_max,bins,DistancesL,Strand1L,Strand2L):
	"""
	This function should take as input:
	window:Maximum distance to look into
	Bins:Number of bins to divide the signal into
	Distances already calculated
	Strands for factors in those distances
	This function should return:
	The strand asymmetries for each bin
	#***There is a bug with an extra bin being generated than those expected***#
	"""
	Range_Bins=binner(window_min,window_max,bins)
	
	Strand1D={i:[] for i in range(min(Range_Bins),max(Range_Bins)+1)}
	Strand2D={i:[] for i in range(min(Range_Bins),max(Range_Bins)+1)}
	DistancesD={i:[] for i in range(min(Range_Bins),max(Range_Bins)+1)}
	for i in range(len(DistancesL)):
		for index,k in enumerate(Range_Bins):
			if k[0] <= DistancesL[i] <k[1]:
				bin_of_distance = index+1
				Strand1D[index+1]+=[Strand1L[i]]
				Strand2D[index+1]+=[Strand2L[i]]
				DistancesD[index+1]+=[DistancesL[i]]

	p_pL=[];m_mL=[];p_mL=[];m_pL=[];same_strandL=[];opposite_strandL=[];convergentL=[];divergentL=[]
	for i in range(1,len(Range_Bins)+1):
		p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent=orientation(Strand1D[i],Strand2D[i])
		p_pL.append(p_p)
		m_mL.append(m_m)
		p_mL.append(p_m)
		m_pL.append(m_p)
		same_strandL.append(same_strand)
		opposite_strandL.append(opposite_strand)
		convergentL.append(convergent)
		divergentL.append(divergent)
	return p_pL,m_mL,p_mL,m_pL,same_strandL,opposite_strandL,convergentL,convergentL

def orientation(sign1L,sign2L):
	"""
	Calculates the orientation combinations
	Should only accept +/- signs
	"""
	p_p=0;m_m=0;p_m=0;m_p=0;
	for index in range(len(sign1L)):
		sign1=sign1L[index]
		sign2=sign2L[index]
		if sign1 in ["+","-"] and sign2 in ["+","-"]:
			if sign1==sign2:
				if sign1=="+":
					p_p+=1
				elif sign1=="-":
					m_m+=1
			else:
				if sign1=="+" and sign2=="-":
					p_m+=1
				elif sign1=="-" and sign2=="+":
					m_p+=1
	same_strand=p_p+m_m;
	opposite_strand=p_m+m_p;
	convergent=p_m;
	divergent=m_p;
	return p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent

def ratio_calc(first_strandI,second_strandI):
	"""
	Calculates ratio of strand asymmetry, accepts integer numbers in two strands
	"""
	if first_strandI+second_strandI!=0:
		Ratio = first_strandI/float(first_strandI+second_strandI)
		return Ratio
	else:
		return np.nan

def statistical_evaluation(occs_strand1,occs_strand2,number_of_files_scanned,expected_asym=0.5):
	"""
	Inputs: Occurrences in each strand and orientation, integers
	Calculates the statistical significance
	Binomial test: p-value
	Bonferonni correction: p-value	
	Ratio of strand asymmetry
	if the user expects a background bias then that should be given as an expected ratio, different than zero (that would alter our bionimal test, and we could just  change 0.5 to a variable)
	"""
	Ratio_strand1_2 = ratio_calc(occs_strand1,occs_strand2)
	p_val_strand1_2 = binom_test(occs_strand1,occs_strand1+occs_strand2,expected_asym)
	p_val_strand1_2_Bonferoni = min(1,p_val_strand1_2*number_of_files_scanned)
	return Ratio_strand1_2,p_val_strand1_2,p_val_strand1_2_Bonferoni


def table_gen(NamesL_pairs,p_pL,m_mL,p_mL,m_pL,p_valsL,p_vals_BonferoniL,RatiosL,p_valsL_divergent_convergent,p_valsL_divergent_convergent_BonferoniL,RatiosL_divergent_convergent):
	"""
	accepts list of  pairs of files / Names (factors) as input, list of p-values, Ratios and returns table
	This is more useful for multiple lists of files.
	"""
	datafile=open("statistics_asymmetry.txt","w")
	datafile.write("Feature 1"+'\t'+"Feature 2"+"\t"+"plus-plus"+'\t'+"minus-minus"+'\t'+"plus-minus"+'\t'+"minus-plus"+'\t'+"p-value same opposite"+'\t'+"p-value same opposite Bonferoni corrected"+'\t'+"Ratio same opposite"+'\t'+"p-value divergent convergent"+'\t'+"p-value divergent convergent Bonferoni corrected"+'\t'+"Ratio divergent convergent"+'\n')
	for i in range(len(NamesL)):
		datafile.write(NamesL_pairs[i][0]+'\t'+NamesL_pairs[i][1]+'\t'+str(p_pL[i])+'\t'+str(m_mL[i])+'\t'+str(p_mL[i])+'\t'+str(m_pL[i])+'\t'+str(p_valsL[i])+'\t'+str(p_vals_BonferoniL[i])+'\t'+str(RatiosL[i])+'\t'+str(p_valsL_divergent_convergent[i])+'\t'+str(p_valsL_divergent_convergent_BonferoniL[i])+'\t'+str(RatiosL_divergent_converg[i])+'\n')
	datafile.close()
	return

# Ensures that code below is not run when this file is imported into another file
if __name__ == "__main__":
	pass
# test area
#works
#print overlap(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"))
#works
#print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,False,False,False)
#works
#print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,True,False,False)
#works
#print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,False,True,False)
#works
#print strand_annotate_third_BED_overlap(read_BED("Myeloid.indels"),read_BED("Ensembl.genes_hg19_TSSs.bed"))
#works
#DataL,ScoreL=read_BED("MCF7_RepliStrand.leading_lagging.bed",True)
#separate_on_score(ScoreL,DataL,10)
# works - minor error with extra bin, needs fixing
#Strand1,Strand2,DistancesL=proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),"G4","hg19_TSS",0,500,False,False,10)
#print len(Strand1),len(Strand2),len(DistancesL)
#print asym_binned(0,500,10,DistancesL,Strand1,Strand2)
#asymmetries_single(read_BED("All_G4.bed"),0,1000,["++--"],0,False,0.0001,"test")
