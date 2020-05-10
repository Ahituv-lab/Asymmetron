## ASYMMETRON ###
import itertools, math
import numpy as np
from scipy.stats import binom_test
try:
    from pybedtools import BedTool
except:
    print("Pybedtools not imported")
try:
    import visualizations
except:
    print("visualisations not imported")


# Google drive linke with the draft https://docs.google.com/document/d/1elnyyHShRcY5X406qk9O2-odU5yb7F15vWtwmSzrmrw/edit
def pairs_generator(pathL1, pathL2, NamesL1, NamesL2):
    """
    Takes as input two sets of lists (comma separated, provided by the user
    Returns all the combinations of two between them
    The idea would be to run something like python asymmetron.py -A path1,path2,path3 -B path4,path5,path6 -Names_A Name1,Name2,Name3 -Names_B Name4,Name5,Name6
    """
    return list(itertools.product(pathL1, pathL2)), list(itertools.product(NamesL1, NamesL2))


def read_BED(path, last_col=False):
    """
    This function reads bed files.
    last_col=True: If an extra column for a score (e.g. replication timing or gene expression) is given
    This function returns a list of lists with the BED file information and the list of scores.
    """
    if last_col == False:
        Data = [];
        with open(path) as f:
            for line in f:
                Data.append(line.strip().split()[:5])
        return Data

    elif last_col == True:
        Data = [];
        Score = [];
        with open(path) as f:
            for line in f:
                Data.append(line.strip().split()[:5])
                Score.append(float(line.strip().split()[-1]))
        return Data, Score
    else:
        print("ERROR")


def readVCFtoBED(path):
    """
    This function converts a VCF file to BED
    Not important, just for compatibility with multiple input types
    """
    DataBED = [];
    with open(path)as f:
        for line in f:
            if line[0] == "#":
                var = line.strip().split()
                chrom = "chr" + var[0].replace("chr", "").replace("Chr", "").replace("CHR", "")
                pos = int(var[1])  # VCF is 0-based coordinate system, BED is 1-based coordinate system
                ID = var[2]
                DataBED.append([chrom, pos - 1, pos, ID])
    return DataBED


def binner(min_size, max_size, bin_no):
    """
    Takes as input the distance range and number of bins and returns the bins in form (min,max) for every bin.
    """
    bin_size = float(max_size - min_size) / float(bin_no)
    Bins = [(min_size + bin_size * k, min_size + bin_size * (k + 1)) for k in range(0, bin_no)]
    return Bins


def separate_on_score(path_score, path, number_of_bins):
    """
    path_score: is the path to the BED file with last column containing the scores.
    path: the path to the second BED file.
    number_of_bins: number of groups to generate based on the score min and max values and in which to perform the strand asymmetry analysis.
    This would be useful e.g. if we want to see if expression levels are associated with mutational strand asymmetry or if replication timing is.
    returns the strand asymmetry scores for same / opposite and convergent / divergent for each bin as two lists.
    """
    # We should check Score to be integer / float in our checks
    DataL, ScoreL = read_BED(path_score, last_col=True)
    DataL2 = read_BED(path, last_col=False)
    StepsL = binner(min(ScoreL), max(ScoreL), number_of_bins)

    # Separates the DataL based on the ScoreL bins intro groups.
    DataStepsL = [];ScoresStepsL = [];
    for step in StepsL:
        DataStep = [];ScoreStep = [];
        for i in range(len(ScoreL)):
            if ScoreL[i] >= step[0] and ScoreL[i] <= step[1]:
                DataStep += [DataL[i]]
                ScoreStep += [ScoreL[i]]
        DataStepsL += [DataStep]
        ScoresStepsL += [ScoreStep]

    # Calculates the asymmetry for same / opposite and convergent / divergent asymmetry for each bin.
    Ratio_Same_Opposite = [];Ratio_Convergent_Divergent=[];
    for step in range(len(StepsL)):
        p_p_step, m_m_step, p_m_step, m_p_step, same_strand_step, opposite_strand_step, convergent_step, divergent_step = overlap(DataStepsL[step], DataL2)
        if same_strand_step + opposite_strand_step != 0:
            Ratio_Same_Opposite.append(same_strand_step / float(same_strand_step + opposite_strand_step))
        else:
            Ratio_Same_Opposite.append(0.5)

        if p_m_step+m_p_step != 0:
            Ratio_Convergent_Divergent.append(p_m_step/ float(p_m_step+ m_p_step))

    return Ratio_Same_Opposite,Ratio_Convergent_Divergent,StepsL


def strand_annotate_third_BED_overlap(unnotated_path, annotated_path):
    """
    For a third file that doesn't have its own annotation e.g. mutation files since mutations are on both strands
    This function enables the strand annotation of such a file
    Using an annotated file as mirror, for overlapping instances
    Returns the unnotated file, with strand annotation
    """
    DataL_unnotated = BedTool(unnotated_path)
    DataL_annotated = BedTool(annotated_path)
    Overlap_strand = DataL_unnotated.intersect(DataL_annotated, wao=True)
    Overlap_strand_df = Overlap_strand.to_dataframe()
    Chromosome = list(Overlap_strand_df.iloc[:, 0])
    Start = list(Overlap_strand_df.iloc[:, 1])
    End = list(Overlap_strand_df.iloc[:, 2])
    ID = list(Overlap_strand_df.iloc[:, 3])
    Strand = list(Overlap_strand_df.iloc[:, -2])
    Chromosome, Start, End, ID, Strand = zip(
        *((chrom, start, end, id_used, strand) for chrom, start, end, id_used, strand in
          zip(Chromosome, Start, End, ID, Strand) if strand in ["+", "-"]))
    DataL = [];
    for i in range(len(Chromosome)):
        DataL.append([Chromosome[i], Start[i], End[i], ID[i], Strand[i]])
    return DataL


def overlap(path1, path2):
    """
    This is the main function of contained_asymmetries.py
    Takes as inputs the paths to the two BED files to compare.
    Uses pybedtools intersect function to find overlapping coordinates between regions and motifs.
    Returns the number of occurrences in each orientation for Regions:BED path1 to Motifs:BED path2.
    """
    DataL1 = BedTool(path1).sort()
    DataL2 = BedTool(path2).sort()
    overlap = DataL1.intersect(DataL2, wao=True)
    Overlap_df = overlap.to_dataframe()
    Strand1 = list(Overlap_df.iloc[:, 4])
    Strand2 = list(Overlap_df.iloc[:, 9])
    p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = orientation(Strand1, Strand2)
    return p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent


def proximal(path1, path2, name1, name2, window_min, window_max, upstream=False, downstream=False, bins=None):
    """
       This is the main function of pairwise_asymmetries.py
       Uses pybedtools closest function to find proximal coordinates
       Then calculates asymmetry through orientation function for proximal pairs
       # the flags it uses from here https://bedtools.readthedocs.io/en/latest/content/tools/closest.html
       if bins==True then return not the counts but the lists of counts to bin them
       """
    # Finds the occurrences within the proximity limits and saves their pairwise orientation.
    DataL1 = BedTool(path1).sort()
    DataL2 = BedTool(path2).sort()
    if upstream == downstream and upstream == True:
        closest = DataL1.closest(DataL2, D='ref')
    elif upstream == True:
        closest = DataL1.closest(DataL2, D='a', id=False, iu=True)
    elif downstream == True:
        closest = DataL1.closest(DataL2, D='b', iu=False, id=True)
    else:
        closest = DataL1.closest(DataL2, D='ref')

    closest_df = closest.to_dataframe()
    Strand1 = list(closest_df.iloc[:, 4])
    Strand2 = list(closest_df.iloc[:, 9])
    Distance = [abs(i) for i in list(closest_df.iloc[:, -1])]
    Distance, Strand1, Strand2 = zip(
        *((dist, strand1, strand2) for dist, strand1, strand2 in zip(Distance, Strand1, Strand2) if
          dist < window_max and dist >= window_min))
    p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = orientation(Strand1, Strand2)

    # Calculate the distance distributions for all orientations
    Distances_orientations = get_distance_orientations(Distance,Strand1,Strand2,window_min,window_max)

    p_pL_bin = [];m_mL_bin = []; # Same orientation
    p_mL_bin = [];m_pL_bin = []; # Opposite orientation
    same_strandL_bin = []; opposite_strandL_bin = []; # Combined same / opposite orientations
    convergentL_bin = [];divergentL_bin = [];

    Bins = [];
    if bins != None:
        # Performs the same analysis for each bin.
        Bins = binner(window_min, window_max, bins)
        for index, bin_i in enumerate(Bins):
            Strand1Bin = [];Strand2Bin = [];
            min_bin, max_bin = bin_i
            for k in range(len(Distance)):
                if Distance[k] >= min_bin and Distance[k] < max_bin:
                    Strand1Bin.append(Strand1[k])
                    Strand2Bin.append(Strand2[k])

            p_p_bin, m_m_bin, p_m_bin, m_p_bin, same_strand_bin, opposite_strand_bin, convergent_bin, divergent_bin = orientation(
                Strand1Bin, Strand2Bin)
            p_pL_bin.append(p_p_bin);m_mL_bin.append(m_m_bin); #Same orientation, per bin
            p_mL_bin.append(p_m_bin);m_pL_bin.append(m_p_bin); #Opposite orientation per bin
            same_strandL_bin.append(same_strand_bin);opposite_strandL_bin.append(opposite_strand_bin);
            convergentL_bin.append(convergent_bin);divergentL_bin.append(divergent_bin)

    return (Distances_orientations, p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent), (Bins,p_pL_bin,m_mL_bin,p_mL_bin,m_pL_bin,same_strandL_bin,opposite_strandL_bin,convergentL_bin,divergentL_bin)


def get_distance_orientations(DistanceL,Strand1L,Strand2L,window_min,window_max):
    same_strandL_distance=[];opposite_strandL_distance=[];
    divergentL_distance=[];convergentL_distance=[];
    for index in range(len(Strand1L)):
        if (DistanceL[index]< window_max and DistanceL[index] >= window_min):
            sign1 = Strand1L[index]
            sign2 = Strand2L[index]
            if sign1 in ["+", "-"] and sign2 in ["+", "-"]:
                if sign1 == sign2:
                     same_strandL_distance.append(DistanceL[index])
                else:
                     opposite_strandL_distance.append(DistanceL[index])

                     if sign1 == "+" and sign2 == "-":
                         convergentL_distance.append(DistanceL[index])
                     elif sign1 == "-" and sign2 == "+":
                         divergentL_distance.append(DistanceL[index])

    return (same_strandL_distance,opposite_strandL_distance,divergentL_distance,convergentL_distance)

def asym_binned(window_min, window_max, bins, DistancesL, Strand1L, Strand2L):
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
    Range_Bins = binner(window_min, window_max, bins)

    Strand1D = {i: [] for i in range(min(Range_Bins), max(Range_Bins) + 1)}
    Strand2D = {i: [] for i in range(min(Range_Bins), max(Range_Bins) + 1)}
    DistancesD = {i: [] for i in range(min(Range_Bins), max(Range_Bins) + 1)}
    for i in range(len(DistancesL)):
        for index, k in enumerate(Range_Bins):
            if k[0] <= DistancesL[i] < k[1]:
                bin_of_distance = index + 1
                Strand1D[index + 1] += [Strand1L[i]]
                Strand2D[index + 1] += [Strand2L[i]]
                DistancesD[index + 1] += [DistancesL[i]]

    p_pL = [];m_mL = []; # Same strand orientation
    p_mL = [];m_pL = []; # Opposite strand orientation
    same_strandL = [];
    opposite_strandL = [];
    convergentL = [];
    divergentL = []
    for i in range(1, len(Range_Bins) + 1):
        p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = orientation(Strand1D[i], Strand2D[i])
        p_pL.append(p_p);m_mL.append(m_m);
        p_mL.append(p_m);m_pL.append(m_p);
        same_strandL.append(same_strand)
        opposite_strandL.append(opposite_strand)
        convergentL.append(convergent)
        divergentL.append(divergent)
    return p_pL, m_mL, p_mL, m_pL, same_strandL, opposite_strandL, convergentL, convergentL


def orientation(sign1L, sign2L):
    """
    Calculates the orientation combinations
    Should only accept +/- signs
    """
    p_p = 0;
    m_m = 0;
    p_m = 0;
    m_p = 0;
    for index in range(len(sign1L)):
        sign1 = sign1L[index]
        sign2 = sign2L[index]
        if sign1 in ["+", "-"] and sign2 in ["+", "-"]:
            if sign1 == sign2:
                if sign1 == "+":
                    p_p += 1
                elif sign1 == "-":
                    m_m += 1
            else:
                if sign1 == "+" and sign2 == "-":
                    p_m += 1
                elif sign1 == "-" and sign2 == "+":
                    m_p += 1
    same_strand = p_p + m_m;
    opposite_strand = p_m + m_p;
    convergent = p_m;
    divergent = m_p;
    return p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent


def ratio_calc(first_strandI, second_strandI):
    """
    Calculates ratio of strand asymmetry, accepts integer numbers in two strands
    """
    if first_strandI + second_strandI != 0:
        Ratio = first_strandI / float(first_strandI + second_strandI)
        return Ratio
    else:
        return np.nan


def statistical_evaluation(occs_strand1, occs_strand2, number_of_files_scanned, expected_asym=0.5):
    """
    Inputs: Occurrences in each strand and orientation, integers
    Calculates the statistical significance
    Binomial test: p-value
    Bonferonni correction: p-value
    Ratio of strand asymmetry
    if the user expects a background bias then that should be given as an expected ratio, different than zero (that would alter our bionimal test, and we could just  change 0.5 to a variable)
    """
    Ratio_strand1_2 = ratio_calc(occs_strand1, occs_strand2)
    p_val_strand1_2 = binom_test(occs_strand1, occs_strand1 + occs_strand2, expected_asym)
    p_val_strand1_2_Bonferoni = min(1, p_val_strand1_2 * number_of_files_scanned)
    return Ratio_strand1_2, p_val_strand1_2, p_val_strand1_2_Bonferoni


def table_gen(NamesL_pairs, p_pL, m_mL, p_mL, m_pL, p_valsL, p_vals_BonferoniL, RatiosL, p_valsL_divergent_convergent,
              p_valsL_divergent_convergent_BonferoniL, RatiosL_divergent_convergent, output_table):
    """
    accepts list of  pairs of files / Names (factors) as input, list of p-values, Ratios and returns table
    This is more useful for multiple lists of files.
    """
    datafile = open(output_table, "w")
    datafile.write(
        "Feature 1" + '\t' + "Feature 2" + "\t" + "plus-plus" + '\t' + "minus-minus" + '\t' + "plus-minus" + '\t' + "minus-plus" + '\t' + "p-value same opposite" + '\t' + "p-value same opposite Bonferoni corrected" + '\t' + "Ratio same opposite" + '\t' + "p-value divergent convergent" + '\t' + "p-value divergent convergent Bonferoni corrected" + '\t' + "Ratio divergent convergent" + '\n')
    for i in range(len(NamesL_pairs)):
        datafile.write(
            NamesL_pairs[i][0] + '\t' + NamesL_pairs[i][1] + '\t' + str(p_pL[i]) + '\t' + str(m_mL[i]) + '\t' + str(
                p_mL[i]) + '\t' + str(m_pL[i]) + '\t' + str(p_valsL[i]) + '\t' + str(p_vals_BonferoniL[i]) + '\t' + str(
                RatiosL[i]) + '\t' + str(p_valsL_divergent_convergent[i]) + '\t' + str(
                p_valsL_divergent_convergent_BonferoniL[i]) + '\t' + str(RatiosL_divergent_convergent[i]) + '\n')
    datafile.close()
    return

def table_bins_gen(NamesL_pairs,Score_names,Ratio_Bins,Ratio_Convergent_Divergent_Bins,output_table):
    datafile = open(output_table, "w")
    datafile.write("Score Range Min"+'\t'+str("Score Range Max")+'\t'+"Feature 1"+'\t'+"Feature 2"+"\t"+"Ratio same opposite"+'\t'+"Ratio divergent convergent" + '\n')
    for i in range(len(NamesL_pairs)):
        datafile.write(str(round(Score_names[i][0],0))+"-"+str(round(Score_names[i][1],0))+'\t'+str(NamesL_pairs[i][0])+'\t'+str(NamesL_pairs[i][1])+'\t'+str(Ratio_Bins[i])+'\t'+str(Ratio_Convergent_Divergent_Bins[i])+'\n')
    datafile.close()
    return

def table_consecutive(ConsecutiveL,namesL,path_out):
    max_consecutive = max([max(k.keys()) for k in ConsecutiveL])
    with open(path_out, 'w') as output:
        output.write(str("Number of consecutive occurrences")+'\t'+'\t'.join([str(x) for x in range(1,max_consecutive+1)])+'\n')

        for i in range(len(ConsecutiveL)):
            ConsecutiveLT = sorted(ConsecutiveL[i].items()) 
            consecutive, times_found = zip(*ConsecutiveLT) 
            consecutiveL = list(consecutive)
            times_foundL = list(times_found)            

            for k in range(1,max_consecutive+1):
               if k not in consecutiveL:
                   consecutiveL=consecutiveL[:k-1]+[k]+consecutiveL[k-1:]
                   times_foundL=times_foundL[:k-1]+[0]+times_foundL[k-1:]
           
            output.write(str(namesL[i])+'\t'+'\t'.join([str(x) for x in times_foundL])+'\n')
    return



def write_BED_out(DataL,path_out):
    # Write significant results in an output file
    with open(path_out, 'w') as output_file:
        for line in DataL:
            output_file.write('\t'.join([str(x) for x in list(line)])+'\n')
    return 


def find_sub_str(my_str, sub_str):
    start = 0
    while True:
        start = my_str.find(sub_str, start)
        if start == -1: return
        yield start
        start += len(sub_str)  # use start += 1 to find overlapping matches


def extract_pattern(DataL, pattern, min_distance, max_distance, threshold):
    """
    Finds all occurences of the pattern in BED file that meet the distance requirements.
    :param DataL: A list, where each element is a list, representing a line of a sorted BED-file
    :param pattern: The pattern to search
    :param min_distance: Minimum distance between two motifs. If two occurences of the motif do not have this minimum
                         distance, they break the pattern
    :param max_distance: Maximum distance between two motifs. If two occurences of the motif do not have this minimum
                         distance, they break the pattern
    :param threshold:
    :return: A list of all occurences of the pattern meeting the minimum distance requirements
    """
    n = len(pattern)

    # Create a string with all signs (5th column in the sorted BED file)
    signs = ""
    for line in DataL:
        signs += line[4]

    # Find all occurences of the pattern in the string of signs without accounting for distances
    occs = list(find_sub_str(signs, pattern))

    # Remove occurences that do not meet the distance criterion
    for i in range(len(occs)):
        index = occs[i]
        for j in range(n - 1):
            distance = max(0, int(DataL[index + j + 1][1]) - int(DataL[index + j][2]))
            if distance < min_distance or distance >= max_distance:
                occs[i] = None
                break
    occs = [x for x in occs if x is not None]

    # Here we translate the probability threshold to consecutive occcurrences 
    number_of_tests = len(DataL)
    total_plus = signs.count("+")
    total_minus = signs.count("-")
    probability={}
    probability["+"]=total_plus/float(total_plus+total_minus)
    probability["-"] = 1-probability["+"]
    probability_pattern = np.prod([probability[k] for k in list(pattern)])
    consecutive_threshold= next(x for x, val in enumerate(range(len(DataL))) if (probability_pattern**x)*number_of_tests < threshold)

    # Filter for number of consecutive occurences that meet the threshold criterion and are within the distance window
    DataL_significant = []
    counter = 1
    DataL_temp = DataL[occs[0]:occs[0]+n];
    from collections import defaultdict
    consecutiveL= defaultdict(int);
    for i in range(1, len(occs)):
        index = occs[i]
        distance = max(0, int(DataL[index][1]) - int(DataL[index -1][2]))
        if occs[i]-occs[i-1] == n and (distance >= min_distance and distance < max_distance):
            counter += 1
            DataL_temp.extend(DataL[index:index+n])
        else:
            consecutiveL[counter]+=1
            if counter >= consecutive_threshold:
                # Could also dump to file here to save on memory
                DataL_significant.extend(DataL_temp)
            DataL_temp.clear()
            counter = 1

    # Add last lines
    if counter >= consecutive_threshold:
        DataL_significant.extend(DataL_temp)
        consecutiveL[counter] += 1

    distancesL = []
    for i in range(len(occs)-1):
        index = occs[i]
        index_next = occs[i+1]
        distance = max(0, int(DataL[occs[i+1]][1]) - int(DataL[occs[i]+n-1][2]))
        distancesL.append(distance)

    return consecutiveL,distancesL,DataL_significant

def asymmetries_single(path, patternsL, min_distance, max_distance, threshold):
    from random import shuffle
    DataL = read_BED(path,False)

    # Here we want to shuffle the strand column
    ColL = [k[4] for k in DataL]
    shuffle(ColL)
    DataL_random=[]
    for i in range(len(DataL)):
        line = DataL[i][:-1] + [ColL[i]]
        DataL_random.append(line)

    consecutiveL=[];occsL=[];DataL_significantL=[];consecutive_controlL=[];occs_controlL=[];DataL_significant_controlL=[];
    for pattern in patternsL:
        consecutive,occs, DataL_significant = extract_pattern(DataL, pattern, min_distance, max_distance, threshold)
        consecutive_control,occs_control, DataL_significant_control = extract_pattern(DataL_random, pattern, min_distance, max_distance, threshold)
        consecutiveL.append(consecutive);occsL.append(occs);consecutive_controlL.append(consecutive_control);occs_controlL.append(occs_control);
        DataL_significantL.append(DataL_significant);
    return DataL_significantL,consecutiveL,occsL,consecutive_controlL,occs_controlL

# Ensures that code below is not run when this file is imported into another file
if __name__ == "__main__":
    with open("test_extract_pattern.bed", "r") as f:
        DataL = []
        for line in f.readlines():
            DataL.append(line.strip().split("\t"))
    out = extract_pattern(DataL, "+-", 0, 3, 2)
    print("The following dictionary includes the number of consecutive appearances of the pattern, e.g. when looking "
          "for +- in +-+-+---+- the result should be {1:1}, {3:1}", out[0])
    print("The distances between consecutive appearances of the pattern are: ", out[1] )
    print("The following lines are part of a sequence of consecutive repetitions of the pattern that meet both the "
          "threshold and distance requirements\n", out[2])

#print(asymmetries_single("test_extract_pattern.bed","+-", 0, 3, 2))
# test area
# works
# print overlap(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"))
# works
# print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,False,False,False)
# works
# print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,True,False,False)
# works
# print proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),0,500,False,True,False)
# works
# print strand_annotate_third_BED_overlap(read_BED("Myeloid.indels"),read_BED("Ensembl.genes_hg19_TSSs.bed"))
# works
# DataL,ScoreL=read_BED("MCF7_RepliStrand.leading_lagging.bed",True)
# separate_on_score(ScoreL,DataL,10)
# works - minor error with extra bin, needs fixing
# Strand1,Strand2,DistancesL=proximal(read_BED("All_G4.bed"),read_BED("Ensembl.genes_hg19_TSSs.bed"),"G4","hg19_TSS",0,500,False,False,10)
# print len(Strand1),len(Strand2),len(DistancesL)
# print asym_binned(0,500,10,DistancesL,Strand1,Strand2)
# asymmetries_single(read_BED("All_G4.bed"),0,1000,["++--"],0,False,0.0001,"test")
