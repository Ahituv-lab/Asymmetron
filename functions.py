## ASYMMETRON ###
import itertools
import numpy as np
import warnings
from scipy.stats import binom_test

try:
    from pybedtools import BedTool
except:
    print("Pybedtools not imported")

try:
    import visualizations
except:
    print("visualizations not imported")

def pairs_generator(pathL1, pathL2, NamesL1, NamesL2):
    """
    :param pathL1: list of paths
    :param pathL2: list of paths
    :param NamesL1: list of names, corresponding to pathsL1
    :param NamesL2: list of names corresponding to pathsL2
    :return: tuple of two lists. List 1 consists of all possible combinations of (pathsL1, pathsL2). List 2 consists
    of the corresponding combinations of names
    """
    return list(itertools.product(pathL1, pathL2)), list(itertools.product(NamesL1, NamesL2))


def read_BED(path, last_col=False):
    """
    This function reads bed files.
    last_col=True: If an extra column for a score (e.g. replication timing or gene expression) is given
    This function returns a list of lists with the BED file information and the list of scores.
    """
    if not last_col:
        Data = []
        with open(path) as f:
            for line in f:
                Data.append(line.strip().split()[:6])
        return Data

    elif last_col:
        Data = []
        Score = []
        with open(path) as f:
            for line in f:
                Data.append(line.strip().split()[:6])
                Score.append(float(line.strip().split()[-1]))
        return Data, Score
    else:
        print("ERROR")


def binner(min_size, max_size, bin_no):
    """
    :param min_size lower bound of the interval to divide into bins
    :param max_size upper bound of the interval to divide into bins
    :param bin_no the interval will be divided into that many bins
    :returns list of tuples, representing the lower and upper limits of each subinterval after dividing into bins

    This function separates the input interval into bins
    """
    bin_size = float(max_size - min_size) / float(bin_no)
    Bins = [(min_size + bin_size * k, min_size + bin_size * (k + 1)) for k in range(0, bin_no)]
    return Bins


def separate_on_score(path_score, path, number_of_bins, number_of_files, expected_asym, expected_asym_c_d):
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
    DataStepsL = []
    ScoresStepsL = []
    for step in StepsL:
        DataStep = []
        ScoreStep = []
        for i in range(len(ScoreL)):
            if step[0] <= ScoreL[i] <= step[1]:
                DataStep += [DataL[i]]
                ScoreStep += [ScoreL[i]]
        DataStepsL += [DataStep]
        ScoresStepsL += [ScoreStep]

    # Calculates the asymmetry for same / opposite and convergent / divergent asymmetry for each bin.
    Ratio_Same_Opposite = []
    Ratio_Convergent_Divergent = []
    Binom_Test_Same_Opposite = []
    Binom_Test_Same_Opposite_Bonferoni = []
    Binom_Test_Convergent_Divergent = []
    Binom_Test_Convergent_Divergent_Bonferoni = []
    for step in range(len(StepsL)):
        p_p_step, m_m_step, p_m_step, m_p_step, same_strand_step, opposite_strand_step, convergent_step, divergent_step = overlap(
            DataStepsL[step], DataL2)
        if same_strand_step + opposite_strand_step != 0:
            Ratio_Same_Opposite.append(same_strand_step / float(same_strand_step + opposite_strand_step))
        else:
            Ratio_Same_Opposite.append(0.5)

        if p_m_step + m_p_step != 0:
            Ratio_Convergent_Divergent.append(p_m_step / float(p_m_step + m_p_step))

        binom_same_opposite = binom_test(same_strand_step, same_strand_step + opposite_strand_step, expected_asym)
        binom_same_opposite_Bonferoni = min(1, binom_same_opposite * number_of_files)

        binom_convergent_divergent = binom_test(p_m_step, p_m_step + m_p_step, expected_asym_c_d)
        binom_convergent_divergent_Bonferoni = min(1, binom_convergent_divergent * number_of_files)

        Binom_Test_Same_Opposite.append(binom_same_opposite);
        Binom_Test_Same_Opposite_Bonferoni.append(binom_same_opposite_Bonferoni);
        Binom_Test_Convergent_Divergent.append(binom_convergent_divergent);
        Binom_Test_Convergent_Divergent_Bonferoni.append(binom_convergent_divergent_Bonferoni)

    return Ratio_Same_Opposite, Ratio_Convergent_Divergent, StepsL, Binom_Test_Same_Opposite, Binom_Test_Same_Opposite_Bonferoni, Binom_Test_Convergent_Divergent, Binom_Test_Convergent_Divergent_Bonferoni


def strand_annotate_third_BED_overlap(unnotated_path, annotated_path):
    """
    For a third file that doesn't have its own annotation e.g. mutation files since mutations are on both strands
    This function enables the strand annotation of such a file
    Using an annotated file as mirror, for overlapping instances
    Returns the unnotated file, with strand annotation
    """
    DataL_unnotated = BedTool(unnotated_path)
    DataL_annotated = read_BED(annotated_path)
    Overlap_strand = DataL_unnotated.intersect(DataL_annotated, wao=True)
    Overlap_strand_df = Overlap_strand.to_dataframe()
    Chromosome = list(Overlap_strand_df.iloc[:, 0])
    Start = list(Overlap_strand_df.iloc[:, 1])
    End = list(Overlap_strand_df.iloc[:, 2])
    ID = list(Overlap_strand_df.iloc[:, 3])
    Strand = list(Overlap_strand_df.iloc[:, -2])
    Start_Annotated = list(Overlap_strand_df.iloc[:, -6])
    End_Annotated = list(Overlap_strand_df.iloc[:, -5])
    Chromosome, Start, End, ID, Strand, Start_Annotated, End_Annotated = zip(
        *((chrom, start, end, id_used, strand, start_annot, end_annot) for chrom, start, end, id_used, strand, start_annot, end_annot in
          zip(Chromosome, Start, End, ID, Strand, Start_Annotated, End_Annotated) if strand in ["+", "-"]))
    DataL = []
    for i in range(len(Chromosome)):
        if Chromosome[i] == DataL[-1][0] and Start[i] == DataL[-1][1] and End[i] == DataL[-1][2]:
            center_previous = (DataL[-1][1] + DataL[-1][2])/2
            center_current = (Start[i]+End[i])/2
            center_annotated = (Start_Annotated[i] + End_Annotated[i]) / 2
            if abs(center_annotated - center_previous) > abs(center_annotated - center_current):
                DataL[-1] = [Chromosome[i], Start[i], End[i], ID[i], ".", Strand[i]]
            else:
                pass
        else:
            DataL.append([Chromosome[i], Start[i], End[i], ID[i], ".", Strand[i]])
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
    Strand1 = list(Overlap_df.iloc[:, 5])
    Strand2 = list(Overlap_df.iloc[:, 11])
    p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = orientation(Strand1, Strand2)
    return p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent


def proximal(path1, path2, window_min, window_max, upstream=False, downstream=False, bins=None):
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
    elif upstream is True:
        closest = DataL1.closest(DataL2, D='ref', id=False, iu=True)
    elif downstream is True:
        closest = DataL1.closest(DataL2, D='ref', iu=False, id=True)
    else:
        closest = DataL1.closest(DataL2, D='ref')

    closest_df = closest.to_dataframe()
    Strand1_init = list(closest_df.iloc[:, 5])
    Strand2_init = list(closest_df.iloc[:, 11])
    Distance_init = [i for i in list(closest_df.iloc[:, -1])]
    Distance1_temp, Strand1, Strand2 = zip(
        *((dist, strand1, strand2) for dist, strand1, strand2 in zip(Distance_init, Strand1_init, Strand2_init) if
          abs(dist) <= window_max and abs(dist) >= window_min and dist >= 0))
    Distance2_temp, Strand1_temp, Strand2_temp = zip(
        *((dist, strand2, strand1) for dist, strand1, strand2 in zip(Distance_init, Strand1_init, Strand2_init) if
          abs(dist) <= window_max and abs(dist) >= window_min and dist < 0 ))
    Distance = list(Distance1_temp)+list(Distance2_temp)
    Strand1 = list(Strand1)+list(Strand1_temp)
    Strand2 = list(Strand2)+list(Strand2_temp)
    p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = orientation(Strand1, Strand2)

    # Calculate the distance distributions for all orientations
    Distances_orientations = get_distance_orientations(Distance, Strand1, Strand2, window_min, window_max)

    p_pL_bin = []
    m_mL_bin = []  # Same orientation
    p_mL_bin = []
    m_pL_bin = []  # Opposite orientation
    same_strandL_bin = []
    opposite_strandL_bin = []  # Combined same / opposite orientations
    convergentL_bin = []
    divergentL_bin = []

    Bins = []
    if bins is not None:
        # Performs the same analysis for each bin.
        Bins = binner(window_min, window_max, bins)
        for index, bin_i in enumerate(Bins):
            Strand1Bin = []
            Strand2Bin = []
            min_bin, max_bin = bin_i
            for k in range(len(Distance)):
                if Distance[k] >= min_bin and Distance[k] < max_bin:
                    Strand1Bin.append(Strand1[k])
                    Strand2Bin.append(Strand2[k])

            p_p_bin, m_m_bin, p_m_bin, m_p_bin, same_strand_bin, opposite_strand_bin, convergent_bin, divergent_bin = orientation(
                Strand1Bin, Strand2Bin)
            p_pL_bin.append(p_p_bin)
            m_mL_bin.append(m_m_bin)  # Same orientation, per bin
            p_mL_bin.append(p_m_bin)
            m_pL_bin.append(m_p_bin)  # Opposite orientation per bin
            same_strandL_bin.append(same_strand_bin)
            opposite_strandL_bin.append(opposite_strand_bin)
            convergentL_bin.append(convergent_bin)
            divergentL_bin.append(divergent_bin)

    return (Distances_orientations, p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent), (
        Bins, p_pL_bin, m_mL_bin, p_mL_bin, m_pL_bin, same_strandL_bin, opposite_strandL_bin, convergentL_bin,
        divergentL_bin)


def get_distance_orientations(DistanceL, Strand1L, Strand2L, window_min, window_max):
    same_strandL_distance = []
    opposite_strandL_distance = []
    divergentL_distance = []
    convergentL_distance = []
    for index in range(len(Strand1L)):
        if (DistanceL[index] <= window_max and DistanceL[index] >= window_min):
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

    return (same_strandL_distance, opposite_strandL_distance, divergentL_distance, convergentL_distance)


def orientation(sign1L, sign2L):
    """
    This function takes as input two lists of signs and calculates their relative position (same strand,
    opposite strand, convergent and divergent)
    :param signs1L list of signs
    :param signs2L list of sings
    :return p_p: number of times both corresponding signs were in + strand
    :return m_m: number of times both corresponding signs were in - strand
    :return p_m: first sign in + strand, second sign in - strand
    :return m_p: first sign in - strand, second sign in + strand
    :return same_strand: sum of p_p+m_m
    :return opposite_strand: sum of p_m + m_p
    :return convergent: same as p_m
    :return divergent: same as m_p
    """
    p_p = 0
    m_m = 0
    p_m = 0
    m_p = 0
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
    same_strand = p_p + m_m
    opposite_strand = p_m + m_p
    convergent = p_m
    divergent = m_p
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
    if the user expects a background bias then that should be given as an expected ratio, different than zero
    (that would alter our bionimal test, and we could just  change 0.5 to a variable)
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
        "Feature_1" + '\t' + "Feature_2" + "\t" + "plus_plus" + '\t' + "minus_minus" + '\t' + "plus_minus" + '\t' + "minus_plus" + '\t' + "p_value_same_opposite" + '\t' + "p-value_same_opposite_Bonferoni_corrected" + '\t' + "Ratio_same_opposite" + '\t' + "p_value_divergent_convergent" + '\t' + "p_value_divergent_convergent Bonferoni corrected" + '\t' + "Ratio divergent convergent" + '\n')
    for i in range(len(NamesL_pairs)):
        datafile.write(
            NamesL_pairs[i][0] + '\t' + NamesL_pairs[i][1] + '\t' + str(p_pL[i]) + '\t' + str(m_mL[i]) + '\t' + str(
                p_mL[i]) + '\t' + str(m_pL[i]) + '\t' + str(p_valsL[i]) + '\t' + str(p_vals_BonferoniL[i]) + '\t' + str(
                RatiosL[i]) + '\t' + str(p_valsL_divergent_convergent[i]) + '\t' + str(
                p_valsL_divergent_convergent_BonferoniL[i]) + '\t' + str(RatiosL_divergent_convergent[i]) + '\n')
    datafile.close()
    return


def table_bins_gen(Score_names, Ratio_Bins, Ratio_Convergent_Divergent_Bins, Binom_Test_Same_Opposite_Bins,
                   Binom_Test_Same_Opposite_Bonferoni_Bins, Binom_Test_Convergent_Divergent_Bins,
                   Binom_Test_Convergent_Divergent_Bonferoni_Bins, output_table):
    datafile = open(output_table, "w")
    datafile.write("Score_Range_Min" + '\t' + str(
        "Score_Range_Max") + '\t' + "Feature_1" + '\t' + "Feature_2" + "\t" + "Ratio_same_opposite" + '\t' + "Ratio_divergent_convergent" + '\n')
    for i in range(len(Score_names)):
        datafile.write(str(round(Score_names[i][0], 0)) + "-" + str(round(Score_names[i][1], 0)) + '\t' + str(
            Ratio_Bins[i]) + '\t' + str(Ratio_Convergent_Divergent_Bins[i]) + '\t' + str(
            Binom_Test_Same_Opposite_Bins[i]) + '\t' + str(Binom_Test_Same_Opposite_Bonferoni_Bins[i]) + '\t' + str(
            Binom_Test_Convergent_Divergent_Bins[i]) + '\t' + str(
            Binom_Test_Convergent_Divergent_Bonferoni_Bins[i]) + '\n')

    datafile.close()
    return


def table_consecutive(ConsecutiveL, namesL, path_out):
    for k in range(len(ConsecutiveL)):
        if not ConsecutiveL[k]:
            ConsecutiveL[k] = {}

    consecutives_all = ([max(k.keys()) for k in ConsecutiveL if k != {}])
    if consecutives_all:
        max_consecutive = max(consecutives_all)
    else:
        max_consecutive = 0

    with open(path_out, 'w') as output:
        output.write(str("Number of consecutive occurrences") + '\t' + '\t'.join(
            [str(x) for x in range(1, max_consecutive + 1)]) + '\n')

        for i in range(len(ConsecutiveL)):
            if ConsecutiveL[i] != {}:
                ConsecutiveLT = sorted(ConsecutiveL[i].items())
                consecutive, times_found = zip(*ConsecutiveLT)
                consecutiveL = list(consecutive)
                times_foundL = list(times_found)
            else:
                consecutiveL = []
                times_foundL = []

            for k in range(1, max_consecutive + 1):
                if k not in consecutiveL:
                    consecutiveL = consecutiveL[:k - 1] + [k] + consecutiveL[k - 1:]
                    times_foundL = times_foundL[:k - 1] + [0] + times_foundL[k - 1:]

            output.write(str(namesL[i]) + '\t' + '\t'.join([str(x) for x in times_foundL]) + '\n')
    return


def write_BED_out(DataL, path_out):
    # Write significant results in an output file
    with open(path_out, 'w') as output_file:
        for line in DataL:
            output_file.write('\t'.join([str(x) for x in list(line)]) + '\n')
    return


def find_sub_str(my_str, sub_str):
    start = 0
    while True:
        start = my_str.find(sub_str, start)
        if start == -1:
            return
        yield start
        start += len(sub_str)


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
    :return consecutiveL: A dictionary where keys are the number of consecutive appearances of the pattern and values
                          how many times they appeared, e.g. {3:2} means that there were there were two instances in the
                          file where the pattern appeared three consecutive times
    :return distancesL: A list which includes the distances between consecutive appearances of the pattern
    :return DataL_significant: A list of list (same format as DataL), which includes all rows of the initial file that
                               were found to be part of a significant appearance of the pattern. e.g. if rows 10-20
                               include 5 conesecutive appearances of the pattern, which gives a p-value higher than
                               the threshold, then those rows will be included in DataL_significant. The
                               calculated p-value is appended at the end of each row.
    """
    from collections import defaultdict
    DataL = [list(line) for line in DataL]
    n = len(pattern)

    # Create a string with all signs (5th column in the sorted BED file)
    signs = ""
    for line in DataL:
        signs += line[5]

    if pattern == "basic":
        same_total=0;opposite_total=0;
        same=0;opposite=0;
        SameL = defaultdict(int); OppositeL= defaultdict(int);
        DistancesL_same=[];DistancesL_opposite=[];DataSignificantL_same=[];DataSignificantL_opposite=[];
        DataL_temp = [DataL[0]]

        for i in range(1, len(DataL)):
            distance = max(0, int(DataL[i][1]) - int(DataL[i-1][2]))
            if (distance >= min_distance and distance <= max_distance):
    
                if DataL[i-1][5]==DataL[i][5] and DataL[i-1][5] in ["+","-"]:
                    same+=1;same_total+=1
                    DistancesL_same.append(distance)

                    if opposite>0:
                        OppositeL[opposite]+=1
                        if 0.5**opposite<threshold:
                            DataSignificantL_opposite+=DataL_temp
                        DataL_temp=[];

                    DataL_temp+=[DataL[i]]
                    opposite=0;
                        
               
                if (DataL[i-1][5]=="+" and DataL[i][5]=="-") or (DataL[i-1][5]=="-" and DataL[i][5]=="+"): 
                    opposite+=1;opposite_total+=1;
                    DistancesL_opposite.append(distance)

                    if same>0:
                        SameL[same]+=1
                        if 0.5**same<threshold:
                            DataSignificantL_same+=DataL_temp
                        DataL_temp=[];
                     
                    DataL_temp+=[DataL[i]]
                    same=0;

            else:
            
                if same > 0:
                    SameL[same]+=1
                    if 0.5**same<threshold:
                        DataSignificantL_same+=DataL_temp

                if opposite > 0:
                    OppositeL[opposite]+=1
                    if 0.5**opposite<threshold:
                        DataSignificantL_opposite+=DataL_temp

                DataL_temp=[];
                same=0;opposite = 0;
     	
        return SameL,OppositeL,DistancesL_same,DistancesL_opposite,DataSignificantL_same,DataSignificantL_opposite,same_total,opposite_total

    # Find all occurences of the pattern in the string of signs without accounting for distances
    occs = list(find_sub_str(signs, pattern))

    # Remove occurences that do not meet the distance criterion
    for i in range(len(occs)):
        index = occs[i]
        for j in range(n - 1):
            distance = max(0, int(DataL[index + j + 1][1]) - int(DataL[index + j][2]))
            if distance < min_distance or distance > max_distance:
                occs[i] = None
                break
    occs = [x for x in occs if x is not None]

    # If no occurence of the pattern was found, return empty lists
    if not occs:
        return [], [], []

    # Here we translate the probability threshold to consecutive occcurrences
    number_of_tests = len(DataL)
    total_plus = signs.count("+")
    total_minus = signs.count("-")
    probability = {}
    probability["+"] = total_plus / float(total_plus + total_minus)
    probability["-"] = 1 - probability["+"]
    probability_pattern = np.prod([probability[k] for k in list(pattern)])
    # Warn the user if probability_pattern == 1
    if probability_pattern == 1:
        msg = "The probability of this pattern  being found at any given point of the file was calculated to be 1. " \
              "This usually would happen in degenerate cases (e.g. looking for pattern \"+\" in a file consisting " \
              "only of + and could result in unpredictable outputs"
        warnings.warn(msg)
    # Lowest consecutive number of patterns that with probability of appearning lower than the threshold (p-value)
    from math import ceil, log
    consecutive_threshold = ceil(log(threshold / number_of_tests) / log(probability_pattern))

    # Filter for number of consecutive occurences that meet the threshold criterion and are within the distance window
    DataL_significant = []
    counter = 1
    DataL_temp = DataL[occs[0]:occs[0] + n];
    consecutiveL = defaultdict(int);
    for i in range(1, len(occs)):
        index = occs[i]
        distance = max(0, int(DataL[index][1]) - int(DataL[index - 1][2]))
        if occs[i] - occs[i - 1] == n and (distance >= min_distance and distance <= max_distance):
            counter += 1
            DataL_temp.extend(DataL[index:index + n])
        else:
            consecutiveL[counter] += 1
            if counter >= consecutive_threshold:
                # Add p_value to list for the number of consecutive occurences of the pattern
                p_value = (probability_pattern ** counter) * number_of_tests
                for line in DataL_temp:
                    line.append(p_value)
                DataL_significant.extend(DataL_temp)
            DataL_temp.clear()
            counter = 1

    # Add last lines
    if counter >= consecutive_threshold:
        p_value = (probability_pattern ** counter) * number_of_tests
        for line in DataL_temp:
            line.append(p_value)
        DataL_significant.extend(DataL_temp)
        consecutiveL[counter] += 1

    distancesL = []
    for i in range(len(occs) - 1):
        index = occs[i]
        index_next = occs[i + 1]
        distance = max(0, int(DataL[occs[i + 1]][1]) - int(DataL[occs[i] + n - 1][2]))
        distancesL.append(distance)

    return consecutiveL, distancesL, DataL_significant


def asymmetries_single(path, patternsL, min_distance, max_distance, threshold,simulations):
    print("consecutive asymmetries running on " +str(simulations)+" simulations")
    from random import shuffle
    DataL = list(read_BED(path))

    # Here we want to shuffle the strand column
    ColL = [k[5] for k in DataL]
    ColL_copy = [k[5] for k in DataL]

    consecutiveL = []
    occsL = []
    DataL_significantL = []
    consecutive_controlL = []
    occs_controlL = []
    DataL_significant_controlL = []

     
    consecutive_controlL=[];occs_controlL=[];DataL_significant_controlL=[];sameL_c=[];oppositeL_c=[];
    if patternsL==["basic"]:


        bias_same=0;bias_opposite=0;
        sameL,oppositeL, DistancesL_same, DistancesL_opposite, DataL_significant_Same,DataL_significant_Opposite,same_total,opposite_total = extract_pattern(DataL, "basic", min_distance, max_distance, threshold)
        consecutiveL.append(sameL);consecutiveL.append(oppositeL);
        occsL.append(DistancesL_same);occsL.append(DistancesL_opposite);
        DataL_significantL.append(DataL_significant_Same)
        DataL_significantL.append(DataL_significant_Opposite)

        Ratio_same_opposite = same_total/float(same_total+opposite_total)

        for simulation in range(simulations):

            shuffle(ColL_copy)
            DataL_random = []
            for i in range(len(DataL)):
                line = DataL[i][:-1] + [ColL_copy[i]]
                DataL_random.append(line)
            
            same_c,opposite_c,Distances_same_c, Distances_opposite_c, Data_significant_c_same,Data_significant_c_opposite,same_total_control,opposite_total_control  = extract_pattern(DataL_random, "basic", min_distance, max_distance, threshold)
            

            # only keep significant control occurrences consecutive for one simulation
            if simulation+1 == simulations:
                consecutive_controlL.append(same_c);consecutive_controlL.append(opposite_c);
                occs_controlL.append(Distances_same_c);occs_controlL.append(Distances_opposite_c);

            Ratio_same_opposite_c = same_total_control/float(same_total_control+opposite_total_control)
            if Ratio_same_opposite>Ratio_same_opposite_c:
                bias_same+=1
            else:
                bias_opposite+=1

        p_val= min(2*min(max(1,bias_same)/float(bias_same+bias_opposite),max(1,bias_opposite)/float(bias_same+bias_opposite)),1)  
        print("empirical p-value is "+str(p_val))

    else:
        for pattern in patternsL:
            consecutive_surplus=0;consecutive_shortage=0;

            consecutive, occs, DataL_significant = extract_pattern(DataL, pattern, min_distance, max_distance, threshold)
            consecutiveL.append(consecutive)
          
            for simulation in range(simulations):
                shuffle(ColL_copy)
                DataL_random = []
                for i in range(len(DataL)):
                    line = DataL[i][:-1] + [ColL_copy[i]]
                    DataL_random.append(line)

                consecutive_control, occs_control, DataL_significant_control = extract_pattern(DataL_random, pattern,min_distance, max_distance,threshold)
                if sum(consecutive_control.values())>sum(consecutive.values()):
                    consecutive_shortage+=1
                else:
                    consecutive_surplus+=1

                if simulation+1 == simulations:
                    occs_controlL.append(occs_control)
                    consecutive_controlL.append(consecutive_control)


            occsL.append(occs)
            DataL_significantL.append(DataL_significant)

            p_val= min(2*min(max(1,consecutive_surplus)/float(consecutive_surplus+consecutive_shortage),max(1,consecutive_shortage)/float(consecutive_surplus+consecutive_shortage)),1)
            print(consecutive_surplus,consecutive_shortage)
            print("empirical p-value for pattern " +str(pattern)+" is "+str(p_val))

    
    return DataL_significantL, consecutiveL, occsL, consecutive_controlL, occs_controlL


# Ensures that code below is not run when this file is imported into another file
if __name__ == "__main__":
     pass
