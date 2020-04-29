import sys
#import functions
import argparse
import wrapper_functions as wf

def fun3(args):
    paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)
    number_of_files= len(motifsA)*len(motifsB)
    # for a pair of files finds the orientations
    if bins==False:
        p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent=functions.proximal(path1,path2,min_distance,max_distance,upstream=upstream_only,downstream=downstream_only,in_parts=bins)
        #same vs opposite analysis
        Ratio_same_opposite,p_val_same_opposite,p_val_same_opposite_Bonferoni=functions.statistical_evaluation(same_strand,opposite_strand,number_of_files,expected_asym=expected_asym)
        # convergent vs divergent analysis
        Ratio_conv_diverg,p_val_conv_diverg,p_val_conv_diver_Bonferoni=functions.statistical_evaluation(p_m,m_p,number_of_files,expected_asym=expected_asym_conv_div)
        # generates table <- this should be done for all pairs.
	if plots:
            # generates histogram same opposite, we need to decide the output1
            functions.barplot_gen(same_strand,opposite_strand,output1)
            # generates historam covergent divergent, we need to decide the output2
            functions.barplot_gen(p_m,m_p,output2)
    return

if __name__ == "__main__":
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("motifsA", help="BED-formatted files. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("motifsB", help="BED-formatted files. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-nA", "--names_A", help="Optional argument. A name for each of the motif A files for more human-readable output. Each name must correspond to a region file path")
    parser.add_argument("-nB", "--names_B", help="Optional argument. A name for each of the motif B files for more human-readable output. Each name must correspond to a motif file path")
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    parser.add_argument("-ea", "--expected_asym", help="Optional argument. The expected asymmetry bias between the pairs of motifs. Default is 0.5", type=float)
    parser.add_argument("-ec", "--expected_asym_conv_div", help="Optional argument. The expected convergent / divergent asymmetry bias between the pairs of motifs Default is 0.5.", type=float)
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the min_distance will not be considered as significant for the purpose of this analysis. Default = 0", type=int)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the max_distance will not be considered as significant for the purpose of this analysis. Default = 100", type=int)
    parser.add_argument("-up", "--upstream_only", help="Perform the analysis only for occurrences of motif A upstream of occurrences of motif B, within the distance limits", type=int)
    parser.add_argument("-down", "--downstream_only", help="Perform the analysis only for occurrences of motif A downstream of occurrences of motif B, within the distance limits", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    args = parser.parse_args()

    fun3(args)
