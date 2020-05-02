import sys,os
import functions
import argparse
import wrapper_functions as wf

def fun3(args):

    paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)
    
    # I think orientation analysis if user points to file(s) should go here before we start the asymmetries estimations


    output_path = wf.output_path("pairwise_asymmetries")

    number_of_files= len(motifsA)*len(motifsB)

    # Generates all pairs between motifsA and motifsB
    motif_pairs,names_pairs=functions.pairs_generator(pathL1,pathL2,NamesL1,NamesL2)

    p_pL=[];m_mL=[];m_pL=[];p_mL=[];p_valsL=[];p_vals_BonferoniL=[];RatiosL=[];p_val_conv_divergL=[];p_val_conv_diver_BonferoniL=[];Ratio_conv_divergL=[];

    # Loops through all combinations
    for i in range(len(motif_pairs)):
            # for a pair of files finds the orientations
            p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent=functions.proximal(path1,path2,name1,name2,min_distance,max_distance,upstream=upstream_only,downstream=downstream_only,in_parts=bins)
            #same vs opposite analysis
            Ratio_same_opposite,p_val_same_opposite,p_val_same_opposite_Bonferoni=functions.statistical_evaluation(same_strand,opposite_strand,number_of_files,expected_asym=expected_asym)
            p_pL.append(p_p);m_mL.append(m_m); # same orientation data
            p_mL.append(p_m);m_pL.append(m_p); # opposite orientation data
            p_valsL.append(p_val_same_opposite);p_val_same_opposite_BonferoniL.append(p_val_same_opposite_Bonferoni); # p-values for same vs opposite
            RatiosL.append(Ratio_same_opposite); # Ratio of asymmetry for same / opposite analysis

            # convergent vs divergent analysis
            Ratio_conv_diverg,p_val_conv_diverg,p_val_conv_diver_Bonferoni=functions.statistical_evaluation(p_m,m_p,number_of_files,expected_asym=expected_asym_conv_div)
            p_val_conv_divergL.append(p_val_conv_diverg);p_val_conv_diver_BonferoniL.append(p_val_conv_diver_BonferoniL);
            Ratio_conv_divergL.append(Ratio_conv_diverg);
	    if plots:
                 # generates histogram same opposite, we need to decide the output1
                 functions.barplot_gen(same_strand,opposite_strand,os.path.join(directory, names_pairs[0]+"_"+names_pairs[1]+ "_same_opposite_orientations.png"))
                 # generates historam covergent divergent, we need to decide the output2
                 functions.barplot_gen(p_m,m_p,os.path.join(directory, names_pairs[0]+"_"+names_pairs[1]+ "_divergent_convergent_orientations.png"))

            # If bins is true I already put in functions.proximal that it generates two barplots. Also consider a table to be generated. Also, we need to put the output of that in the same directory as outputs_pairwise_asymmetries
        if bins:
            # Here we need to decide what is the outputs we want to provide since they can be too many and complicated or focus on the plots and a small table
            pass

    # generates table <- this should be done for all pairs together.
    functions.table_gen(names_pairs,p_pL,m_mL,p_mL,m_pL,p_valsL,p_vals_BonferoniL,RatiosL,p_val_conv_divergL,p_val_conv_diver_BonferoniL,Ratio_conv_divergL)

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
    direction =parser.add_mutually_exclusive_group()
    direction.add_argument("-up", "--upstream_only", help="Perform the analysis only for occurrences of motif A upstream of occurrences of motif B, within the distance limits", type=int)
    direction.add_argument("-down", "--downstream_only", help="Perform the analysis only for occurrences of motif A downstream of occurrences of motif B, within the distance limits", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    args = parser.parse_args()

    fun3(args)
