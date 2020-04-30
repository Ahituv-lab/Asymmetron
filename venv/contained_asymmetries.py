import sys
import functions
import argparse
import wrapper_functions as wf

def fun2(args):

    paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)
    number_of_files= len(motifs)*len(regions)

    directory = "outputs_contained_asymmetries"
    if not os.path.exists(directory):
        os.makedirs(directory)
  
    motif_region_pairs,names_pairs=functions.pairs_generator(pathL1,pathL2,NamesL1,NamesL2)

    p_pL=[];m_mL=[];p_mL=[];m_pL=[];same_strandL=[];opposite_strandL=[];convergentL=[];divergentL=[];
    p_valsL=[];p_vals_BonferoniL=[];RatiosL=[];p_val_conv_divergL=[];p_val_conv_diver_BonferoniL=[];Ratio_conv_divergL=[];
    for i in range(len(motif_region_pairs)): 

        p_p,m_m,p_m,m_p,same_strand,opposite_strand,convergent,divergent=functions.overlap(motifs,regions)
        p_pL.apend(p_p);m_mL.append(m_m);p_mL.append(p_m);m_pL.append(m_pL);
        same_strandL.append(same_strand);opposite_strandL.append(opposite_strand);convergentL.append(convergent);divergentL.append(divergent);
        #same vs opposite analysis
        Ratio_same_opposite,p_val_same_opposite,p_val_same_opposite_Bonferoni=functions.statistical_evaluation(same_strand,opposite_strand,number_of_files,expected_asym=expected_asym)
        # convergent vs divergent analysis
        Ratio_conv_diverg,p_val_conv_diverg,p_val_conv_diver_Bonferoni=functions.statistical_evaluation(p_m,m_p,number_of_files,expected_asym=expected_asym_conv_div)
        p_val_conv_divergL.append(p_val_conv_diverg);p_val_conv_divergL.append(p_val_conv_diver_Bonferoni);Ratio_conv_divergL.append(Ratio_conv_diverg);

        if plots:
            # generates histogram same opposite, we need to decide the output1
            functions.barplot_gen(same_strand,opposite_strand,os.path.join(directory, names_pairs[0]+"_"+names_pairs[1]+ "_same_opposite_orientations.png"))
            # generates historam covergent divergent, we need to decide the output2
            functions.barplot_gen(p_m,m_p,os.path.join(directory, names_pairs[0]+"_"+names_pairs[1]+ "_divergent_convergent_orientations.png"))

    # generates table <- this should be done for all pairs.
    functions.table_gen(names_pairs,p_pL,m_mL,p_mL,m_pL,p_valsL,p_vals_BonferoniL,RatiosL,p_val_conv_divergL,p_val_conv_diver_BonferoniL,Ratio_conv_divergL)
   
    return

if __name__ == "__main__":
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("regions", help="BED-formatted files, containing the regions within which to estimate motif asymmetries. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("path", help="Enter the path of the file for each of which the asymmetries are calculated. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-nr", "--names_regions", help="Optional argument. A name for each of the region files for more human-readable output. Each name must correspond to a region file path")
    parser.add_argument("-nm", "--names_motifs", help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a motif file path")
    parser.add_argument("-or", "--orientation_regions", help="Optional argument. Orient file(s) relative to annotated BED-formated region file(s) and perform the analysis for the un-annoated file with the new annotations. ")
    parser.add_argument("-om", "--orientation_motifs", help="Optional argument. Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations. ")
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    parser.add_argument("-ea", "--expected_asym", help="Optional argument. The expected asymmetry bias between the regions and the motifs. Default is 0.5", type=float)
    parser.add_argument("-ec", "--expected_asym_conv_div", help="Optional argument. The expected convergent / divergent asymmetry bias between the regions and the motifs. Default is 0.5.", type=float)
    parser.add_argument("-s", "--score", help="Optional flag. If provided, assumes the last column of the region files is a scoring metric and uses it to subdivide the analysis into quartiles", action="store_true")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Number of bins to subdivide the results into. Only runs when --score is provided. Default value is 10.", type=int)
    args = parser.parse_args()

    fun2(args)
