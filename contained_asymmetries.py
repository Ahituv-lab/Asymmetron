import functions
import orientation
import argparse
import wrapper_functions as wf
try:
    import visualizations
except:
    print("visualisations not imported")


def fun2(args):

        regionsL=wf.path_checker(args.regions)
        motifsL=wf.path_checker(args.motifs)
        motifsL_names=wf.name_splitter(args.names_motifs, motifsL)
        regionsL_names=wf.name_splitter(args.names_regions, regionsL)
        expected_asym=args.expected_asym
        expected_asym_conv_div = args.expected_asym_conv_div
        plots=args.plots
        score = args.score
        bins_score = args.bins_score

        # Orientation of third file using the motifs
        orientation_motifs = args.orientation_motifs
        if orientation_motifs!= None:
            wf.path_checker(orientation_motifs)  # Test if orientation_motifs path exists
            paths_after_orientation=[]
            for path in motifsL:
                 name_orientation= orientation.fun4(orientation, path)
                 paths_after_orientation.append([name_orientation])
            motifsL = paths_after_orientation

    
        number_of_files = len(motifsL) * len(regionsL)
	# All possible pairs between region files and motif files
        motif_region_pairs, names_pairs = functions.pairs_generator(motifsL, regionsL, motifsL_names, regionsL_names)

	# Save the results for the final table
        p_pL = [];m_mL = [] # Same orientation
        p_mL = [];m_pL = []  # Opposite orientation
        same_strandL = []; opposite_strandL = [] # Summed
        convergentL = []; divergentL = []
        p_val_same_oppositeL = []; p_val_same_opposite_BonferoniL = []
        Ratio_same_oppositeL = [];Ratio_conv_divergL = []
        p_val_conv_divergL = [];p_val_conv_diver_BonferoniL = []

        # Perform all comparisons of each pair
        for i in range(len(motif_region_pairs)):

            p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(motif_region_pairs[i][0], motif_region_pairs[i][1])

            p_pL.append(p_p);m_mL.append(m_m);  # Same strand orientation
            p_mL.append(p_m);m_pL.append(m_p);  # Opposite strand orientation

            same_strandL.append(same_strand);opposite_strandL.append(opposite_strand)
            convergentL.append(convergent);divergentL.append(divergent)

            # same vs opposite analysis
            Ratio_same_opposite, p_val_same_opposite, p_val_same_opposite_Bonferoni = functions.statistical_evaluation(same_strand, opposite_strand, number_of_files, expected_asym = expected_asym)
            Ratio_same_oppositeL.append(Ratio_same_opposite)
            p_val_same_oppositeL.append(p_val_same_opposite)
            p_val_same_opposite_BonferoniL.append(p_val_same_opposite_Bonferoni)

            # convergent vs divergent analysis
            Ratio_conv_diverg, p_val_conv_diverg, p_val_conv_diver_Bonferoni = functions.statistical_evaluation(p_m, m_p,number_of_files,expected_asym=expected_asym_conv_div)
            p_val_conv_divergL.append(p_val_conv_diverg)
            Ratio_conv_divergL.append(Ratio_conv_diverg)
            p_val_conv_diver_BonferoniL.append(p_val_conv_diver_Bonferoni)

            if plots:
                # generates histogram same opposite
                visualizations.barplot_gen(same_strand, opposite_strand,"Same","Opposite", wf.output_path("contained_asymmetries","png", "same_opposite_orientation", names_pairs[i][0],names_pairs[i][1]))

                # generates historam covergent divergent
                visualizations.barplot_gen(p_m, m_p,"Convergent","Divergent", wf.output_path("contained_asymmetries", "png","convergent_divergent_orientation", names_pairs[i][0],names_pairs[i][1]))

            if score:
                Ratio_Bins,Ratio_Convergent_Divergent_Bins,Score_names,Binom_Test_Same_Opposite,Binom_Test_Same_Opposite_Bonferoni,Binom_Test_Convergent_Divergent,Binom_Test_Convergent_Divergent_Bonferoni = functions.separate_on_score(motif_region_pairs[i][1], motif_region_pairs[i][0], bins_score,number_of_files,expected_asym,expected_asym_conv_div)
                functions.table_bins_gen(Score_names,Ratio_Bins,Ratio_Convergent_Divergent_Bins,Binom_Test_Same_Opposite,Binom_Test_Same_Opposite_Bonferoni,Binom_Test_Convergent_Divergent,Binom_Test_Convergent_Divergent_Bonferoni,wf.output_path("contained_asymmetries","txt", "Table_Strand_Asymmetries_Scores", names_pairs[i][0],names_pairs[i][1]))
                 
                visualizations.barplot_single_gen(Ratio_Bins, [(int(round(score_name[0],0)),int(round(score_name[1],0))) for score_name in Score_names], "Strand Asymmetry","Score", wf.output_path("contained_asymmetries","png", "same_opposite_orientation_separated_score", names_pairs[i][0],names_pairs[i][1]))
                visualizations.barplot_single_gen(Ratio_Convergent_Divergent_Bins, [(int(round(score_name[0],0)),int(round(score_name[1],0))) for score_name in Score_names], "Strand Asymmetry","Score", wf.output_path("contained_asymmetries","png", "convergent_divergent_orientation_separated_score", names_pairs[i][0],names_pairs[i][1]))

        # generates table <- this should be done for all pairs.
        functions.table_gen(names_pairs, p_pL, m_mL, p_mL, m_pL, p_val_same_oppositeL, p_val_same_opposite_BonferoniL, Ratio_same_oppositeL, p_val_conv_divergL,p_val_conv_diver_BonferoniL,Ratio_conv_divergL,wf.output_path("contained_asymmetries","txt","all_assymetries"))

        return

def contained_asymmetries_parser():
	# Note: Optional arguments have a - or -- in front of them
        parser = argparse.ArgumentParser()
        parser.add_argument("regions",
	                    help="BED-formatted files, containing the regions within which to estimate motif asymmetries. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
        parser.add_argument("motifs",
	                    help="Enter the path of the file for each of which the asymmetries are calculated. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
        parser.add_argument("-nr", "--names_regions",
	                    help="Optional argument. A name for each of the region files for more human-readable output. Each name must correspond to a region file path")
        parser.add_argument("-nm", "--names_motifs",
	                    help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a motif file path")
        parser.add_argument("-om", "--orientation_motifs",
	                    help="Optional argument. Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations. ")
        parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
        parser.add_argument("-ea", "--expected_asym",
	                    help="Optional argument. The expected asymmetry bias between the regions and the motifs. Default is 0.5",
	                    type=float, default=0.5)
        parser.add_argument("-ec", "--expected_asym_conv_div",
	                    help="Optional argument. The expected convergent / divergent asymmetry bias between the regions and the motifs. Default is 0.5.",
	                    type=float, default=0.5)
        parser.add_argument("-s", "--score",
	                    help="Optional flag. If provided, assumes the last column of the region files is a scoring metric and uses it to subdivide the analysis into quartiles",
	                    action="store_true")
        parser.add_argument("-bs", "--bins_score",
                            help="Optional flag. Number of bins to separate the score into",
                            type=wf.check_positive_int, default=10)
        args = parser.parse_args()
        return args

if __name__ == "__main__":
    args = contained_asymmetries_parser()
    fun2(args)


