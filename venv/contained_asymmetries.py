import sys,os
import functions
import argparse
import wrapper_functions as wf
import visualizations


def fun2(args):
	#paths, orientation_paths, names = wf.sanitize(args.motifs, args.orientation, args.names)
	# this is used for Bonferoni correction

	# Missing link to fun4 if provided

	# Folder to save all outputs
	directory = "outputs_contained_asymmetries"
	if not os.path.exists(directory):
		os.makedirs(directory)

        motifsL=args.motifs.split(",")
        regionsL=args.regions.split(",")

	motifsL_names=args.names_motifs
        if motifsL_names==None:
            motifsL_names=[k.split("/")[-1] for k in motifsL]
        else:
            motifsL_names = motifsL_names.split(",")

        regionsL_names=args.names_regions
        if regionsL_names==None:
            regionsL_names=[k.split("/")[-1] for k in regionsL]
        else:
            regionsL_names=regionsL_names.split(",")

        expected_asym=args.expected_asym
        if expected_asym==None:
            expected_asym=0.5
 
        number_of_files = len(motifsL) * len(regionsL)
	# All possible pairs between region files and motif files
	motif_region_pairs, names_pairs = functions.pairs_generator(motifsL, regionsL, motifsL_names, regionsL_names)

	# Save the results for the final table
	p_pL = [];
	m_mL = [];
	p_mL = [];
	m_pL = [];
	same_strandL = [];
	opposite_strandL = [];
	convergentL = [];
	divergentL = [];
	p_val_same_oppositeL = [];
	p_val_same_opposite_BonferoniL = [];
	Ratio_same_oppositeL = [];
	p_val_conv_divergL = [];
	p_val_conv_diver_BonferoniL = [];
	Ratio_conv_divergL = [];

	# Perform all comparisons of each pair
	for i in range(len(motif_region_pairs)):

		p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = functions.overlap(motif_region_pairs[i][0], motif_region_pairs[i][1])

		p_pL.append(p_p);
		m_mL.append(m_m);  # Same strand orientation
		p_mL.append(p_m);
		m_pL.append(m_pL);  # Opposite strand orientation

		same_strandL.append(same_strand);
		opposite_strandL.append(opposite_strand);
		convergentL.append(convergent);
		divergentL.append(divergent);

		# same vs opposite analysis
		Ratio_same_opposite, p_val_same_opposite, p_val_same_opposite_Bonferoni = functions.statistical_evaluation(
			same_strand, opposite_strand, number_of_files, expected_asym=expected_asym)
		Ratio_same_oppositeL.append(Ratio_same_opposite);
		p_val_same_oppositeL.append(p_val_same_opposite);
		p_val_same_opposite_BonferoniL.append(p_val_same_opposite_Bonferoni)

		# convergent vs divergent analysis
		Ratio_conv_diverg, p_val_conv_diverg, p_val_conv_diver_Bonferoni = functions.statistical_evaluation(p_m, m_p,
		                                                                                                    number_of_files,
		                                                                                                    expected_asym=expected_asym_conv_div)
		p_val_conv_divergL.append(p_val_conv_diverg);
		p_val_conv_divergL.append(p_val_conv_diver_Bonferoni);
		Ratio_conv_divergL.append(Ratio_conv_diverg);

		if plots:
			# generates histogram same opposite, we need to decide the output1
                        visualizations.barplot_gen(same_strand, opposite_strand, wf.output_path("contained_asymmetries", "same_opposite_orientation.png", names_pairs[0],names_pairs[1]))
			# generates historam covergent divergent, we need to decide the output2
                        visualizations.barplot_gen(p_m, m_p, os.path.join(directory, names_pairs[0] + "_" + names_pairs[1] + "_divergent_convergent_orientations.png"))

		if score:
			# Here we need to decide if we want to include the score for both -regions and -motifs and perform the analyses separately, score needs to go with number of score_bins
			if score_regions != False:
				Ratio_Same_Opposite, Score_names = separate_on_score(regions, motifs, number_of_bins)
				visualizations.barplot_single_gen(Ratio_Same_Opposite, Score_names, output_plot)
			if score_motifs != False:
				Ratio_Same_Opposite, Score_names = separate_on_score(motifs, regions, number_of_bins)
				visualizations.barplot_single_gen(Ratio_Same_Opposite, Score_names, output_plot)

		# I think we don't need bins here. Only type of bins would have been spliting the regions in sub-parts and doing the analysis in those but I don't think it adds much.
		# if bins:
		#    pass

	# generates table <- this should be done for all pairs.
	functions.table_gen(names_pairs, p_pL, m_mL, p_mL, m_pL, p_valsL, p_vals_BonferoniL, RatiosL, p_val_conv_divergL,
	                    p_val_conv_diver_BonferoniL, Ratio_conv_divergL)

	return


if __name__ == "__main__":
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
	parser.add_argument("-or", "--orientation_regions",
	                    help="Optional argument. Orient file(s) relative to annotated BED-formated region file(s) and perform the analysis for the un-annoated file with the new annotations. ")
	parser.add_argument("-om", "--orientation_motifs",
	                    help="Optional argument. Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations. ")
	parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
	parser.add_argument("-ea", "--expected_asym",
	                    help="Optional argument. The expected asymmetry bias between the regions and the motifs. Default is 0.5",
	                    type=float)
	parser.add_argument("-ec", "--expected_asym_conv_div",
	                    help="Optional argument. The expected convergent / divergent asymmetry bias between the regions and the motifs. Default is 0.5.",
	                    type=float)
	parser.add_argument("-s", "--score",
	                    help="Optional flag. If provided, assumes the last column of the region files is a scoring metric and uses it to subdivide the analysis into quartiles",
	                    action="store_true")
	parser.add_argument("-t", "--threshold",
	                    help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.",
	                    type=int)
	parser.add_argument("-b", "--bins",
	                    help="Optional argument. Number of bins to subdivide the results into. Only runs when --score is provided. Default value is 10.",
	                    type=wf.check_positive_int)
	args = parser.parse_args()
	print(args.bins)

	fun2(args)
