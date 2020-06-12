import sys, os
import functions
import argparse
import wrapper_functions as wf
import visualizations


def fun3(args):
    # arguments
    motifsAL = wf.path_checker(args.motifsA)
    motifsBL = wf.path_checker(args.motifsB)
    motifsAL_names = wf.name_splitter(args.names_A, motifsAL)
    motifsBL_names = wf.name_splitter(args.names_B, motifsBL)
    min_distance = args.min_distance
    max_distance = args.max_distance
    expected_asym = args.expected_asym
    expected_asym_conv_div = args.expected_asym_conv_div
    upstream_only = args.upstream_only
    downstream_only = args.downstream_only
    plots = args.plots
    bins = args.bins

    orientation = args.orientation
    if orientation != None:
        paths_after_orientation = []
        for path in motifsAL:
            name_orientation = orientation.fun4(orientation, path)
            paths_after_orientation.append([name_orientation])
        motifsAL = paths_after_orientation

    # I think orientation analysis if user points to file(s) should go here before we start the asymmetries estimations

    number_of_files = len(motifsAL) * len(motifsBL)

    # Generates all pairs between motifsA and motifsB
    motif_pairs, names_pairs = functions.pairs_generator(motifsAL, motifsBL, motifsAL_names, motifsBL_names)

    p_pL = []
    m_mL = []
    m_pL = []
    p_mL = []
    p_val_same_oppositeL = []
    p_val_same_opposite_BonferoniL = []
    RatiosL = []
    p_val_conv_divergL = []
    p_val_conv_diver_BonferoniL = []
    Ratio_conv_divergL = []

    # Loops through all combinations
    for i in range(len(motif_pairs)):
        motifA = motif_pairs[i][0]
        motifB = motif_pairs[i][1]

        # for a pair of files finds the orientations
        total_asymmetries, per_bin_asymmetries = functions.proximal(motifA, motifB, names_pairs[i][0],
                                                                    names_pairs[i][1], min_distance, max_distance,
                                                                    upstream=upstream_only, downstream=downstream_only,
                                                                    bins=bins)
        Distances_orientations, p_p, m_m, p_m, m_p, same_strand, opposite_strand, convergent, divergent = total_asymmetries

        # same vs opposite analysis
        Ratio_same_opposite, p_val_same_opposite, p_val_same_opposite_Bonferoni = functions.statistical_evaluation(
            same_strand, opposite_strand, number_of_files, expected_asym=expected_asym)
        p_pL.append(p_p)
        m_mL.append(m_m)  # same orientation data
        p_mL.append(p_m)
        m_pL.append(m_p)  # opposite orientation data
        p_val_same_oppositeL.append(p_val_same_opposite)
        p_val_same_opposite_BonferoniL.append(p_val_same_opposite_Bonferoni)  # p-values for same vs opposite
        RatiosL.append(Ratio_same_opposite);  # Ratio of asymmetry for same / opposite analysis

        # convergent vs divergent analysis
        Ratio_conv_diverg, p_val_conv_diverg, p_val_conv_diver_Bonferoni = functions.statistical_evaluation(p_m, m_p,
                                                                                                            number_of_files,
                                                                                                            expected_asym=expected_asym_conv_div)
        p_val_conv_divergL.append(p_val_conv_diverg)
        p_val_conv_diver_BonferoniL.append(p_val_conv_diver_Bonferoni)
        Ratio_conv_divergL.append(Ratio_conv_diverg)

        if plots:
            # generates histogram same opposite
            visualizations.barplot_gen(same_strand, opposite_strand, "Same", "Opposite",
                                       wf.output_path("pairwise_asymmetries", "png", "same_opposite_orientation",
                                                      names_pairs[i][0], names_pairs[i][1]))
            # generates historam covergent divergent
            visualizations.barplot_gen(p_m, m_p, "Convergent", "Divergent",
                                       wf.output_path("pairwise_asymmetries", "png", "convergent_divergent_orientation",
                                                      names_pairs[i][0], names_pairs[i][1]))

            same_strandL_distance, opposite_strandL_distance, divergentL_distance, convergentL_distance = Distances_orientations
            print(same_strandL_distance)
            visualizations.distnace_distribution_gen(same_strandL_distance, opposite_strandL_distance, "Same",
                                                     "Opposite", min_distance, max_distance,
                                                     wf.output_path("pairwise_asymmetries", "png",
                                                                    "ditribution_same_opposite", names_pairs[i][0],
                                                                    names_pairs[i][1]))
            visualizations.distnace_distribution_gen(convergentL_distance, divergentL_distance, "Convergent",
                                                     "Divergent", min_distance, max_distance,
                                                     wf.output_path("pairwise_asymmetries", "png",
                                                                    "ditribution_convergent_divergent",
                                                                    names_pairs[i][0], names_pairs[i][1]))

        if bins:
            BinsL, p_p_binsL, m_m_binsL, p_m_binsL, m_p_binsL, same_strand_binsL, opposite_strand_binsL, convergent_binsL, divergent_binsL = per_bin_asymmetries
            # Here we need to decide what is the outputs we want to provide since they can be too many and complicated or focus on the plots and a small table
            # Same Opposite orientation
            visualizations.barplot_pair_lists_gen([(round(binned[0],0), round(binned[1],0)) for binned in BinsL], same_strand_binsL, opposite_strand_binsL, "Same", "Opposite",
                                                  "Distance", "Strand Orientation",
                                                  wf.output_path("pairwise_asymmetries", "png",
                                                                 "bins_same_opposite_orientation", names_pairs[i][0],
                                                                 names_pairs[i][1]))
            # Convergent Divergent orientation
            visualizations.barplot_pair_lists_gen([(round(binned[0],0), round(binned[1],0)) for binned in BinsL], convergent_binsL, divergent_binsL, "Convergent", "Divergent",
                                                  "Distance", "Strand Orientation",
                                                  wf.output_path("pairwise_asymmetries", "png",
                                                                 "bins_convergent_divergent_orientation",
                                                                 names_pairs[i][0], names_pairs[i][1]))

    # generates table <- this should be done for all pairs together.
    functions.table_gen(names_pairs, p_pL, m_mL, p_mL, m_pL, p_val_same_oppositeL, p_val_same_opposite_BonferoniL,
                        RatiosL, p_val_conv_divergL, p_val_conv_diver_BonferoniL, Ratio_conv_divergL,
                        wf.output_path("pairwise_asymmetries", "txt", "all_assymetries"))
    return


def pairwise_asymmetries_parser():
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("motifsA",
                        help="BED-formatted files. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("motifsB",
                        help="BED-formatted files. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-nA", "--names_A",
                        help="Optional argument. A name for each of the motif A files for more human-readable output. Each name must correspond to a region file path")
    parser.add_argument("-nB", "--names_B",
                        help="Optional argument. A name for each of the motif B files for more human-readable output. Each name must correspond to a motif file path")
    parser.add_argument("-o", "--orientation",
                        help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-p", "--plots", action="store_true", help="Optional flag. Display output plots")
    parser.add_argument("-ea", "--expected_asym", type=wf.check_valid_probability, default=0.5,
                        help="Optional argument. The expected asymmetry bias between the pairs of motifs. Default is 0.5")
    parser.add_argument("-ec", "--expected_asym_conv_div", type=wf.check_valid_probability, default=0.5,
                        help="Optional argument. The expected convergent / divergent asymmetry bias between the pairs of motifs Default is 0.5.")
    parser.add_argument("-min", "--min_distance", type=wf.check_positive_int, default=1,
                        help="Two consecutive motifs with distance lower than the min_distance will not be considered as significant for the purpose of this analysis. Default = 1")
    parser.add_argument("-max", "--max_distance", type=wf.check_positive_int, default=100,
                        help="Two consecutive motifs with distance higher than the max_distance will not be considered as significant for the purpose of this analysis. Default = 100")
    parser.add_argument("-b", "--bins", type=wf.check_positive_int, default=1,
                        help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1")
    direction = parser.add_mutually_exclusive_group()
    direction.add_argument("-up", "--upstream_only", action="store_true",
                           help="Perform the analysis only for occurrences of motif A upstream of occurrences of motif B, within the distance limits")
    direction.add_argument("-down", "--downstream_only", action="store_true",
                           help="Perform the analysis only for occurrences of motif A downstream of occurrences of motif B, within the distance limits")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = pairwise_asymmetries_parser()
    fun3(args)
