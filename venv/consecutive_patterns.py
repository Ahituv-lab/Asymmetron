import sys
import functions
import argparse
import wrapper_functions as wf


def fun1(args):
    paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)
    for path in paths:
        Counter_consecutive_real,Counter_consecutive_control=asymmetries_single(path,min_distance,max_distance,patterns,bins=0,plot=plots,threshold)
 
        # Orientation link is missing / fun4 to be used here

    return


if __name__ == "__main__":
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Enter the path of the file to analyze. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-n", "--names", help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a path")	   # --patterns eg.. ++/+-+
    parser.add_argument("-pt","--patterns", help = "Patterns to search, comma separated. Default is ++,--,+-,-+.")
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the min_distance will not be considered as significant for the purpose of this analysis. Default = 0", type=int)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the max_distance will not be considered as significant for the purpose of this analysis. Default = 100", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    # parser.add_argument("-dc", "--expected_div_conv", help="Optional argument. Expected divergence convergence. Default = 0.5", type=float) # Not sure if necessary
    # parser.add_argument("-sp", "--expected_same_op", help="Optional argument. Expected ................ Default = 0.5", type=float)
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.", type=int)
    args = parser.parse_args()

    fun1(args)





