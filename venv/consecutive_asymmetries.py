import sys
#import functions
import argparse
import wrapper_functions as wf

def fun2(args):
    paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)
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