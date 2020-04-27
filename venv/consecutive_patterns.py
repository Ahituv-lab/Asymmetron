import sys
#import functions
import argparse
import os

def path_checker(paths):
    """:param paths Input given by the user as a string of paths, comma seperated
       :return a list of all paths given by the user (a list with one element if only one path was given)

       This functions takes the input paths given by the user, verifies that they exist and returns a list of paths
    """
    pathsL = [x.strip() for x in paths.split(',')]
    for path in pathsL:
        if not os.path.exists(path):
            err= "File \"" + path + "\" was not found"
            raise IOError(err)
    return pathsL

def fun1():
    return


if __name__ == "__main__":
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Enter the path of the file to analyze. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-n", "--names", help="The name of each of the motif files")
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the min_distance will not be considered as significant for the purpose of this analysis. Default = 0", type=int)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the max_distance will not be considered as significant for the purpose of this analysis. Default = infinite", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    # parser.add_argument("-dc", "--expected_div_conv", help="Optional argument. Expected divergence convergence. Default = 0.5", type=float) # Not sure if necessary
    # parser.add_argument("-sp", "--expected_same_op", help="Optional argument. Expected ................ Default = 0.5", type=float)
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations.")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.")
    args = parser.parse_args()
    paths = path_checker(args.path)  # Converts the input to a list of paths. List can include only one element, if one path is given by the user
    # Converts the input to a list of path for the orientation files. If this optional argument was not given, the variable is set to None
    try:
        orientation_paths = path_checker(args.orientation)
    except AttributeError:
        orientation_paths = None




