import sys
import functions
import argparse
import wrapper_functions as wf
import visualizations
from collections import Counter

def fun1(args):
    #paths, orientation_paths, names = wf.sanitize (args.path, args.orientation, args.names)

    paths = args.paths.split(",")
    names = args.names
    if names == None:
        names = paths
    else: 
       names = names.split(",")

    min_distance = args.min_distance
    if min_distance == None:
        min_distance = 0;

    max_distance = args.max_distance
    if max_distance == None:
        max_distance = 100;

    patterns = args.patterns

    bins = args.bins
    plots = args.plots
    threshold = args.threshold

    if patterns == None:
        patterns = ["++","--","+-","-+"]

    if bins == None:
        bins = 1

    if threshold == None:
        threshold = 0.05
	
    for index,path in enumerate(paths):
        name = names[index]
        consecutiveL,occsL,consecutive_controlL,occs_controlL =functions.asymmetries_single(path,patterns,min_distance,max_distance,threshold)
        Counter_consecutive_realL=[Counter(consecutive_pattern) for consecutive_pattern in consecutiveL]

        if plots==True:
            for i in range(len(Counter_consecutive_realL)):
                 consecutive, times_found = zip(*Counter_consecutive_realL[i].items())
                 print(consecutive, times_found)
                 ConsecutiveD = dict(zip(consecutive, times_found))
                 TimesFullList=[];
                 for k in range(1,max(consecutive)):
                     if k in ConsecutiveD.keys():
                         TimesFullList.append(ConsecutiveD[k])
                     else:
                         TimesFullList.append(0)
                 visualizations.barplot_single_gen(range(1,len(TimesFullList)+1),TimesFullList,wf.output_path("consecutive_patterns", ".png", ''))

                 if bins>1:
                     Bins=functions.binner(min_distance,max_distance,bins)
                     for min_bin,max_bin in Bins:
                         consecutiveL_bin,occsL_bin,consecutive_controlL_bin,occs_controlL_bin =functions.asymmetries_single(path,patternsL,min_bin,max_bin,threshold)
                     # Plot barplot of occs consecutive in each bin
                     visualizations.barplot_single_gen(OccsL,OccsL,wf.output_path("consecutive_patterns", ".png", ''))

            # Need to add here vizualization as heatmap for all patterns and number of consecutive


        # Orientation link is missing / fun4 to be used here

    return


if __name__ == "__main__":
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", help="Enter the path of the file to analyze. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-n", "--names", help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a path")	   # --patterns eg.. ++/+-+
    parser.add_argument("-pt","--patterns", help = "Patterns to search, comma separated. Default is ++,--,+-,-+.")
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the min_distance will not be considered as significant for the purpose of this analysis. Default = 0", type=int)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the max_distance will not be considered as significant for the purpose of this analysis. Default = 100", type=int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.", type=float)
    args = parser.parse_args()

    fun1(args)

