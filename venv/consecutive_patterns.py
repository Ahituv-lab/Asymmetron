import sys,os
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
       names = wf.name_splitter(names)

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
        patterns = ["++","--","+-","-+"];

    if bins == None:
        bins = 1;

    if threshold == None:
        threshold = 0.05;

    ConsecutiveD_Total=[]
    # A simple way to integrate orientation in the flags here. Probably needs a lot of improvement though.
    orientation = args.orientation
    if orientation != None:
        paths_after_orientation=[];
        for path in paths:
            os.system("python orientation.py "  + path + " " + orientation)
            paths_after_orientation.append(["ASYMMETRON"+os.path.basename(path)+"_ANNOTATED_"+os.path.basename(orientation)])
        paths = paths_after_orientation
    
    for index,path in enumerate(paths):
        name = names[index]
        DataL_significant,consecutiveL,occsL,consecutive_controlL,occs_controlL = functions.asymmetries_single(path,patterns,min_distance,max_distance,threshold)
        Counter_consecutive_realL=[Counter(consecutive_pattern) for consecutive_pattern in consecutiveL]
        Counter_consecutive_controlL=[Counter(consecutive_pattern_control) for consecutive_pattern_control in consecutive_controlL]

        # Write significant results in an output file
        with open(wf.output_path("consecutive_patterns","bed",path.split("/")[-1],"statistically_siginificant"), 'w') as output_file:
            for line in DataL_significant:
                output_file.write('\t'.join([str(x) for x in line])+'\n')

        for i in range(len(Counter_consecutive_realL)):
                 consecutive, times_found = zip(*Counter_consecutive_realL[i].items())
                 consecutive_control, times_found_control = zip(*Counter_consecutive_controlL[i].items())
                 ConsecutiveD = dict(zip(consecutive, times_found))

                 ConsecutiveD_Total.append(ConsecutiveD)

                 ConsecutiveD_control = dict(zip(consecutive_control, times_found_control))
                 TimesFullList=[];TimesFullList_control=[];
                 for k in range(1,max(max(consecutive),max(consecutive_control))):
                     if k in ConsecutiveD.keys():
                         TimesFullList.append(ConsecutiveD[k])
                     else:
                         TimesFullList.append(0)

                     if k in ConsecutiveD_control.keys():
                          TimesFullList_control.append(ConsecutiveD_control[k])
                     else:
                          TimesFullList_control.append(0)
                 
                 if plots==True:
                     visualizations.barplot_single_gen(TimesFullList,range(1,len(TimesFullList)+1),"Occurrences","Consecutive occurrences",wf.output_path("consecutive_patterns","png",path.split("/")[-1],str(patterns[i])))
                     visualizations.barplot_pair_lists_gen(range(1,len(TimesFullList)+1),TimesFullList_control,TimesFullList,"Expected","Observed","Consecutive occurrences",'',wf.output_path("consecutive_patterns","png",path.split("/")[-1]+"_with_controls",str(patterns[i])))
                     # We want to show biases in distances of consecutive
                     visualizations.distribution_gen(occsL[i],occs_controlL[i],wf.output_path("consecutive_patterns","png",path.split("/")[-1],"distances_inconsecutive_pattern_"+str(patterns[i])))

                 # I think instead of Bins here it can be gradient of distances or something like that
                 if bins>1:
                     consecutiveLL_bin=[];occsLL_bin=[];consecutive_controlLL_bin=[];occs_controlLL_bin=[];
                     Bins=functions.binner(min_distance,max_distance,bins)
                     for min_bin,max_bin in Bins:
                         DataL_significant_bin,consecutiveL_bin,occsL_bin,consecutive_controlL_bin,occs_controlL_bin = functions.asymmetries_single(path,patterns[i],min_bin,max_bin,threshold)
                         consecutiveLL_bin.append(consecutiveL_bin);
                         occsLL_bin.append(occsL_bin);
                         consecutive_controlLL_bin.append(consecutive_controlL_bin);
                         occs_controlLL_bin.append(occs_controlL_bin);
                     
                     consecutiveLL_binT = np.array(consecutiveLL_bin).T.tolist();occsLL_binT = np.array(occsLL_bin).T.tolist();
                     consecutive_controlLL_binT = np.arrray(consecutive_controlLL_bin).T.tolist(); occs_controlLL_binT = np.array(occs_controlLL_bin).T.tolist()
                     # Plot barplot of occs consecutive in each bin
                     visualizations.barplot_single_gen(OccsL,OccsL,"Occurrences",wf.output_path("consecutive_patterns", "png", 'distribution_distances'))

                 # Table with all the outputs for all strands
                 functions.table_consecutive(ConsecutiveD_Total,patterns,wf.output_path("consecutive_patterns","txt",path.split("/")[-1],"_Consecutive_Patterns_Total",str(patterns[i])))

        # Need to add here vizualization as heatmap for all patterns and number of consecutive
        #ConsecutiveD_Total[i]


    return

def consecutive_patterns_parser():
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", help="Enter the path of the file to analyze. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-n", "--names", help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a path")	   # --patterns eg.. ++/+-+
    parser.add_argument("-pt","--patterns", help = "Patterns to search, comma separated. Default is ++,--,+-,-+.",
                        type=wf.check_valid_pattern)
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the "
                                                       "min_distance will not be considered as significant for the "
                                                       "purpose of this analysis. Default = 0", type=wf.check_positive_int)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the "
                                                       "max_distance will not be considered as significant for the "
                                                       "purpose of this analysis. Default = 100",
                        type=wf.check_positive_int)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int)  # Needs rephrasing
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.", type=float)
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = consecutive_patterns_parser()
    fun1(args)

