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
    else:
        min_distance = int(min_distance)

    max_distance = args.max_distance
    if max_distance == None:
        max_distance = 100;
    else:
        max_distance = int(max_distance)

    patterns = args.patterns

    bins = args.bins
    plots = args.plots
    threshold = args.threshold

    if patterns == None:
        patterns = ["++","--","+-","-+"];

    if bins == None:
        bins = 1;
    else:
        bins = int(bins);

    if threshold == None:
        threshold = 0.05;

    # A simple way to integrate orientation in the flags here. Probably needs a lot of improvement though.
    orientation = args.orientation
    if orientation != None:
        paths_after_orientation=[];
        for path in paths:
            os.system("python orientation.py "  + path + " " + orientation)
            paths_after_orientation.append(["ASYMMETRON"+path.split("/")[-1]+"_ANNOTATED_"+orientation.split("/")[-1]])
        paths = paths_after_orientation
    
    for index,path in enumerate(paths):
        name = names[index]
        DataL_significant,consecutiveL,occsL,consecutive_controlL,occs_controlL = functions.asymmetries_single(path,patterns,min_distance,max_distance,threshold)
        Counter_consecutive_realL=[Counter(consecutive_pattern) for consecutive_pattern in consecutiveL]
        Counter_consecutive_controlL=[Counter(consecutive_pattern_control) for consecutive_pattern_control in consecutive_controlL]

        # Write significant results in an output file
        with open(wf.output_path("consecutive_patterns",path.split("/")[-1],"bed","statistically_siginificant"), 'w') as output_file:
            for line in DataL_significant:
                output_file.write('\t'.join([str(x) for x in line])+'\n')


        for i in range(len(Counter_consecutive_realL)):
                 consecutive, times_found = zip(*Counter_consecutive_realL[i].items())
                 consecutive_control, times_found_control = zip(*Counter_consecutive_controlL[i].items())
                 ConsecutiveD = dict(zip(consecutive, times_found))
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
                     visualizations.barplot_single_gen(TimesFullList,range(1,len(TimesFullList)+1),"Consecutive occurrences",wf.output_path("consecutive_pattern_"+str(patterns[i]), "png", ''))
                     visualizations.barplot_pair_lists_gen(range(1,len(TimesFullList)+1),TimesFullList_control,TimesFullList,"Expected","Observed","Consecutive occurrences",'',wf.output_path("consecutive_pattern_with_control_"+str(patterns[i]),"png", ''))

                 # I think instead of Bins here it can be gradient of distances or something like that
                 if bins>1:
                     consecutiveLL_bin=[];occsLL_bin=[];consecutive_controlLL_bin=[];occs_controlLL_bin=[];
                     Bins=functions.binner(min_distance,max_distance,bins)
                     for min_bin,max_bin in Bins:
                         DataL_significant_bin,consecutiveL_bin,occsL_bin,consecutive_controlL_bin,occs_controlL_bin = functions.asymmetries_single(path,patterns,min_bin,max_bin,threshold)
                         consecutiveLL_bin.append(consecutiveL_bin);
                         occsLL_bin.append(occsL_bin);
                         consecutive_controlLL_bin.append(consecutive_controlL_bin);
                         occs_controlLL_bin.append(occs_controlL_bin);
                     
                     consecutiveLL_binT = np.array(consecutiveLL_bin).T.tolist();occsLL_binT = np.array(occsLL_bin).T.tolist();
                     consecutive_controlLL_binT = np.arrray(consecutive_controlLL_bin).T.tolist(); occs_controlLL_binT = np.array(occs_controlLL_bin).T.tolist()
                     print(consecutive_controlLL_binT,"bins")
                     # Plot barplot of occs consecutive in each bin
                     #visualizations.barplot_single_gen(OccsL,OccsL,wf.output_path("consecutive_patterns", ".png", ''))

                # Need to add here vizualization as heatmap for all patterns and number of consecutive
                

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

