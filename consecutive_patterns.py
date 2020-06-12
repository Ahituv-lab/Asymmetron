import sys,os
import numpy as np
import functions
import orientation
import argparse
import wrapper_functions as wf
try:
    import visualizations
except:
    print("visualisations not imported")

from collections import Counter

def fun1(args):

    paths = wf.path_checker(args.paths)
    names = wf.name_splitter(args.names, paths)


    min_distance = args.min_distance
    max_distance = args.max_distance
    plots = args.plots
    patterns = [x.strip() for x in args.patterns.split(',')]
    bins = args.bins
    threshold = args.threshold

    ConsecutiveD_Total=[]
    # A simple way to integrate orientation in the flags here. Probably needs a lot of improvement though.
    ort = args.orientation

    if ort is not None:
        paths_after_orientation=[]
        for path in paths:
            name_orientation=orientation.fun4(path,ort)
            paths_after_orientation.append([name_orientation])  # Are the square brackets here needed?
        paths = paths_after_orientation
    
    for index,path in enumerate(paths):
        name = names[index]
        DataL_significant,consecutiveL,occsL,consecutive_controlL,occs_controlL = functions.asymmetries_single(path,patterns,min_distance,max_distance,threshold)
        
        # Table with all the outputs for all strands, number of consecutive occurrences found for each strand pattern
        functions.table_consecutive(consecutiveL,patterns,wf.output_path("consecutive_patterns","txt",os.path.basename(path),"_Consecutive_Patterns_Total"))

        for i in range(len(patterns)):

            # Write significant results in an output file
            functions.write_BED_out(DataL_significant[i],wf.output_path("consecutive_patterns","bed",os.path.basename(path),"statistically_siginificant",patterns[i]))

      
            
            consecutiveLT = sorted(consecutiveL[i].items())
            if consecutiveLT!=[]:
                consecutive,times_found = zip(*consecutiveLT) 
            else:
                consecutive=[];times_found=[]

            consecutive_controlLT = sorted(consecutive_controlL[i].items())
            if consecutive_controlLT!=[]:
                consecutive_control,times_found_control = zip(*consecutive_controlLT)
            else:
                consecutive_control=[];times_found_control=[]

            if plots==True:
                     # Plot number of consecutive occurrences of the pattern
                     consecutive_total = list(consecutive)
                     times_found_total = list(times_found)

                     consecutive_control_total = list(consecutive_control)
                     times_found_control_total = list(times_found_control)

                     print(consecutive_control_total,consecutive_total)
                     if consecutive_control_total!=[] and consecutive_total!=[]:
                         max_consecutive = max(max(consecutive_control_total),max(consecutive_total));
                     elif consecutive_total!=[]:
                         max_consecutive = max(consecutive_total);
                     elif consecutive_control_total!=[]:
                         max_consecutive = max(consecutive_control_total)
                     else:
                         max_consecutive =0;

                     # Adding the consecutive occurrences not found
                     for k in range(1,max_consecutive+1):
                         if k not in consecutive_total:
                             consecutive_total = consecutive_total[:k-1]+[k] + consecutive_total[k-1:]
                             times_found_total = times_found_total[:k-1]+[0] + times_found_total[k-1:]
 
                         if k not in consecutive_control_total:
                             consecutive_control_total = consecutive_control_total[:k-1] + [k]+consecutive_control_total[k-1:]
                             times_found_control_total = times_found_control_total[:k-1] + [0]+times_found_control_total[k-1:]

                     visualizations.barplot_single_gen(times_found_total,consecutive_total,"Occurrences","Consecutive occurrences",wf.output_path("consecutive_patterns","png",os.path.basename(path),str(patterns[i])))

                     # Plot number  of consecutive occurrences of the pattern in the real data and in the scrambled version
                     visualizations.barplot_pair_lists_gen(consecutive_total,times_found_control_total,times_found_total,"Expected","Observed","Consecutive occurrences",'',wf.output_path("consecutive_patterns","png",os.path.basename(path)+"_with_controls",str(patterns[i])))

                     # We want to show biases in distances of consecutive
                     occsL_filtered=[];occs_controlL_filtered=[]
                     for occ in occsL[i]:
                         if min_distance<= occ <= max_distance:
                             occsL_filtered.append(occ)

                     for occ_c in occs_controlL[i]:
                         if min_distance<= occ_c <= max_distance:
                             occs_controlL_filtered.append(occ_c)

                     visualizations.distribution_gen(occs_controlL_filtered,occsL_filtered,wf.output_path("consecutive_patterns","png",os.path.basename(path),"distances_inconsecutive_pattern_"+str(patterns[i])))

        
        if bins>1:
            Bins=functions.binner(min_distance,max_distance,bins)
            for pat in patterns:
                consecutiveLL_bin=[];occsLL_bin=[];consecutive_controlLL_bin=[];occs_controlLL_bin=[];
                for min_bin,max_bin in Bins:
                    DataL_significant_bin,consecutiveL_bin,occsL_bin,consecutive_controlL_bin,occs_controlL_bin = functions.asymmetries_single(path,[pat],min_bin,max_bin,threshold)
                    consecutiveLL_bin.append(consecutiveL_bin[0])
                    occsLL_bin.append(occsL_bin)
                    consecutive_controlLL_bin.append(consecutive_controlL_bin[0]);
                    occs_controlLL_bin.append(occs_controlL_bin)
                     
                consecutiveLL_binT = np.array(consecutiveLL_bin).T.tolist();occsLL_binT = np.array(occsLL_bin).T.tolist();
                consecutive_controlLL_binT = np.array(consecutive_controlLL_bin).T.tolist(); occs_controlLL_binT = np.array(occs_controlLL_bin).T.tolist()

                functions.table_consecutive(consecutiveLL_bin,["Bin:"+str(bin_range[0])+"-"+str(bin_range[1]) for bin_range in Bins],wf.output_path("consecutive_patterns","txt",os.path.basename(path),"_Consecutive_Patterns_bins",str(pat)))

                if plots:
                    # Plot barplot of consecutive in each bin  
                    visualizations.heatmap_gen(consecutiveLL_bin,consecutive_controlLL_bin,Bins,wf.output_path("consecutive_patterns", "png", 'distribution_distances',pat))

    return

def consecutive_patterns_parser():
    # Note: Optional arguments have a - or -- in front of them
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", help="Enter the path of the file to analyze. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-n", "--names", help="Optional argument. A name for each of the motif files for more human-readable output. Each name must correspond to a path")	   # --patterns eg.. ++/+-+
    parser.add_argument("-pt","--patterns", help = "Patterns to search, comma separated. Default is ++,--,+-,-+.",
                        type=wf.check_valid_pattern, default="++,--,+-,-+")
    parser.add_argument("-min", "--min_distance", help="Two consecutive motifs with distance lower than the "
                                                       "min_distance will not be considered as significant for the "
                                                       "purpose of this analysis. Default = 0",
                        type=wf.check_positive_int, default=0)
    parser.add_argument("-max", "--max_distance", help="Two consecutive motifs with distance higher than the "
                                                       "max_distance will not be considered as significant for the "
                                                       "purpose of this analysis. Default = 100",
                        type=wf.check_positive_int, default=100)
    parser.add_argument("-b", "--bins", help="Optional argument. Split output data and graphs in the specified number of bins. Default = 1", type=int, default=1)  # Needs rephrasing
    parser.add_argument("-p", "--plots", help="Optional flag. Display output plots", action="store_true")
    parser.add_argument("-o", "--orientation", help="Optional argument. Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations. Can enter multiple paths as a comma separated string, e.g. \"path1, path2\"")
    parser.add_argument("-t", "--threshold", help="Optional argument. Threshold of p-value of consecutive patterns to save in new BED file.", type=wf.check_valid_probability, default=0.05)
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = consecutive_patterns_parser()
    fun1(args)

