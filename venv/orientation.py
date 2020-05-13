import sys
import functions
import argparse
import wrapper_functions as wf

def fun4(args):
    """
    This function generates the annotated file and saves it as an output"
    """
    motif_no_annotation = wf.path_checker(args.motif_no_annotation)
    motif_annotation = wf.path_checker(args.motif_annotation)
    for path in motif_no_annotation:
        for path2 in motif_annotation:
            Annotation_data=functions.strand_annotate_third_BED_overlap(path,path2)
            path_out = "ASYMMETRON"+path+"_ANNOTATED_"+path2
            with open(path_out, "w") as f:
                for i in Annotation_data:
                    f.write('\t'.join([str(x) for x in i])+'\n')
    return path_out

def orientation_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("motif_no_annotation", help="BED-formatted file")
    parser.add_argument("motif_annotation", help="BED-formatted file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    orientation_parser()
    fun4(args)
