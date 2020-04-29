import sys
#import functions
import argparse
import wrapper_functions as wf

def fun4(args):
    """
    This function generates the annotated file and saves it as an output"
    """
    motif_no_annotation = wf.path_checker(args.motif_no_annotation)
    motif_annotation = wf.path_checker(args.motif.annotation)
    Annotation_data=strand_annotate_third_BED_overlap(motif_no_annotation,motif_annotation)
    datafile=open(output,"w")
    for i in Annotation_data:
        datafile.write('\t'.join([str(x) for x in i])+'\n')
    datafile.close()
    return


if __name__ == "__main":
    parser.add_argument("motif_no_annotation", help="BED-formatted file")
    parser.add_argument("motif_annotation", help="BED-formatted file")
    args = parser.parse_args()
    fun4(args)
