import os
import functions
import argparse
import wrapper_functions as wf

def fun4(motif_no_annotation, motif_annotation):
    """
    This function generates the annotated file and saves it as an output"
    """
    Annotation_data=functions.strand_annotate_third_BED_overlap(motif_no_annotation,motif_annotation)
    path_out = wf.output_path("orientation", "bed", os.path.basename(motif_no_annotation), os.path.basename(motif_annotation))
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

def orientation_wrapper(args):
    """Wrapper calling fun4 with the appropriate arguments"""
    motif_no_annotation = wf.path_checker(args.motif_no_annotation)
    motif_annotation = wf.path_checker(args.motif_annotation)
    for path1 in motif_no_annotation:
        for path2 in motif_annotation:
            fun4(path1, path2)


if __name__ == "__main__":
    args = orientation_parser()
    orientation_wrapper(args)

