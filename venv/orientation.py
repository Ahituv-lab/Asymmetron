import sys
#import functions
import argparse
import wrapper_functions as wf

def fun4(args):
    motif_no_annotation = wf.path_checker(args.motif_no_annotation)
    motif_annotation = wf.path_checker(args.motif.annotation)

if __name__ == "__main":
    parser.add_argument("motif_no_annotation", help="BED-formatted file")
    parser.add_argument("motif_annotation", help="BED-formatted file")
    args = parser.parse_args()
    fun4(args)