import itertools,math
import numpy as np
from scipy.stats import binom_test
from pybedtools import BedTool
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def histogram_gen(strand1L,strand2L,bins_used,output):
        """ 
        Accepts the strand asymmetries  per bin and generates histograms
        """
        RatiosL=[ratio_calc(strand1L[i],strand2L[i]) for i in range(len(strand1L))]
        plt.bar(range(1,len(strand1L)+1,1),RatiosL,align="center",color="lightblue")
        plt.xticks(range(len(bins)),bins)
        plt.xlabel("Bins (Distance)")
        plt.ylabel("Strand Asymmetry Ratio")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return


def barplot_gen(strand1,strand2,name1,name2,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        ax = plt.subplot(111)
        plt.bar(range(1,3),[strand1,strand2],align="center")
        plt.xticks(range(1,3),[name1,name2])
        plt.ylabel("Occurrences")
        plt.xlabel("Strand Orientation")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return

def barplot_pair_lists_gen(bin_sizes_rangeL,List1,List2,name1,name2,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        ax = plt.subplot(111)
        plt.bar(range(1,len(List1)*3+1,3),List1,label=name1,align="center")
        plt.bar(range(2,len(List2)*3+1,3),List2,label=name2,align="center")
        plt.xticks(range(1,len(List1)*3+1,3),[str(bin_sizes_rangeL[k][0])+"-"+str(bin_sizes_rangeL[k][1]) for k in range(len(bin_sizes_rangeL))])
        plt.ylabel("Occurrences")
        plt.xlabel("Bins")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        plt.legend(frameon=False,title="Orientation")
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return

def barplot_single_gen(List1,List1_names,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        ax = plt.subplot(111)
        plt.bar(range(1,len(List1)*1+1,1),List1,align="center")
        plt.xticks(range(1,len(List1)*1+1,1),List1_names)
        plt.ylabel("Occurrences")
        plt.xlabel("Bins")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return


