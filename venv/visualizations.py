import itertools,math
import numpy as np
from scipy.stats import binom_test
from pybedtools import BedTool
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_styler():
        ax = plt.subplot(111)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        return

def histogram_gen(strand1L,strand2L,bins_used,output):
        """ 
        Accepts the strand asymmetries  per bin and generates histograms
        """
        plot_styler()
        RatiosL=[ratio_calc(strand1L[i],strand2L[i]) for i in range(len(strand1L))]
        plt.bar(range(1,len(strand1L)+1,1),RatiosL,align="center",color="lightblue")
        plt.xticks(range(len(bins)),bins)
        plt.xlabel("Bins (Distance)")
        plt.ylabel("Strand Asymmetry Ratio")
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return


def barplot_gen(strand1,strand2,name1,name2,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        plot_styler()
        plt.bar(range(1,3),[strand1,strand2],align="center")
        plt.xticks(range(1,3),[name1,name2])
        plt.ylabel("Occurrences")
        plt.xlabel("Strand Orientation")
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return

def barplot_pair_lists_gen(x_tickL,List1,List2,name1,name2,x_label,title_legend,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        plot_styler()
        plt.bar(range(1,len(List1)*3+1,3),List1,label=name1,align="center")
        plt.bar(range(2,len(List2)*3+1,3),List2,label=name2,align="center")
        plt.xticks(range(1,len(List1)*3+1,3),x_tickL)
        plt.ylabel("Occurrences")
        plt.xlabel(x_label)
        plt.legend(frameon=False,title=title_legend)
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return

def barplot_single_gen(List1,List1_names,x_label,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        plot_styler()
        plt.bar(range(1,len(List1)*1+1,1),List1,align="center")
        plt.xticks(range(1,len(List1)*1+1,1),List1_names)
        plt.ylabel("Occurrences")
        plt.xlabel(x_label)
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return


def heatmap_gen(DataLL,output):
       import seaborn as sns
       import pandas as pd
       df = pd.DataFrame(np.array(DataLL))
       sns.heatmap(df)
       sns.savefig(output)
       sns.close()
       return

def distribution_gen(occsL,occsL_control,output):
       from collections import Counter
       plot_styler()
       occsD=Counter(occsL)
       occsD_control=Counter(occsL_control)

       occsD_sorted=sorted(occsD.items(), key=lambda k: -k[0]) 
       occsD_sorted_control=sorted(occsD_control.items(), key=lambda k: -k[0])

       plt.plot([k[0] for k in occsD_sorted],[m[1] for m in occsD_sorted],"o")
       plt.plot([k[0] for k in occsD_sorted_control],[m[1] for m in occsD_sorted_control],"o")
       plt.savefig(output)
       plt.close()
