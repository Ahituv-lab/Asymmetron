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
        xtcik_number = 10
        plt.bar(range(1,len(List1)*3+1,3),List1,label=name1,align="center")
        plt.bar(range(2,len(List2)*3+1,3),List2,label=name2,align="center")
        plt.xticks(range(1,len(List1)*3+1,3*xtcik_number),[x_tickL[k] for k in range(1,len(x_tickL)*1+1,xtcik_number)],rotation=90,fontsize=9)
        plt.ylabel("Occurrences")
        plt.xlabel(x_label)
        plt.legend(frameon=False,title=title_legend)
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return

def barplot_single_gen(List1,List1_names,y_label,x_label,output):
        """ 
        This should be an option for the user if he wants to generate vizualizations too.
        """
        plot_styler()
        plt.bar(range(1,len(List1)*1+1,1),List1,align="center")
        xtick_number = 10
        plt.xticks(range(1,len(List1)*1+1,xtick_number),[List1_names[k] for k in range(1,len(List1_names)*1+1,xtick_number)],rotation=90)
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
        return


def heatmap_gen(DataLL,DataLL_control,BinsL,output):
       try:
           import seaborn as sns
       except:
            print("seaborn not imported")

       try:
           import pandas as pd
       except: 
            print("pandas not imported")

       for k in range(len(DataLL)):
          if DataLL[k]==[]:
             DataLL[k]={};
       for k in range(len(DataLL_control)):
         if DataLL_control[k]==[]:
             DataLL_control[k]={};

       all_cons = [max(k.keys()) for k in DataLL if k!={}]

       if all_cons==[]:
          return

       else:
          max_cons = max(all_cons)
       RatioLL=[]
       for i in range(len(DataLL)):
           RatioL=[];
           for k in range(1,max_cons+1):
               if k in DataLL_control[i].keys():
                   if float(DataLL_control[i][k])!=0 and float(DataLL[i][k])!=0:
                       RatioL.append(math.log10(DataLL[i][k]/float(DataLL_control[i][k])))
                   else:
                       RatioL.append(np.nan)
               else:
                   RatioL.append(np.nan)
           RatioLL.append(RatioL)

       df = pd.DataFrame(np.array(RatioLL),index=BinsL)
       mask = df.isnull()
       sns.heatmap(df,cbar_kws={'label': 'log(Enrichment)'}, mask=mask)
       plt.xlabel("Consecutive occurrences")
       plt.ylabel("Distance bins")
       plt.tight_layout()
       plt.savefig(output)
       plt.close()
       return

def distribution_gen(occsL,occsL_control,output):
       """
       This function calculates the distance of consecutive patterns and plots it for both the real data and the controls.
       """
       from collections import Counter
       plot_styler()
       Distances_consecutiveL = [occsL[k+1]-occsL[k] for k in range(len(occsL)-1)]
       Distances_consecutive_controlL = [occsL_control[k+1]-occsL_control[k] for k in range(len(occsL_control)-1)]

       Distances_consecutiveD = Counter(Distances_consecutiveL).most_common()
       Distances_consecutive_controlD = Counter(Distances_consecutive_controlL).most_common()

       plt.plot([k[0] for k in Distances_consecutiveD],[m[1] for m in Distances_consecutiveD],"o",markersize=2,label="Observed")
       plt.plot([k[0] for k in Distances_consecutive_controlD],[m[1] for m in Distances_consecutive_controlD],"o",markersize=2,label="Expected")
       plt.xlabel("Distance of consecutive")
       plt.ylabel("Occurrences")
       plt.legend(frameon=False)
       plt.savefig(output)
       plt.close()
       return

def distnace_distribution_gen(same_strandL_distance,opposite_strandL_distance,name1,name2,min_dist,max_dist,output):
      try:
          import seaborn as sns
      except: 
           print("seaborn not imported")

      plot_styler()
      plt.hist(same_strandL_distance,50,histtype='step',label=name1)
      plt.hist(opposite_strandL_distance,50,histtype='step',label=name2)
      plt.xlabel("Distance")
      plt.ylabel("Occurrences")
      plt.xlim(min_dist,max_dist)
      plt.legend(frameon=False,title="Orientation")
      plt.savefig(output)
      plt.close()
