#Author: Tess Marvin (tmarvin@nd.edu)
#Usage: python tubule_boxes.py
#Purpose: Crate box and whisker plot of tubulation data for Vaughan Lab
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from argparse import ArgumentParser
import numpy as np
from scipy import stats
#from statsmodels.stats.multicomp import pairwise_tukeyhsd
from collections import Counter

def boxplot_hero(tubdata):
    '''
    So this function makes box plots for the three treatment groups
    '''
    c="black"
    bp= tubdata.boxplot(grid = False, boxprops=dict(color=c), capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c))
    clevels = np.linspace(0., 1., tubdata.shape[0])
    columns=list(tubdata)
    '''
    This was for a random jitter along the x axis for each box -- didn't look as good as I had hoped

    for i in range(1,tubdata.shape[1]+1):
        y=tubdata[columns[i-1]]
        x=np.random.normal(i, 0.07, size=len(y))
        plt.scatter(x,y,alpha=0.5, color="black")
    '''
    '''
    This is a non-random jitter on the x-axis -- hopefully it looks good
    '''
    for i in range(1,tubdata.shape[1]+1):
        y=sorted(tubdata[columns[i-1]].tolist())
        coun = Counter(y)
        x=[]
        k=0.04
        pos = True
        j=0
        while j < len(y):
            if coun[y[j]] == 1:
                x.append(i)
                j=j+1
            else:
                    for n in range(coun[y[j]]):
                        if(pos):
                            x.append(i+(n*k))
                            pos= not pos
                        else:
                            x.append(i-(n*k))
                            pos= not pos
                    j = j+coun[y[j]]
        plt.scatter(x,y,alpha=0.5, color="black")
    #Edit Needed: make this a command line arg
    bp.set(ylabel='Tubule Count')
    #Edit Needed: make this a command line arg
    bp.set(title='Effect of STARD9 Knockdown on Lysosomal Membrane Tubulation')
    #Here we want to do an ANOVA because there are more than two groups
    #A list to hold the list of data points for each cell type
    li = []
    #Fill up the empty list so it looks like this [[data for cell type],[data for next cell type], [data for last]]
    for c in columns:
        li.append(tubdata[c].tolist())
    f,p = stats.f_oneway(*[li[l] for l in range(len(li))])
    print(f,p)
    plt.show()
if __name__ == '__main__':
    raw_data= './tubule_data.csv'
    tubdata = pd.read_csv(raw_data,
                          sep = ',')
    boxplot_hero(tubdata)
