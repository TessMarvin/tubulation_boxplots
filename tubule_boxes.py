#Author: Tess Marvin (tmarvin@nd.edu)
#Usage:
#Purpose: Crate box and whisker plot of tubulation data for Vaughan Lab
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from argparse import ArgumentParser
import numpy as np
from scipy import stats

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
    for i in range(1,tubdata.shape[1]+1):
        y=tubdata[columns[i-1]]
        x=np.random.normal(i, 0.04, size=len(y))
        plt.scatter(x,y,alpha=0.5, color="black")
    #Here we want to do an ANOVA
    for c in columns:
        print(tubdata[c])
    f,p = stats.f_oneway(tubdata[c] for c in columns)
    plt.show()
if __name__ == '__main__':
    raw_data= './tubule_data.csv'
    tubdata = pd.read_csv(raw_data,
                          sep = ',')
    boxplot_hero(tubdata)
