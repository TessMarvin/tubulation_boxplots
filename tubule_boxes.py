#Author: Tess Marvin (tmarvin@nd.edu)
#Usage: pythonw tubule_boxes.py
#Purpose: Crate box and whisker plot of tubulation data for Vaughan Lab
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from argparse import ArgumentParser
import numpy as np
from scipy import stats
from gooey import Gooey
from gooey import GooeyParser
from argparse import ArgumentParser
from collections import Counter
#from statsmodels.stats.multicomp import pairwise_tukeyhsd

def boxplot_hero(tubdata, yaxis= 'Tubule Count', fig_title= 'Effect of STARD9 Knockdown on Lysosomal Membrane Tubulation'):
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
    bp.set(ylabel= yaxis)
    #Edit Needed: make this a command line arg
    bp.set(title= fig_title)
    #Here we want to do an ANOVA because there are more than two groups
    #A list to hold the list of data points for each cell type
    li = []
    #Fill up the empty list so it looks like this [[data for cell type],[data for next cell type], [data for last]]
    for c in columns:
        li.append(tubdata[c].tolist())
    f,p = stats.f_oneway(*[li[l] for l in range(len(li))])
    print(f,p)
    plt.show()
#So this turns the command line arguments into a beautiful GUI
#Here I built out a File Menu with an About Menu
@Gooey(
    program_name='Tubulation Analysis',
    menu=[{
    'name':'File',
    'items': [{
            'type': 'AboutDialog',
            'menuTitle': 'About',
            'name': 'Tubulation Analysis and Figure Generation',
            'description': 'A tool to create visuals of data concerning lysosomal tubulation',
            'version': '1.0',
            'copyright': '2020',
            'website': 'https://github.com/TessMarvin',
            'developer': 'Tess Marvin (tmarvin@nd.edu)',
            'license': 'University of Notre Dame'
    }]
    }]
)
def main():
    #So first we will handle the arguments that are "required" -- the files and the gene of interest
    #Here we give our GUI a title
    parser = GooeyParser(description="Dashboard for Lysosomal Tubulation Analysis")
    #Here we allow for the selection of the data files to analyze
    #Because there is no - infront of file_chooser it is required!
    parser.add_argument("file_chooser", help = 'Choose the csv file to analyze.', widget='FileChooser')
    parser.add_argument("-Graph_title", help="Choose the label for the Figure.")
    parser.add_argument("-y_axis_title", help="Choose the label for the y-axis.")
    #Now we parse all these arguments
    args = parser.parse_args()
    raw_data= args.file_chooser
    yaxis=args.y_axis_title
    fig_title=args.Graph_title
    #if no CSV file is found and we want to visualize, fail gracefully and request the CSV file be provided
    if(len(raw_data) == 0):
        print("Please ensure that you select files to analyze")
        return(None)
    #if the user would like to provide y-axis and figure title, use it, otherwise use the default 
    else:
        tubdata = pd.read_csv(raw_data, sep = ',')
        if(yaxis is not None):
            if(fig_title is not None):
                boxplot_hero(tubdata,yaxis,fig_title)
            else:
                boxplot_hero(tubdata,yaxis)
        elif(fig_title is not None):
            boxplot_hero(tubdata,fig_title=fig_title)
        else:
            boxplot_hero(tubdata)

if __name__ == '__main__':
    main()
