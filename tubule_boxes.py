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
from statsmodels.stats.multicomp import pairwise_tukeyhsd

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
    #Now we need to check if p is less than or equal to 0.05 -- if so -- run a TUKEY TukeyHSDResults
    if(p <= 0.05):
        #We need to build a new dataframe
        df_t= pd.DataFrame(columns=['Cell_Type','Tubulation'])
        columns = list(tubdata)
        for c in range(len(columns)):
            for j in range(tubdata.shape[0]):
                new_row={'Cell_Type':str(columns[c]), 'Tubulation':int(tubdata.iloc[j,c])}
                df_t=df_t.append(new_row, ignore_index=True)
        #print(df_t)
        res = pairwise_tukeyhsd(pd.to_numeric(df_t['Tubulation']),df_t['Cell_Type'], alpha=0.05)
        #now save the results as a dataframe:
        df_res = pd.DataFrame(data=res._results_table.data[1:], columns=res._results_table.data[0])
        print(df_res)
        #So now we want to put these statistical annotations on the graph
        for c in range(len(columns)):
            rand = 0
            num=[3,7]
            for j in range(c+1,len(columns)):
                p_val = (df_res.loc[((df_res['group1']==columns[c]) | (df_res['group1']==columns[j])) \
                & ((df_res['group2']== columns[j]) | (df_res['group2']== columns[c])), ['p-adj']])
                #So here we are saying the bar should be above which boxes
                x1, x2 = c+1, j+1
                max1, max2 = max(tubdata[columns[c]].tolist()), max(tubdata[columns[j]].tolist())
                abs_max=max(max1,max2)+num[rand]
                rand=rand+1
                #y is the y coordinate of the bar annotation (20 units above the highest data point in the subplot)
                #h is how far down to draw the tips of the bar downward towards the data (it looks like a line with a small taper down on each side)
                y, h, col = abs_max, .5, 'k'
                #This part draws the bar annotation above the boxplots
                bp.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                final_p = p_val.iloc[0].to_list()
                if(final_p[0] > 0.05):
                    bp.text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color=col)
                #if t test is significant indicate with * and indicate the p-value is less than 0.05
                elif(final_p[0]<=0.05 and final_p[0] > 0.01):
                    if(final_p[0]==0.05):
                        bp.text((x1+x2)*.5, y+h, "* (p = 0.05)", ha='center', va='bottom', color=col)
                    else:
                        bp.text((x1+x2)*.5, y+h, "* (p < 0.05)", ha='center', va='bottom', color=col)
                #if t test is very significant indicate with ** and indicate the p-value is less than 0.01
                elif(final_p[0]<=0.01):
                    if(final_p[0]==0.01):
                        bp.text((x1+x2)*.5, y+h, "** (p = 0.01)", ha='center', va='bottom', color=col)
                    else:
                        bp.text((x1+x2)*.5, y+h, "** (p < 0.01)", ha='center', va='bottom', color=col)
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
