
# coding: utf-8

# In[9]:

## Check frequencies and distributions of the severity scores
## Import dependencies

import pandas as pd
import re
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import warnings
import statsmodels.api as sm
import statsmodels.formula.api as smf
import tabulate as tb
import sys
import datetime as dt
import os

## Force the use of backend on cluster.

if sys.stdin.isatty():
    plt.switch_backend('agg')

## Define a function that takes two inputs: dat(dataframe) and col (colname). Assumes first column is ID variable 
## Function does the following:
## 1) Check the frequencies of mild, moderate and severe anxiety (GADsev) as an absolute count 
## and a proportion of total
## 2) Prints these, and summaries of total individuals in dataset, with any anxiety, and with any anxiety severity( >5)
## NOTE: this is set up quite specifically for GADsev scores. Should ignore the last count (endorse >5 items) if using
## Any other target columns

def sumstats(dat,col):
  
    ## define current date, and current working directory 
    
    now = dt.date.today()
    current = os.getcwd()
    print current    
    ## Create paths for plot directories (containing a date folder for all plot sets)
    
    pltdir = 'plots/'
    datedir = 'plots/{0}/'.format(now)
    pathplt = os.path.join(current,pltdir)
    

    ## Check if there is already a plot folder, and a date folder within that. If not, create both. 
    ## Only do this if being run from the command line
    
    if sys.stdin.isatty():
        fname = 'NewPlot'
        scriptdir = os.path.join(current, r'plots/{0}/{1}'.format(now,fname))
        path1 = os.path.join(current, r'plots/{0}/'.format(now))
        path2 = os.path.join(current, r'plots/{0}/{1}'.format(now,fname))
        if not os.path.exists(pathplt):
            os.mkdir(pltdir)
            os.mkdir(datedir)
            os.mkdir(scriptdir)
        elif not os.path.exists(path1):
            os.mkdir(path1)
            os.mkdir(path2)
        elif not os.path.exists(path2): 
            os.mkdir(path2)

    ## Try to make this flexible about running interactively versus command line. Identify if interactive node, or not
    ## If interactive, use plot.show(), if not then do not try to interactively show plots. 
    ## Save plots within current working directory plots folder, within folder named with todays date
    
    def show(title):
        if sys.stdin.isatty():
            newfig = os.path.join(path2,title)
            plt.savefig(newfig,bbox_inches='tight')
            plt.clf()                                                       ## Reset the plot area (clear figure)
        else:
            # running interactively
            return plt.show()
        return
        
    warnings.filterwarnings('ignore')

    ## Set up font dictionaries for later plots. Title and Axis dictionaries
        
    font_title = {'family': 'sans-serif',
                  'color': 'black',
                  'weight': 'bold',
                  'size': 'x-large'}
    
    font_ax = {'family': 'sans-serif',
               'color': 'black',
               'weight': 'bold',
               'style':'oblique',
               'size': 'large'}
    
    grp = range((dat[col].nunique() + 1))                          ## Count unique values.

  
    ## Get the counts and proportion of total for all items in col

    ## This option gives actual counts if False
    count = pd.Series.value_counts(dat[col], normalize = False,    
                                             sort = False,
                                             bins = None)
    ## Gives proportion of the total if True

    proportion = pd.Series.value_counts(dat[col], normalize = True, 
                                             bins = None) 
    
    ## Calculate interesting counts (total participants, total endorsing at least 1 items, total endorsing >5)
    totn = dat.iloc[:,0].count()
    totann = dat[col].count()
   
    if col == 'GADsev':                                             ## If we are working with severity, take 1 to =
        totsevn = dat[col][dat[col] >= 1].count()                   ## > 5 item endorsement
    else:
        totsevn = dat[col][dat[col] >= 5].count()
                

    ## Set up our on screen feedback for count and proportion variables
    
    cnt = "Counts by {0} score =\n\n{1}\n\n".format(col,count)
    prp = "As a proportion of total =\n\n{0}\n\n".format(proportion)
    tot = "THE TOTAL RECORDS IN THIS DATA SET ARE {0:>37}\n". format(totn)
    totsev = "the total records in this data set endorsing 5 or more items{0:>15}\n\n\n".format(totsevn)

    ## Calculate the summary stats for chosen col
    
    low = min(dat[col])
    high = max(dat[col])
    mean = dat[col].mean()
    SD = dat[col].std()
    med = dat[col].median()
    skew = dat[col].skew()
    kurt = dat[col].kurt()
    rang = '{0} - {1}'.format(low,high)
       
    ## Set up our on screen feedback for the summary statistics 
    
    sumstats = tb.tabulate([['Mean', mean],['Std Dev', SD],['Range',rang],['Median', med],['Skew',skew], 
                            ['Kurtosis',kurt]],
                            headers=["Test","Result"],
                            tablefmt="rst")
    
    sums = "\n\nSUMMARY STATISTICS FOR {0}: \n{1}\n\n".format(col,sumstats)
    
     
    ## Verbose summary of data set

    ## Print counts and proportions on the screen. Only print group counts if there are less than 6 groups
    ## If group < 6, assume group level data and draw a graph to show the counts of each group as a clear 
    ## proportion of the total number
    ## Only do this if there is a small amount of groups (i.e. group level data likely)
    
    print '\n######################################'
    print '\n'
    print '      Summary of {0}\n            {1}'.format(col,now)
    print '\n'
    print '######################################\n\n'


    
    if len(grp) > 5:
        print '\nThis appears to be continuous data. Summary statistics, histograms and QQplots will be created\n\n'
        print tot,sums
        if sys.stdin.isatty():
            print 'Histogram, histogram (log scale) and qqplot saved in: current working directory/plots/<date>.'
        else:
            print 'Histogram, histogram (log scale) and qqplot displayed below:'
    else:
        print '\nThis appears to be group level data. Groups will be shown as proportion of the total dataset\n'
        print tot,cnt,prp
        if sys.stdin.isatty():
            print 'Bar plot showing group counts saved in current working directory.'
        else:
            print 'Bar plot showing group counts displayed below:'


    
        labels = []                                                    ## Create an empty list for variable labels

        for i in grp:                                                  ## Loop through the unique values
            n = i + 1                                                  ## start from 1 (instead of zero)
            while n <= (len(grp) -1):                                  ## Stop at the last true value
                labels.append('group {0}'.format(n))                   ## Add each new value to the label list
                break                                                  ## Stop when n > last value in grp list


        bartot = [totann]*(len(grp)-1)                                 ## this creates a variable to create a bg colour 
                                                                       ## for the total cases behind each bin
        ymin = 0
        ymax = totann

        xmin = (min(grp) - 1)   
        xmax = (max(grp) +2)                                      

        plt.bar(grp[1:], count, color='b')                              ## Make the group count bars in blue (exc 0 index)
        plt.bar(grp[1:], bartot, color='r', bottom=count)               ## Overlay this on the red total bars

        plt.ylabel('Count', fontdict=font_ax)
        plt.title('{0} group counts as proportion of N'.format(col),fontdict=font_title)
        plt.ylim(ymin,ymax)                                             ## Set the y limit to the maximum individuals
        plt.xlim(xmin,xmax)                                             ## This allows padding in the margins in this case
        plt.xticks(grp[1:], labels, rotation = 'vertical')              ## Label the x axis groups
       
        figtitle = '{0}group_proportion_BarPlot{1}.pdf'.format(col,now)
             
        show(figtitle)
        
        print '\n\n     _*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*___END___*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_\n\n'    
    
    ## Now check variable distribution as a histogram and qq plot only if NOT group level data
    
    if len(grp) > 5:
        hmin = low                                    ## Set new plot limits based on values not groups
        hmax = high

        plt.hist(dat[col],
                range = [hmin,hmax],
                log = False)
        plt.title('{0} Distribution of scores'.format(col),fontdict=font_title)
        
        figtitle = '{0}_Histogram_{1}.pdf'.format(col,now)

        show(figtitle)

        plt.hist(dat[col],
                 range = [hmin,hmax],
                 log = True)
        plt.title('{0} Distribution of scores log scale'.format(col),fontdict=font_title)

        figtitle = '{0}_Histogram_log_{1}.pdf'.format(col,now)
        
        show(figtitle)
    
        sm.qqplot(dat[col],fit=True, line='45')
        plt.title('{0} qqplot'.format(col), fontdict=font_title)
        plt.suptitle('NOTE: compares against a standard normal distribution. Line at 45 degrees')
 
        figtitle = '{0}_qqplot_{1}.pdf'.format(col,now)

        show(figtitle)

        print '\n\n     _*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*___END___*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_\n\n'    

    

    return 


# In[1]:

"""Function to create a table comparing basic count summaries for different variables in a dataset"""

def SumTable(df,*args):
    
    print '\n####################################################################################'
    print '###'
    print "###   Function creating a summary table of case/control counts and percentages"
    print "###        input is a dataframe, and some list of variables to include"
    print '###'
    print "####################################################################################\n"
    print "NOTE: percent of total dataset shown in table - not excluding individuals dropped for missingness etc\n\n"
    
    
    
    totsamp = df.iloc[:,0].count()
    
    head = ["Definition","Count","Percent"] 
    template = '{:<25} {:<14} {:<6}'                           ## Create template for table
    print template.format(*head)                               ## Add headers
    print template.replace(':', ':-').format('', '', '')       ## Insert dashed line
     
    
    counts = []
    percs = []

    for i in list(*args):
        setup = list(np.where(df[i] == 1,1,0))
        co = float(setup.count(1))
        per = (co/totsamp)*100
        print template.format(i,co,per)

        
    
    return


# In[ ]:



