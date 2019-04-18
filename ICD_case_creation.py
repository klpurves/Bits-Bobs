
# coding: utf-8


# Creating phenotype file for ICD code cases and controls for medication proxy paper analysis

## klp
## 04/19

## When running on cluster, inline input required is: (1) Script
## (2) File containing ICD main and secondary diagnoses columns 
## (3)  File containing IDs 

## Using this will result in file with cases set to 1, controls set to 0, 
# and all others set to NA (missing)

##Output will be 1 phenotype  files containing 2 columns (iid, phen)
## 


# In[2]:

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
import patsy
from IPython.display import HTML
import DataSummary as DS
import math
import sys
import datetime as dt
import argparse
#from rpy2 import importr
#import rpy2
#import pandas.rpy.commom as com


# In[ ]:

sys.stderr.write('\n\n####################################################################################################\n')
sys.stderr.write('###\n')
sys.stderr.write('###    Create 3 levels of anxiety phenotype file for iPSYCH anxiety and stress disorder replication\n')
sys.stderr.write('###\n')
sys.stderr.write('###                                 klp 2018\n')
sys.stderr.write('###\n')
sys.stderr.write('#######################################################################################################\n\n')


# In[ ]:

## Add argument flags to parse from command line

parser = argparse.ArgumentParser()


parser.add_argument('--ICD','--icd',action='store',required=True)
parser.add_argument('--ID','--id',action='store',required=True)

# parse the command line arguments and store as a dictionary
args = vars(parser.parse_args())


# In[2]:

print 'The following will be included'
print(args)


# In[ ]:

icd = args.get("ICD") 
id = args.get("ID")
sys.stderr.write("The icd file is set to:{0}".format(icd))


iid = pd.read_csv(id,sep=' ')

""" Create a function that takes a row in a data frame, 
and a list of possible case or control categories and identifies cases / controls """

def case(Row,list_of_categories):
    if len(set(Row).intersection(list_of_categories)) > 0:
        return 1
    else:
        return np.nan
    


# In[ ]:

""" Create a function that eliminates anyone from cases if they meet exclusion criteria """

def ExcludeCase(Row,list_of_categories):
    if len(set(Row).intersection(list_of_categories)) > 0:
        return 1
    else:
        return np.nan
    


# In[ ]:

""" Create a function that takes a row in a data frame, 
and a list of possible case or control categories and identifies cases / controls """


def control(Row,list_of_categories):
    if len(set(Row).intersection(list_of_categories)) < 1:
        return 1
    else:
        return np.nan

""" ICD case, exclusion and control lists. """

## Correspond to ICD for al;l anxiety and depression disorders

icdcase = ['"F320"','"F321"','"F322"','"F323"','"F328"','"F329"','"F330"',
'"F331"','"F332"','"F333"','"F334"','"F338"','"F339"','"F340"','"F341"',
'"F400"','"F401"','"F402"','"F403"','"F408"','"F409"', '"F410"','"F411"',
'"F412"','"F413"','"F418"','"F419"']
 
## Correspond to ICD ASD, ADHD, ED, BiP,SCZ, Pain 

icdexc = ['"F300"','"F308"','"F311"','"F314"','"F317"','"F380"','"F301"',
'"F309"','"F312"','"F315"','"F318"','"F348"','"F381"','"F302"','"F310"',
'"F313"','"F316"','"F319"',
'"F349"','"F388"',   
'"R520"', '"R521"', '"R522"', '"R529"',
'"F200"','"F201"','"F202"','"F203"','"F204"','"F205"','"F206"','"F208"','"F209"',
'"F220"','"F228"','"F229"','"F230"','"F231"','"F232"','"F233"','"F238"','"F239"',
'"F250"','"F251"','"F252"','"F258"','"F259"',
'"F840"','"F841"','"F842"','"F843"','"F844"','"F845"','"F846"','"F847"','"F848"','"F849"',
'"F500"', '"F501"', '"F502"', '"F503"', '"F504"', '"F505"', '"F508"', '"F509"',
'"F900"', '"F901"', '"F908"', '"F909"']

## All of the above + OCD and stress

icdcon = ['F320"','"F321"','"F322"','"F323"','"F328"','"F329"','"F330"',
'"F331"','"F332"','"F333"','"F334"','"F338"','"F339"','"F340"','"F341"',
'"F400"','"F401"','"F402"','"F403"','"F408"','"F409"', '"F410"','"F411"',
'"F412"','"F413"','"F418"','"F419"',
'"F300"','"F308"','"F311"','"F314"','"F317"','"F380"','"F39"','"F301"',
'"F309"','"F312"','"F315"','"F318"','"F348"','"F381"','"F302"','"F310"',
'"F313"','"F316"','"F319"',
'"F349"','"F388"',   
'"R520"', '"R521"', '"R522"', '"R529"',
'"F20"','"F200"','"F201"','"F202"','"F203"','"F204"','"F205"','"F206"',
'"F208"','"F209"',
'"F21"','"F220"','"F228"','"F229"',
'"F230"','"F231"','"F232"','"F233"','"F238"','"F239"',
'"F250"','"F251"','"F252"','"F258"','"F259"',
'"F840"','"F841"','"F842"','"F843"','"F844"','"F845"','"F846"',
'"F847"','"F848"','"F849"',
'"F500"', '"F501"', '"F502"', '"F503"', '"F504"', '"F505"', '"F508"', '"F509"',
'"F900"', '"F901"', '"F908"', '"F909"',
'"F420"', '"F421"', '"F422"', '"F428"', '"F429"',
'"F430"','"F431"','"F432"','"F433"','"F438"','"F439"'] 



print "                         *******                             \n\n"
print 'UKBB ICD 10 diagnostic codes:', icdcase 
print 'UKBB ICD 10 diagnostic autism case codes:', icdexc
print 'UKBB ICD 10 diagnostic codes for control exclusion:', icdcon 

sys.stderr.write('\n MHQ inclusion/exclusion and control list created \n Check output for UKBB data codes used \n')


# In[ ]:

"""Open the ICD dataset in a memory friendly way using open / readlines 
(this will turn each line into an easily aprsed list)"""

df = open(icd,"r")
d = df.readlines()


# In[ ]:

## Create the new phenotype column containing probable cases by 
## splitting line into component columns and performing function on each line on by one. 
## This is a memory friendly alternative to loading in and operating on entire file at once. 

icdcase_list = []

for n in range(1,len(d)):
    x = re.split(' ',d[n])
    icdcase_list.append(case(x,icdcase))

sys.stderr.write('\n ICD case list created. 1 indicates case status \n')


# In[ ]:

## Exclude any cases with autism

icdexc_list = []

for n in range(1,len(d)):
    x = re.split(' ',d[n])
    icdexc_list.append(ExcludeCase(x,icdexc))
    
    


# In[ ]:

## identify controls

icdcont_list = []

for n in range(1,len(d)):
    x = re.split(' ',d[n])
    icdcont_list.append(control(x,icdcon))

sys.stderr.write('\n ICD control list created. 1 indicates control status \n')


# In[ ]:

""" Provide basic summary of newly created data"""

# Calculate how many cases we have for anxiety.
# If there are screened controls, coutn how many of these, and total sample size.
# Also calculate how many are excluded due to comorbid Anorexia, ASD and ADHD"""

cases = icdcase_list.count(1)
exclu = icdexc_list.count(1)
conts = icdcont_list.count(1)


print "                         *******                             \n\n"

print 'For full descriptives use the data summary sumstats function\n\n'


sys.stderr.write('\n Check output for count of cases and controls from ICD \n')


## close the ICD data to save memory

df.close()
""" add ICD columns to new dataset"""

## turn them into series
icca = pd.Series(icdcase_list)
icex = pd.Series(icdexc_list)
iccon = pd.Series(icdcont_list)


sys.stderr.write('\n ICD lists now series \n')


# create dataframe
## concatenare the case and ex series as columns in a dataframe.

AllData = pd.concat([icca,icex,iccon],axis=1) 

## And add the controls

AllData.columns = ['icd.case','icd.caseexc','icd.cons']
#AllData["icd.caseexc"] = icex
#AllData["icd.cons"] = iccon


sys.stderr.write('\n ICD criteria included as columns to main dataset \n')

## And now create a column indicating if they are cases or controls

AllData["ICDcase"] = np.nan

AllData["ICDcase"] = np.where(AllData["icd.case"] == 1,1,AllData["ICDcase"])

AllData["ICDcase"] = np.where(AllData["icd.caseexc"] == 1,np.nan,AllData["ICDcase"])

AllData["ICDcase"] = np.where(AllData["icd.cons"] == 1,0, AllData["ICDcase"])

## add ID columns 

AllData["IID"] = iid

sys.stderr.write('\n ICD final case column created. See output for summary of cases and controls to ensure this matches expectations \n')


# In[ ]:

""" Provide basic summary of newly created data"""

# Calculate how many cases we have for anxiety.
# If there are screened controls, coutn how many of these, and total sample size.
# Also calculate how many are excluded due to comorbid Anorexia, ASD and ADHD"""

cases = list(np.where((AllData['ICDcase'] == 1),1,0)).count(1)
conts = list(np.where((AllData['ICDcase'] == 0),1,0)).count(1)


print "                         *******                             \n\n"
print  cases, "ICD 10 Cases (excluding autism) \n"
print  conts, "ICD 10 Controls \n "

print 'For full descriptives use the data summary sumstats function\n\n'

PhenFin = AllData[['IID','ICDcase']]

PhenFin.to_csv("/mnt/lustre/groups/ukbiobank/usr/megan/proxyukb/ICDdiagnoses",
                sep = ' ',
               header = True,
               index = False,
	       na_rep='NA',
               float_format='%.0f')


# In[ ]:

# Indicate all processes complete

sys.stderr.write("\n\nAll files successfulyl created. Process terminating...\n\n\n")

