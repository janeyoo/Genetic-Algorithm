import csv
import glob
import numpy as np

"""
DETERMINES THE AVERAGE INITIAL FITNESS AT GENERATION ZERO
WITH N RUNS
"""
path =r'GA Results Data'
allFiles = glob.glob(path + "/*.csv")
initialfitlist = []
for filename in allFiles:
    file_ = open(filename, 'rb')
    mycsv = csv.reader(file_)
    mycsv = list(mycsv)
    initialfitlist.append(int(mycsv[0][1]))

"""
DETERMINES THE NUMBER OF GENERATIONS (ON AVERAGE) TO
GET TO A ZERO FITNESS, WITH N RUNS
"""

path =r'GA Results Data'
allFiles = glob.glob(path + "/*.csv")
lastgenlist = []
for filename in allFiles:
    file_ = open(filename, 'rb')
    mycsv = csv.reader(file_)
    mycsv = list(mycsv)
    lastgenlist.append(int(mycsv[-1][0]))

print "Number of GA Runs: %i" % (len(allFiles))
print "Average initial fitness of randomly generated DNA string: %f \n" % (np.mean(initialfitlist))
print "On average, it will take %f generations to reach the optimal protein." % (np.mean(lastgenlist))