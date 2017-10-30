#Weka Simulator
#Uses python to open a CSV file of attributes, and shuffles the 23rd attribute (rec)
#Saves CSV File, opens CSV File in Java.weka
#Runs Random Forrest Classifier on the data, outputs the accuracy and matrixes to a file.
#Rinse and Repeat

import os
import sys
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd 
import numpy as np
import argparse


parser = argparse.ArgumentParser(prog="WekaChrSim.py", usage="WekaChrSim.py <Chr>", 
                                 description='Performs a 100tree RF simulation and stores the output to calculate P values')

parser.add_argument('chrvar', type=str,
                    help='Chromosome? e.g. 2L, 2R, X, etc..')
parser.add_argument("-S", type=int, default=10,
                    help='Simulation Number, Default = 10')
args = parser.parse_args()
global chrvar
chrvar = args.chrvar
simulationTimes = args.S

#chrvar = str(raw_input("Chromosome? e.g. 2L, 2R, X, etc.."))
#os.system("export CLASSPATH=$CLASSPATH:/home/Andrew/weka/weka-3-6-12/weka.jar")
filename = "chr"+chrvar+"_Data.csv"

def mainfunction(times):
    counter = 0 
    datafile = pd.read_csv(filename)
    recData = datafile[datafile.columns[0]] # col 22= rec_class
    while counter < times:
        #print ' Shuffling CSV File... '
        datafile = datafile.drop('RG_class',axis=1)
        df = np.random.permutation(recData)
        datafile['RG_class'] = df
        os.system("rm "+chrvar+"_data_randomized.csv")
        datafile.to_csv(chrvar+"_data_randomized.csv",index=False)
        #print 'CSV Written, proceeding to RF simulation'
        
        os.system("java -Xmx1024m weka.classifiers.trees.RandomForest -t "+chrvar+"_data_randomized.csv -o -I 100 -K 0 -S 1 -i | grep -A 27 'Stratified' >> "+chrvar+"_sims.out")
        counter += 1
        if counter % 1000 == 0:
            print( str(counter) + ' simulations complete on '+chrvar)
    print' ==== All Simulations Complete! === \a\a\a \n'

mainfunction(simulationTimes)