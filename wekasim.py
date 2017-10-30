#Weka Simulator
#Uses python to open a CSV file of attributes, and shuffles the 23rd attribute (rec)
#Saves CSV File, opens CSV File in Java.weka
#Runs Random Forrest Classifier on the data, outputs the accuracy and matrixes to a file.
#Rinse and Repeat
#export CLASSPATH=$CLASSPATH:/home/Andrew/weka/weka-3-6-12/weka.jar

import os
import sys
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd 
import numpy as np


def mainfunction(times):
    counter = 0 
    datafile = pd.read_csv('rf_data_randomized.csv')
    recData = datafile[datafile.columns[22]] # col 22= rec_class
    while counter < times:
        #print ' Shuffling CSV File... '
        datafile = datafile.drop('RG_b50.results',axis=1)
        df = np.random.permutation(recData)
        datafile['RG_b50.results'] = df
        os.system("rm rf_data_randomized.csv")
        datafile.to_csv("rf_data_randomized.csv",index=False)
        #print 'CSV Written, proceeding to RF simulation'
        
        os.system("java -Xmx1024m weka.classifiers.trees.RandomForest -t newdata.csv -o -I 250 -K 0 -S 1 -i | grep -A 27 'Stratified' >> weka.out")
        counter += 1
        if counter % 100 == 0:
            print( str(counter) + ' simulations complete...')
    print' ==== All Simulations Complete! === \a\a\a \n'

mainfunction(25000)