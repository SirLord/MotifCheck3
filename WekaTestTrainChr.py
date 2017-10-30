#Weka Simulator
#Uses python to open a CSV file of attributes, and shuffles the 23rd attribute (rec)
#Saves CSV File, opens CSV File in Java.weka
#Runs Random Forrest Classifier on the data, outputs the accuracy and matrixes to a file.
#Rinse and Repeat

import os
import sys
import warnings
#with warnings.catch_warnings():
    #warnings.simplefilter("ignore")
    #import pandas as pd 
#import numpy as np

#os.system("export CLASSPATH=$CLASSPATH:/home/Andrew/weka/weka-3-6-12/weka.jar")
TrainList = ['allChr_motifs.csv.arff','Autosome_motif.csv.arff','Xmotif.csv.arff','2Lmotif.csv.arff','2Rmotif.csv.arff','3Lmotif.csv.arff','3Rmotif.csv.arff','4motif.csv.arff']
TestList = TrainList

for train in TrainList:
    for test in TestList:
		print(' testing '+test+' on model '+train)
		os.system("echo '"+train+" as train, "+test+" as test' >> RFMatrix.out")
		os.system("java -Xmx1024m weka.classifiers.trees.RandomForest -t "+train+" -T "+test+" -o -I 1000 -x 10 -S 777 -i -v| grep -A 2 'Weighted' >> RFMatrix2.out")

        
print' ==== All Simulations Complete! === \a\a\a \n'

#TrainList = ['chr2L_Data.arff','chr2R_Data.arff','chr3L_Data.arff','chr3R_Data.arff','X_Data.arff','chr4_Data.arff']
#TestList = ['chr2L_Data.arff','chr2R_Data.arff','chr3L_Data.arff','chr3R_Data.arff','X_Data.arff','chr4_Data.arff']

# "java -Xmx1024m weka.classifiers.trees.RandomForest -t 4motit.csv.arff -o -I 1000 -x 10 -S 777 -i -v| grep -A 2 'Weighted' 