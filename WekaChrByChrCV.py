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
#for test in TestList:
	    print(' testing '+train+' on model '+train)
	    os.system("echo '"+train+" as train, "+train+" as test (10-Fold CV) ' >> RFMatrix.out")
	    os.system("java -Xmx1024m weka.classifiers.trees.RandomForest -t "+train+" -o -I 1000 -x 10 -S 777 -i -v| grep -A 2 'Weighted' >> ChrByChr.out")

        
print' ==== All Simulations Complete! === \a\a\a \n'

#TrainList = ['chr2L_Data.arff','chr2R_Data.arff','chr3L_Data.arff','chr3R_Data.arff','X_Data.arff','chr4_Data.arff']
#TestList = ['chr2L_Data.arff','chr2R_Data.arff','chr3L_Data.arff','chr3R_Data.arff','X_Data.arff','chr4_Data.arff']

# "java -Xmx1024m weka.classifiers.trees.RandomForest -t 4motit.csv.arff -o -I 1000 -x 10 -S 777 -i -v| grep -A 2 'Weighted' 

"""
General options:

-h or -help
        Output help information.
-synopsis or -info
        Output synopsis for classifier (use in conjunction  with -h)
-t <name of training file>
        Sets training file.
-T <name of test file>
        Sets test file. If missing, a cross-validation will be performed
        on the training data.
-c <class index>
        Sets index of class attribute (default: last).
-x <number of folds>
        Sets number of folds for cross-validation (default: 10).
-no-cv
        Do not perform any cross validation.
-split-percentage <percentage>
        Sets the percentage for the train/test set split, e.g., 66.
-preserve-order
        Preserves the order in the percentage split.
-s <random number seed>
        Sets random number seed for cross-validation or percentage split
        (default: 1).
-m <name of file with cost matrix>
        Sets file with cost matrix.
-l <name of input file>
        Sets model input file. In case the filename ends with '.xml',
        a PMML file is loaded or, if that fails, options are loaded
        from the XML file.
-d <name of output file>
        Sets model output file. In case the filename ends with '.xml',
        only the options are saved to the XML file, not the model.
-v
        Outputs no statistics for training data.
-o
        Outputs statistics only, not the classifier.
-i
        Outputs detailed information-retrieval statistics for each class.
-k
        Outputs information-theoretic statistics.
-p <attribute range>
        Only outputs predictions for test instances (or the train
        instances if no test instances provided and -no-cv is used),
        along with attributes (0 for none).
-distribution
        Outputs the distribution instead of only the prediction
        in conjunction with the '-p' option (only nominal classes).
-r
        Only outputs cumulative margin distribution.
-xml filename | xml-string
        Retrieves the options from the XML-data instead of the command line.
-threshold-file <file>
        The file to save the threshold data to.
        The format is determined by the extensions, e.g., '.arff' for ARFF
        format or '.csv' for CSV.
-threshold-label <label>
        The class label to determine the threshold data for
        (default is the first label)

Options specific to weka.classifiers.trees.RandomForest:

-I <number of trees>
        Number of trees to build.
        (default 100)
-K <number of features>
        Number of features to consider (<1=int(log_2(#predictors)+1)).
        (default 0)
-S
        Seed for random number generator.
        (default 1)
-depth <num>
        The maximum depth of the trees, 0 for unlimited.
        (default 0)
-D
        If set, classifier is run in debug mode and
        may output additional info to the console
"""