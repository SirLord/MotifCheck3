
#Andrew Adrian. MotifLandscape2.py
''' Takes output from MotifCheck, and provides a graphical representation and computes frequency.
    Does a smoothing with 1/2 window overlap (can be disabled). Only works for Dm3 genome. Will need modification for others
    Needs R and GGplot2 package to be installed'''
#Last update: 11-12-14 20:27p
# Log: Fixed a bug that doubled the actual outcome by not smoothing correctly
# Log: Got rid of the prefix option.
# Added CAT Chromosome 
import os
import sys
import argparse

#Code to run commandline


parser = argparse.ArgumentParser(prog="MotifLandscape.py", usage="MotifLandscape.py <FromMotifCheck.txt>", 
                                 description='Tallys motif hits and generates basic landscapes for downstream analysis.')

parser.add_argument('motifOut', type=str,
                    help='Input MotifCheck file to tally motifs')

parser.add_argument("-chr", type=str, default=False,
                    help='Specify chromosome.')

parser.add_argument("-windowsize", type=int, default=100000,
                    help='Alter windowsize from 100kb')

parser.add_argument("-Prefix", type=str, default=False,
                    help='Rename the outfile to have a different prefix')

parser.add_argument("-nosmooth", action="store_false", default=True,
                    help='disables 1/2 window smoothing.')

args = parser.parse_args()

motifOut = args.motifOut

if args.Prefix:
    outGraphPrefix = args.Prefix
else:
    outGraphPrefix = motifOut.strip(".txt")
    
if args.chr:
    chromosome = args.chr
else:
    chromosome = outGraphPrefix
window = args.windowsize

if args.nosmooth:
    smoothDivide = 2.0
else:
    smoothDivide = 1.0
    
#pick right size for D.mel chr
if "2L" in chromosome: 
    chrLen = 23011544
    print("Chr 2L detected...")
elif "2R" in chromosome:
    chrLen = 21146708
    print("Chr 2R detected...")
elif "3L" in chromosome:
    chrLen = 24543557
    print("Chr 3L detected...")
elif "3R" in chromosome:
    chrLen = 27905053
    print("Chr 3R detected...")
elif "chr4" in chromosome:
    chrLen = 1378900
elif "X" in chromosome:
    chrLen = 22422827
    print("Chr X detected...")
elif "CAT" in chromosome:
    chrLen = 119029700
    print("Chr CAT detected...")
#elif "YOUR CHROMOSOME" in chromosome:
    #chrLen = LENGTH OF SEQUENCE
    #print("YOUR CHROMOSOME detected...")    
else:
    raise TypeError("Unrecognized Chromosome. Use -chr with 2L,2R,3L,3R,chr4,chrX,CAT")

#infile named TestOut
test = open(motifOut, 'r')
testData = test.readlines()
test.close()
#Find the item in the list of testData with a dictionary
index = 0
while index < len(testData):
    if "{" in testData[index]:
        print "data located"
        break #breaks the while loop, using present index
    elif index == len(testData):
        raise IOError('Incompatable data type')
    else:
        index += 1
      
data = eval(testData[index].strip())
dataKeys = data.keys()
dataKeys.append(chrLen)
dataKeys.sort()

#initialize output csv
outfile = open("TEMP.csv", 'w')
outfile.write("windowEnd,Count\n")
cursor = int(window/smoothDivide)
itemcount = 0
posCountDict = {}
intervals = cursor
#while intervals < chrLen:
    #posCountDict[intervals] = 0
    #intervals += cursor
#print posCountDict
i = 0
while i < len(dataKeys):
    if int(dataKeys[i]) <= int(cursor):
        #print str(dataKeys[i]) +" below the cursor value of: "+str(cursor)
        #print "Scanning region "+str(cursor-window) + " to "+ str(cursor)
        itemcount += 1
        i += 1
        #print str(item) +" is less than "+str(cursor)
        #print "Itemcount = " +str(itemcount)
    else:
        #print "item greater than cursor, storing previous itemcount"
        #outfile.write(str(cursor) +","+str(itemcount) +"\n") #Uncomment to restore non-Overlap 
        posCountDict[int(cursor)] = int(itemcount)
        itemcount = 0
        cursor += int(window/smoothDivide)
        i = i #keep same i
        
    #posCountDict correct at this point.
    
if cursor < chrLen:  #Need to add extra values to the end of chrLen
    windowsExtra = int((chrLen-cursor)/window)
    while windowsExtra >= 0:
        #outfile.write(str(cursor) +","+ '0' +"\n")
        windowsExtra -= 1
        
keysSort = posCountDict.keys()
keysSort.sort() #Sorted keys in a list for calling in order.

index = 0
while index < (len(keysSort)-1):
    if smoothDivide == 2:
        smoothVal = ((posCountDict[keysSort[index]] + posCountDict[keysSort[index + 1]])/smoothDivide)
        outfile.write(str(int((keysSort[index]+keysSort[index+1])/smoothDivide)) +","+str(smoothVal) +"\n")
        index += 1
    else:
        CountVal = ((posCountDict[keysSort[index]]))
        outfile.write(str(int((keysSort[index]))) +","+str(CountVal) +"\n")
        index += 1        
outfile.close()

    
#Form an R Script, and load the graph we just made

R_Outfile = open("tempFile.R", 'w')
# change to the new directory
#R_Outfile.write("setwd(\"/home/Andrew/Desktop\")\n")
# load the necessary libraries
R_Outfile.write("library(ggplot2)\n")
# set the output file
R_Outfile.write("sink(\"temp.out\")\n")
# load the dataset
R_Outfile.write("data <- read.csv(\"TEMP.csv\", header=TRUE)\n")
#create the plot
R_Outfile.write("rplot <- ggplot(data, aes(x=windowEnd,y=Count)) + geom_bar(stat=\"identity\")\n")
R_Outfile.write("rplot <- rplot + ggtitle(\"" + str(outGraphPrefix) + "\") + theme_bw()\n")
R_Outfile.write("pdf(\"" + str(outGraphPrefix) + "_Plot.pdf\")\n")
R_Outfile.write("rplot\n")
R_Outfile.write("dev.off()\n")
# close the output file
R_Outfile.write("sink()\n")
R_Outfile.close()



#Run the R script from command line
os.system("Rscript tempFile.R")

#Clean up files
var = 'y' #deactivating raw input
if var == "y":
    os.system("mv TEMP.csv " + outGraphPrefix + ".csv" )
    print ("Data saved as "+ str(outGraphPrefix) + ".csv" )
else:
    print("Cleaning up temporary files...")
    os.system("rm TEMP.csv")
os.system("rm tempFile.R")

try:
    with open(str(outGraphPrefix) + "_Plot.pdf"):
        print("Process completed. Graph saved as: "+str(outGraphPrefix)+"_Plot.pdf")
except IOError:
    print 'Process Failed. Check Input Data'


