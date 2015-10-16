#Get Motif Alignment
#Re-Generate Position Frequency Matrix
#Display Sequence Logo for Motif
'''REQUIRES R and seqLogo package!'''

import sys
import os
import math
import argparse

##Code to run commandline
#Code for running commandline and parsing arguments

parser = argparse.ArgumentParser(prog="MotifAlign.py", usage="MotifAlign.py <fromMotifCheck2.txt> -options", 
                                 description='Counts significant motif hits from a given input matrix file in a given input fasta sequence file.')

parser.add_argument('inputfile', type=str,
                    help='input file to generate motif logo')
parser.add_argument("-PFM", action="store_true", 
                    help='Return the PFM for the Logo')
parser.add_argument("-RC", action="store_true", 
                    help='Returns a Reverse Compliment Logo')
args = parser.parse_args()

inputfile = args.inputfile
PFM = args.PFM
RC = args.RC

#Read the posistions from inputfile (from MotifCheck2.py)
def readPositions(inputfile):
    
    #infile named TestOut
    test = open(inputfile, 'r')
    testData = test.readlines()
    test.close()
    CurrentMotif = inputfile[0:inputfile.find("_")]
    #Find the item in the list of testData with a dictionary
    index = 0
    while index < len(testData):
        if "{" in testData[index]:
            print "data located for  " + str(inputfile)
            break #breaks the while loop, using present index
        elif index == len(testData):
            raise IOError('Incompatable data file \a')
        else:
            index += 1
    data = eval(testData[index].strip())
    return data

data = readPositions(inputfile)

dataList = []
AvgProbList = []

global motifLen 
motifLen = len(data[data.keys()[0]][1])

for item in data.keys():
    AvgProbList.append(data[item][0])
    if data[item][2] == "+" and len(data[item][1]) == motifLen :
       # outfile.write(data[item][1]+"\n")
        dataList.append(data[item][1])
        
    elif data[item][2] == "-":
        Reverse = data[item][1][::-1]
        reverseCompliment = ""
        for letter in Reverse:
            if letter == "A":
                reverseCompliment += "T"
            elif letter == "T":
                reverseCompliment += "A"
            elif letter == "G":
                reverseCompliment += "C"
            elif letter == "C":
                reverseCompliment += "G"
        #outfile.write(reverseCompliment+"\n")
        dataList.append(reverseCompliment)

#Make a pwm 
def makePWM(dataList): 
    # c("ACGT")
    outmatrix = [] # [[0,0,0,0,],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for i in range(0,motifLen):
        outmatrix.append([0,0,0,0])
    for motif in dataList:
        index = 0
        for letter in motif:
            if letter == "A":
                outmatrix[index][0] += 1
            elif letter == "T":
                outmatrix[index][3] += 1
            elif letter == "G":
                outmatrix[index][2] += 1
            elif letter == "C":
                outmatrix[index][1] += 1
            index += 1
            
    return outmatrix
    
totals = makePWM(dataList)
posIndex = 0
for position in totals:
    total = 0.0
    index = 0
    for base in position:
        total += base
    for base in position:
        totals[posIndex][index] = round(base/total,4)
        index += 1
    posIndex += 1
x = totals
rString = ""
for item in x:
    for number in item:
        rString += str(number) + ","
if(PFM):
    outLine = ""
    for line in totals:
        for position in line:
            outLine += str(position) + " "
        outLine += '\n'
    print ">"+inputfile.strip(".txt")
    print outLine + "\n"

tempfile = open("temp.R","w")
if RC:
    tempfile.write("library(PWMEnrich)\n")
    
tempfile.write("library(seqLogo)\n")
tempfile.write("pdf(\"Logo_"+inputfile.strip(".txt")+".pdf\")\n")
if RC:
    tempfile.write("seqLogo(makePWM(reverseComplement(matrix( c(" +rString[0:-1]+ "), nrow=4, ncol=" +str(motifLen)+ "),alphabet=\"DNA\")))\n")
else:
    tempfile.write("seqLogo(makePWM(matrix( c(" +rString[0:-1]+ "), nrow=4, ncol=" +str(motifLen)+ "),alphabet=\"DNA\"))\n")
tempfile.write("dev.off()")
tempfile.write("\n")
tempfile.close()
os.system("Rscript temp.R")
#os.system("rm temp.R")
#os.system("open Logo."+inputfile.strip(".txt")+".pdf")    
print ("Sequence Logo saved as: Logo_"+inputfile.strip(".txt")+".pdf")
leng = len(AvgProbList)
    
print "Avg probability of "+str(leng)+" entries is "+ str(float(sum(AvgProbList)/leng))

#outfile.close()


