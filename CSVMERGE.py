#CSVMERGE
import os
import sys

#Code to run commandline
if (len(sys.argv) < 4):
    print 'Usage: CSVMERGE.py <OUTNAME.csv> <file1.csv> <file2.csv> <...>'
    sys.exit(1)
numfiles = len(sys.argv) - 1 
OUTFILE = open(str(sys.argv[1]),"w")

#numfiles = 2 - 1  
#OUTFILE = open("testcsvOUT.csv",'w')
#argv = ["testing_.csv", "testing_2.csv"]
outDict = {1:""} # Line:Values separated by commas
i= 2
#i = 0
while i <= numfiles:
    csvFile = open(str(sys.argv[i]),"r")
   # csvFile = open(argv[i],"r")
    print 'reading file ' + str(sys.argv[i])
    csvList = csvFile.readlines() #['windowEnd,Count\n', '75000,140\n']
    line = 1
    for entry in csvList:
        print "reading line " + str(line)
        if line not in outDict.keys():
            outDict[line] = "" #intialize dictionary entry 
        text = entry.strip("\n").split(",") #['windowEnd', 'Count']
        if text == "Count":
            text = sys.argv[i][0:sys.argv[i].find(".csv")]
        else:
            outDict[line] =  outDict[line] + text[1] + ","
        line += 1
    #print str(sys.argv[i]) + " added to file"
    csvFile.close()
    i+= 1
    line = 1
#Finish the dictionary 
for entry in outDict:
    outDict[entry] = outDict[entry][:-1] + "\n"
    OUTFILE.write(outDict[entry])
OUTFILE.close()
