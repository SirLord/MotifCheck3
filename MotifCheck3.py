#!/usr/bin/env python
# LAST EDIT: 10-16-15 11:22
# MotifCheck3.py    
"""Module summary
This program is designed to count significant motif hits from a given matrix within a given D. mel sequence.
When given a matrix and sequence it will concatenate the Fasta sequence strings, shuffle them, and
calculate a FDR probability to use for that particular sequence and motif (also calculates a reverse
compliment motif matrix). Then, the program will scan overlapping windows of the input sequence and
assign a probability for each position. Probabilities below the FDR are removed, overlapping motifs 
dealt with, and finally you are given a summary of the significant motif counts. The positions and 
scores of detected motifs are returned. MotifCheck3 is designed to be used with MotifLandscape2.py, 
MotifLogo.py, RLassoMaker.py. Copyright Andrew Adrian March 2013
"""

import math
import sys
import random
import argparse
from time import *
from Bio import SeqIO

#Code for running commandline and parsing arguments

parser = argparse.ArgumentParser(prog="MotifCheck3.py", usage="MotifCheck3.py <targetsequence.fasta> <matrix.txt> <outfile.txt>", 
                                 description='Counts significant motif hits from a given input matrix file in a given input fasta sequence file.')

parser.add_argument('fastaFile', type=str,
                    help='input fasta file to search for motifs')
parser.add_argument('matrixFile', type=str,
                    help='input matrix file')
parser.add_argument('outputdirectory',type=str,
                    help='output text file')
parser.add_argument("-FDR", action="store_true", 
                    help='Turn on FDR Reporting. Default OFF (MEMORY INTENSIVE!)')
parser.add_argument("-FDR_Detect", action="store_true", 
                    help='Turn on FDR File Detection to speed repeat analyses.')
parser.add_argument("-FDR_Percent", type=int, default=1,
                    help='Change FDR Percentage. Default = 1')
parser.add_argument("-pseudo", type=float, default=1E-10,
                    help='Change matrix salting psuedocount. Default = 1E-10')
parser.add_argument("-exact", action="store_true", 
                    help='Require that motifs are exact words')
parser.add_argument("-threshold", type=float, default=1E-30,
                    help='Basal p below which potential motifs are ignored. Default = 1E-30 ')
parser.add_argument("-SupplyP", type=str, default=False,
                    help='Supply a FDR probability file')
parser.add_argument("-MS", action="store_true",
                    help='activate multishuffle. (Suffle Dinuclotide twice for FDR calcs)')
#parser.print_help()


args = parser.parse_args()

#get a succinct name for FDR File
global shortName
global fastaFile
fastaFile = args.fastaFile
if '/' in fastaFile:
    shortName = fastaFile[name[::-1].find('/')-1:-4].strip(".fasta")
else:
    shortName = fastaFile.strip(".fasta")
global matrixFile
matrixFile = args.matrixFile
global outputdirectory
outputdirectory = args.outputdirectory
global output_shortname
output_shortname = outputdirectory.strip(".txt")
global FDR_percent
FDR_percent = args.FDR_Percent
global pseudo
pseudo = args.pseudo
global exact
exact = args.exact
global THRESHOLD
THRESHOLD = args.threshold 
global SupplyP
SupplyP = args.SupplyP
#open outfile

outFile = open(outputdirectory,'w')
outFile.write("Process initiated " + str(asctime()) + '\n')
print("\t---MotifCheck.py---\nProcess initiated " + str(asctime()) )


#matrix file parser
def parseMatrixFile(inputfile):
    fileA = open(inputfile, "r")
    file = fileA.readlines()
    #pseudocount seems to work best at 1e-10 High sensitivity, low noise.
    #pseudo = float(raw_input("Pseudocount? (0,.001,.0001,.00001): "))
    print("Processing Matrix...")
    global motif_len
    global maxValue
    maxValue = 1.0
    index = 0
    del file[0]
    for item in file:
        file[index] = (item.strip(' \r\n')).split("  ")
        index += 1
    index = 0
    while index < len(file): #turn string into float
        x = 0
        while x < len(file[index]):
            file[index][x] = float(file[index][x]) + pseudo #add pseudocount to matrix
            x+=1
        index +=1
    fileA.close()
    motif_len = len(file)
    print("Matrix of size " + str(motif_len) + " detected")
    outFile.write("Motif Length: " +str(motif_len)+ "\n")
    
    for line in file:
        highval = 0.0
        for item in line:
            if item > highval:
                highval = item
        maxValue *= highval
    print(' maximum P is ' + str(maxValue))
    return file # Matrix as a list

#takes an input fasta file to search against.
def parseFastaFile(inputfile):
    global all_fasta
    all_fasta = ""
    
    #Check to see if we need parse fasta file
    FastaDetect = 0
    try:
        with open(shortName+'_PARSE.txt'):
            FastaDetect = 1
    except IOError:
        pass
   
    if FastaDetect > 0:
        print 'Previous FASTA files detected. Using: '+str(shortName)+'_PARSE.txt'
        fasta = open(shortName+'_PARSE.txt', 'r')
        all_fasta = fasta.readlines()[0]
        fasta.close()   
        
    elif FastaDetect == 0:
        print 'No Previous Fasta files detected.'
        ffile = open(inputfile, "r")
        print("Parsing FASTA input...")

        for seq_record in SeqIO.parse(ffile, "fasta"):
            all_fasta += str(seq_record.seq).upper()
            
        outFile.write("Scanned " +str(len(all_fasta))+ " positions\n")
        print("Saving FASTA File as: " +str(shortName)+'_PARSE.txt') #save fasta parsed
        fastaOutFile = open(shortName+'_PARSE.txt', 'w') 
        fastaOutFile.write(all_fasta)
        ffile.close()
        fastaOutFile.close()
        print("Parsing FASTA Complete...")
  
    return all_fasta
  
    
    
def singleShuffle(sequence):
    print("Nucleotide Shuffle in progress...")
    ShuffleOut =  ''.join(random.sample(sequence,len(sequence)))
    SingleShuffleFile = open(shortName+'SHUFFLE.txt', 'w')
    SingleShuffleFile.write(str(ShuffleOut))
    SingleShuffleFile.close()
    return ShuffleOut


def doubleShuffle(seq):
    print("Dinucleotide Shuffle in progress...")
    listSeq = []
    i = 1
    while i < len(seq):
        if i%2 == 0:
            listSeq.append(seq[i] + seq[i-1])
        i += 1
    DoubleShuffleFile = open(shortName+'DUB_SHUFFLE.txt', 'w')
    DoubleShuffleOut = "".join(random.sample(listSeq,len(listSeq)))
    DoubleShuffleFile.write(str(DoubleShuffleOut))
    DoubleShuffleFile.close()
    return DoubleShuffleOut

#Calculates an FDR threshold to test against. 
def FDR_Calculate(mat,seq):
    global FDR_value
    FDR_value = 0.00
    if exact: #Exact match = TRUE?
            print ("\texact = TRUE Requiring that probability == 1! Not recommended!")
            FDR_value = maxValue
            return
    #Check to see if we need to shuffle
    Detect = 0
    try:
        with open(shortName+'SHUFFLE.txt'):
            Detect = 1
    except IOError:
        print 'Except: Shuffle file not found'
        
    #set up multishuffle for rarer motifs. added 10-14-15 
    if args.MS:
        print ("Multishuffle enabled, please wait...")
        doubleShuffleSeq = doubleShuffle(seq) 
        doubleShuffleSeq += doubleShuffle(seq) #two shuffles!
        print("Calculating "+ str(FDR_percent)+ "% FDR probability...")
        singleShuffleSeq = False # so I can ignore this
   
    elif Detect > 0:
        print 'Previous Shuffle files detected. Using: '+str(shortName)+'SHUFFLE.txt'
        SingleShuffleFile = open(shortName+'SHUFFLE.txt', 'r')
        singleShuffleSeq = SingleShuffleFile.readline()
        SingleShuffleFile.close()
        DoubleShuffleFile = open(shortName+'DUB_SHUFFLE.txt', 'r')
        doubleShuffleSeq = DoubleShuffleFile.readline()
        DoubleShuffleFile.close()

    elif Detect == 0:
        print 'No Previous Shuffle files detected.'
        #Okay, gonna shuffle I guess :(
        singleShuffleSeq = singleShuffle(seq)
        doubleShuffleSeq = doubleShuffle(seq)
        print("Calculating "+ str(FDR_percent)+ "% FDR probability...")
        
   #always do double shuffle     
    FDR_List_Double = motifPAM(mat,doubleShuffleSeq)
    FDR_List_Double.sort()
    i = 0
    belowDouble = 0
    for item in FDR_List_Double:
        if item < THRESHOLD:
            belowDouble += 1
        else:
            i += 1
    FDR_List_Double = FDR_List_Double[belowDouble::]
    print (str(belowDouble) + " items below threshold removed from Double Shuffle") 
    checkNoDouble = int(((((len(FDR_List_Double)) / 100) * FDR_percent) / (motif_len) ) )    
    
    # Do single if MS = False
    if args.MS != True: 
        
        FDR_List_Single = motifPAM(mat,singleShuffleSeq)
        FDR_List_Single.sort()
        
        i = 0
        belowSingle = 0
        for item in FDR_List_Single:
            if item < THRESHOLD:
                belowSingle += 1
            else:
                i += 1
        FDR_List_Single = FDR_List_Single[belowSingle::]
        print (str(belowSingle) + " items below threshold removed from Single Shuffle")
        checkNoSingle = int(((((len(FDR_List_Single)) / 100) * FDR_percent) / (motif_len) ) ) #-1) # can change the 20 to math.fabs(math.log(THRESHOLD,10)
    else:
        FDR_List_Single = FDR_List_Double
        checkNoSingle = int(((((len(FDR_List_Double)) / 100) * FDR_percent) / (motif_len) ) )
            
    
    if args.FDR: #rename FDR File to something unique
        FDRFile = (open(str(output_shortname)+'_VerboseFDR.txt','w'))
    else:
        FDRFile = open(matrixFile + "FDR_" +str(shortName) + '.txt','w') #Make FDR File so I dont have to remake it
    #pick the conservative FDR value. Defaults to Double if MS on    
    if FDR_List_Single[-checkNoSingle] > FDR_List_Double[-checkNoDouble]:
        print "Item "+str(checkNoSingle)+ " Chosen for FDR"
        FDR_value = float(FDR_List_Single[-checkNoSingle])
        if FDR_value > THRESHOLD:
            alternate = False
        else:
            FDR_value = maxValue
            alternate = True
        if(alternate):
            print("\a!Alert! Minimum FDR p Threshold met. Using "+str(maxValue)+" instead...")
        else:
            print str(FDR_value) + " chosen for FDR (Single)"
        FDRFile.write(str(FDR_value)) 
        if args.FDR: #LOTS OF MEM NEEDED . avg file for a chromosome is 300mb
            print "Writing Verbose FDR File...."
            FDRFile.write('\n'+str(FDR_List_Single))    
    else:
        FDR_value = float(FDR_List_Double[-checkNoDouble])
        print "Item "+str(checkNoDouble)+ " Chosen for FDR"
        if FDR_value > THRESHOLD:
            alternate = False
            pass
        else:
            FDR_value = maxValue
            alternate = True
        if(alternate):
            print("\a !Alert! Minimum FDR p Threshold met. Using "+str(maxValue)+" instead...")
        else:
            print str(FDR_value) + " chosen for FDR (Dinucleotide)"      
        
        FDRFile.write(str(FDR_value)) 
        if args.FDR: #LOTS OF MEM NEEDED
            print "Writing Verbose FDR File...."
            FDRFile.write('\n'+str(FDR_List_Double))
    FDRFile.close()

def inspectMatrix(mat): #Allows easy inspecting of a matrix
    print( " A C G T \n" )
    for position in mat: 
        print position
        
def compliment(inputlist): #compliments a sublist for reverseComplimentMatrix()
    resultlist = [inputlist[3], inputlist[2], inputlist[1], inputlist[0]]
    return resultlist

def reverseComplimentMatrix(originalMat): #Provides a reverse compliment matrix.
    RCMatrix = []
    mat = originalMat[::-1]
    for subList in mat:
        RCMatrix.append(compliment(subList))
    return RCMatrix

#Main Function here!
def motifPAM(mat,seq):
    windowsize = len(mat)
    motifScore = 1
    RCmotifScore = 1
    motifPosition = 0 # beginning of sequence is 0
    seqScore = []
    global frameList 
    frameList = []
    bpIndex = 0
    reversible = mat
    # Make Variable RC_Matrix 
    RC_Matrix = reverseComplimentMatrix(reversible)
    
    #for each letter of a sequcnce
    #calculate the probability of it being a motif by multiplying the probs associated with each letter in a matrix
    #within a window. Then iterate forward and repeat.
    
    while bpIndex < (len(seq)-len(mat)): #goes through each bp and...
        # calculate probability for the window size len(mat)
        #print 'calculating start'
        start = bpIndex
        stop = bpIndex + len(mat) 
        windowSeq = seq[start:stop]
        motifScore = 1.0
        RCmotifScore = 1.0
        motifPosition = 0
        
        for letter in windowSeq: #go through each bp of the window
           # print 'entering for loop'
            if letter == "A":
                motifScore *= mat[motifPosition][0] # 0 is A
                motifPosition += 1
            elif letter == "C":
                motifScore *= mat[motifPosition][1] # 1 is C
                motifPosition += 1
            elif letter == "G":
                motifScore *= mat[motifPosition][2] # 2 is G
                motifPosition += 1
            elif letter == "T":
                motifScore *= mat[motifPosition][3] # 3 is T
                motifPosition += 1
            else:
                motifScore *= 1e-20
                
        motifPosition = 0
        for letter in windowSeq: #REVERSE COMPLIMENT CHECK
            #print 'entering for loop'
            if letter == "A":
                RCmotifScore *= RC_Matrix[motifPosition][0] # 0 is A
                motifPosition += 1
            elif letter == "C":
                RCmotifScore *= RC_Matrix[motifPosition][1] # 1 is C
                motifPosition += 1
            elif letter == "G":
                RCmotifScore *= RC_Matrix[motifPosition][2] # 2 is G
                motifPosition += 1
            elif letter == "T":
                RCmotifScore *= RC_Matrix[motifPosition][3] # 3 is T
                motifPosition += 1
            else:
                RCmotifScore *= 1e-20
                
                
        score   = float((motifScore   ))#*1000) 
        RCscore = float((RCmotifScore ))#*1000)

        if score > RCscore:
            seqScore.append(score)
            frameList.append("+")
        else:
            seqScore.append(RCscore)
            frameList.append("-")
   
        bpIndex += 1 #iterate forward
        
    return seqScore

# Now, let's identify how many motifs exist above the FDR and deal with overlaps (Rpt)
def interpretList(inputlist):
    print "Using FDR: " + str(FDR_value)
    index = 0 
    for item in inputlist: # makes values None if below FDR_value
        if item <= FDR_value:
            inputlist[index] = None
        index += 1  
        
    index = 0
    while ((index + motif_len) < len(inputlist)):
        if inputlist[index] != None:
            motifIndex = 1
            topVal = inputlist[index]
            topValIndex = index
            while (motifIndex < motif_len):
               # print 'inputlist[index+motifIndex]='+str(inputlist[(index+motifIndex)])
                if (inputlist[index+motifIndex] > topVal) and (inputlist[index+motifIndex] != ""):
                    topVal = inputlist[index+motifIndex]
                    inputlist[topValIndex] = 'Rpt'
                    topValIndex = index + motifIndex
                elif (inputlist[index+motifIndex] <= topVal) and (inputlist[index+motifIndex] != ""):
                    inputlist[index+motifIndex] = 'Rpt'
                    
                motifIndex += 1
                
            index += motif_len - 1
            
        index += 1
    
    rptCount = 0    
    count = 0  
    index = 0
    probDict = {} #{Position: [Value, Motif, Frame]}       #Frame value added may 27 2014 for making Logos
    for item in inputlist:
        if item == "Rpt":
            rptCount += 1
        elif item != None:
            count += 1
            probDict[index] = [item, all_fasta[index:index+motif_len], frameList[index]]
        index += 1
    print( "Found " +str(count)+ " motif hits with "+str(( count - rptCount/motif_len))+" overlapping motif hits\n")
    outFile.write("Found " +str(count)+ " motif hits with "+str(( count- rptCount/motif_len))+" overlapping motif hits\n")
    return probDict


all_fasta = parseFastaFile(fastaFile)
matrix = parseMatrixFile(matrixFile)

if args.FDR_Detect:
    try: 
        with open(matrixFile + 'FDR_' + str(shortName) + '.txt'):
            fileDetect = 1
            print 'Previous FDR File detected.'
    except IOError:
        print 'Exception: No FDR File!'
else: 
    print "\t FDR File Detection Disabled. Use -FDR_Detect"
    fileDetect = 0
    
if SupplyP:
    try: 
        with open(SupplyP):
            fileDetect = 1
            print ' FDR Probability Supplied...'
            FDRFile = open(SupplyP, 'r')
            FDR_value = float(FDRFile.readline().strip('\n')) 
            FDRFile.close()
            motifOut = motifPAM(matrix, all_fasta)
            probDict = interpretList(motifOut)
    except IOError:
                print 'Exception: '+str(SupplyP)+' is invalid!'
                print str(FDR_value) 
                
   
if fileDetect > 0 and args.FDR_Detect:
    print ' Using: ' + matrixFile + 'FDR_' + str(shortName) + '.txt'
    FDRFile = open(matrixFile + 'FDR_' + str(shortName) + '.txt', 'r')
    FDR_value = float(FDRFile.readline().strip('\n')) 
    FDRFile.close()
    motifOut = motifPAM(matrix, all_fasta)
    probDict = interpretList(motifOut)
    
elif fileDetect == 0:
    print 'No previous FDR file detected. Constructing FDR File'
    FDR_Calculate(matrix, all_fasta)
    motifOut = motifPAM(matrix, all_fasta)
    probDict = interpretList(motifOut)


#Recordkeeping

for item in matrix:
    outFile.write(str(item) + "\n\n")
outFile.write("FDR Threshold: " + str(FDR_value) + "\n")
#outFile.write(str(motifOut) + "\n") #Commented Out; For testing only
#print (str(motifOut) + "\n")
outFile.write("\n\n Position: Probabilities, Motif \n" + str(probDict))
outFile.write("\n\nProcess completed " + str(asctime()))
outFile.close()


####Additional Stuff For Testing
#testSeq = "TCATACGTAGATCACATCATATGATGGTACACACATGATGATCATCATACGATAGTGTGACCCGCTGTACATCATCAGGAGTGTGACACA" 
##matrix is in the format A C G T
#testMat =([[0,.9,0,.1], [.7,.1,.05,.15],[.0,.0,0,1],[.0,1,.05,.05],[.9,.0,.1,.0], [.0,0,.0,1]])

