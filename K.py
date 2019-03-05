#!/usr/bin/env python
import argparse,subprocess,os,sys
from datetime import datetime
from subprocess import call
import shlex
from subprocess import Popen
from subprocess import PIPE
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import csv
import timeit
from Bio import SeqIO
import resource
import shutil

def FileLength(file):
    with open(file) as fin:
        for i, l in enumerate(fin):
            pass
    return i + 1

def collapse(species): 
    dup = set()
    return [i for i in species if i not in dup and not dup.add(i)]


sys.stdout.flush()
parser = argparse.ArgumentParser(description='Metagenomic Analysis of bacterial, viral and archael communities at the species level from NON-human associated habitat')
parser.add_argument('-o',action='store',help='Output directory')
parser.add_argument('-t',action='store',help='Number of threads [1]',type=int,default=1)
parser.add_argument('-infile1',action='store',help='Path to gzip fastq 1')
parser.add_argument('-infile2',action='store',help='Path to gzip fastq 2')
parser.add_argument('-outputprefix',action='store',help='Output prefix')
parser.add_argument('-KallistoMarkerIdx',action='store',help='Kallisto markers index',default='./data/Marker_kallistoIdx')
parser.add_argument('-FullGenome',action='store',help='Full genome fasta',default='./data/MicrobialCompleteGenome.fa')
parser.add_argument('-SpeciesList',action='store',help='Path to species list',default='./data/Microbial_ID')

args = parser.parse_args()

######################################################################
# Defining arguments
OUTPUTDIR = args.o
N = args.t
FILE1 = args.infile1
FILE2 = args.infile2
OUTPUTPREFIX = args.outputprefix
MIdx = args.KallistoMarkerIdx
FG = args.FullGenome
Species = args.SpeciesList

#######################################################################
# Check if output directory exists
if not os.path.exists(OUTPUTDIR):
    os.makedirs(OUTPUTDIR)
if not os.path.exists(OUTPUTDIR+"/"+OUTPUTPREFIX):
    os.makedirs(OUTPUTDIR+"/"+OUTPUTPREFIX)
#ensure_dir(OUTPUTDIR)
########################################################################
# Calculate total reads in raw fastq

print datetime.now(), "..... Counting total reads"
foutstat=open(os.path.join(OUTPUTDIR,OUTPUTPREFIX+"_StatFile"),"ab")

foutstat.write("Total reads in raw fastq" +"\t"+ str(FileLength(FILE1)/4)+"\n")

#######################################################################
# run kallisto on markers
Presence = OUTPUTDIR + "/"  + OUTPUTPREFIX + "marker"
os.makedirs(Presence)
print datetime.now(), "..... Running Kallisto (markers)"
marker = subprocess.Popen(("kallisto quant -i %s -o %s -t %s %s %s" % (MIdx , Presence , N , FILE1 , FILE2)), shell = True)
marker.communicate()

########################################################################
# Getting list of markers genome from kallisto


markerlist=[]
number = 0
fin1 = open(Presence+"/"+"abundance.tsv",'r')
for line in fin1:
    line = line.strip('\n').strip('\r')
    if number>0:
        line = line.split('\t')
        if float(line[3]) > 0:
            markerlist.append(line[0])
    number+=1
  
fin1.close()

# Getting species name from each ID number (for full genome and marker genome)
List = open(Species, 'r')
Fulllist=dict()
Markerlist={}
Fulllist2={}
for Line in List:
    Line = Line.strip('\n').strip('\r')
    Line = Line.split('\t')
    FullName=Line[0]
    MarkerName = Line[1]
    SpeciesName = Line[2]
    Fulllist2[FullName]=[SpeciesName]
    Markerlist[MarkerName] = FullName
List.close()

# Collapse duplicates in a list

UniqueSpecies=[]
# Get full genome ID based on species name from unique species marker
for uniquespecieslist in markerlist:

    UniqueSpecies.append(Markerlist[uniquespecieslist]) 

UniqueFullName=collapse(UniqueSpecies)

#######################################################################
# Get detected species fasta from full genome fasta file

fasta_file = FG  # Input fasta file
result_file = OUTPUTDIR+"/"+OUTPUTPREFIX+ "detected.fa" # Output fasta file


extraction = set(UniqueFullName)
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in extraction:
            SeqIO.write([seq], f, "fasta") 

print datetime.now(), "..... Building Kallisto index (full)"


#######################################################################
# Build kallisto index based on extracted full fasta (FF)

index = subprocess.Popen(("kallisto index -i %s %s"%(OUTPUTDIR+"/"+"detectedFFKI" , result_file)), shell = True)
index.communicate()


# Run full kallisto
FIdx=OUTPUTDIR+"/"+"detectedFFKI"

print datetime.now(), "..... Running Kallisto (full)"
full = subprocess.Popen(("kallisto quant -i %s -o %s -t %s %s %s" \
  % (FIdx , OUTPUTDIR+"/"+OUTPUTPREFIX , N , FILE1 , FILE2)), shell = True)
full.communicate()


#######################################################################
# Getting counts of full genome pseudoalignment 
fulldict={}
LineNumber = 0
fin2 = open(OUTPUTDIR+"/"+OUTPUTPREFIX+"/abundance.tsv")
for line in fin2:
    line = line.strip('\n').strip('\r')
    if LineNumber>0:
        line = line.split('\t')
        if float(line[3]) > 0:
            fulldict[line[0]]=line[3]
    LineNumber+=1

fin2.close()


# Find the species name associated with the ID number
FinalDict={}
for FullID in fulldict:
    if FullID in Fulllist2:
        if tuple(Fulllist2[FullID]) in FinalDict:
            FinalDict[tuple(Fulllist2[FullID])]+=float(fulldict[FullID])
        else:
            FinalDict[tuple(Fulllist2[FullID])]=float(fulldict[FullID])

#print fulldict
fout=open(os.path.join(OUTPUTDIR,OUTPUTPREFIX+"_FinalCount.txt"),"w")
for x in FinalDict:
    fout.write('%s\t%s\n' % (str(x[0]),FinalDict[x]))
fout.close()
foutstat.close()

#########################################################################
# removing intermediate files
os.remove(FIdx)
os.remove(result_file)
shutil.rmtree(Presence)