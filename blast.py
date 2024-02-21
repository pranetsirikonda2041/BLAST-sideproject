from Bio import SeqIO
import sys
from sys import argv
from kaestimate import *
x=0
i=0
j=0
names = []
length = []
start = time.time()
if len(argv) != 3:
    print("Usage: python blast.py file.fasta file.fasta")
    sys.exit(1)
while x<3:
    for record in SeqIO.parse(argv[x], "fasta"):
        print("%s %i" % (record.id, len(record)))
        length.append(len(record))
        names.append(record.id)
        print(length)
    x+=1


def e_valcalc():
# K × m × n × e^λS =  E Value
# K and λ are the Karlin-Altschul parameters. They can be estimated from large sets of random sequence alignments.
#The λ parameter normalizes the alignment score, while the K parameter scales the E-value based on the database
#and sequence lengths.
#S is the alignment score. It is calculated based on the selected scoring matrix and the given sequence alignment.
#The score reflects the sum of substitution and gap scores for the aligned residues.
#n is the query length(base pairs)
#m is the length of the database (i.e., the sum of all the lengths of all the sequences in the database).
    m=0
    n = length[0]
    print(m)
    for lengths in length:
        m = m + lengths
    print(n)


e_valcalc()
def lambdacalc():
    #alphabet
    alph = "ACGT"

    #Scoring Scheme
    match=1
    mismatch=-1
    gapopen=-1
    gapextend=-1

    #Computational Effort
    threads=48
    maxtime=60*60*3 #60*60*13

    ####################################
    #END PARAMETERS
    ####################################
    distr=getBestScoreDistribution(m,n,alph,lambda a,b: getLocalAlignmentScore(a,b,match,mismatch,gapopen,gapextend),threads,maxtime)
        