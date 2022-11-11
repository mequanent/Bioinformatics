# sumOfPairsAlignmentScore.py
# ©©© Mequanent Argaw Muluneh ©©©

# Copy the following command to run the code, with the respective path changes
# python C:\Users\USER\Desktop\BioInfoHW2\hw2_110761504.py --input C:\Users\USER\Desktop\BioInfoHW2\test.fasta --score C:\Users\USER\Desktop\BioInfoHW2\pam250.txt --gopen -10 --gextend -2
# The arguments (--input test.fasta --score pam250.txt --gopen -10 --gextend -2) can be shuffled

import sys
import argparse # For parsing arguments when they are shuffled
import numpy as np
import pandas as pd

# Library installed to handle parsing with Biostrings in python
# Library installed to handle parsing with Biostrings in python
#To install the Bio library into my IDLE python editor I used the following two steps.
#    1. cd C:\Users\USER\AppData\Local\Programs\Python\Python39\Scripts
#    2. pip install Bio
#pip install Bio  

from Bio import SeqIO

# *****************************************************************************
# ******** Arguments' parsing section with shuffling possibility **************

if(len(sys.argv) < 4): # There should be four arguments including the .py file.
  print("You entered {} argument(s) which is less than required, please retry!".format(len(sys.argv)))
  sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "test.fasta")
parser.add_argument("--score", help = "pamx.txt")
parser.add_argument("--gopen", help = "open gap", type = int)
parser.add_argument("--gextend", help = "extension gap", type = int)

args = parser.parse_args()

input_path = args.input
score_path = args.score
gop = args.gopen
gex = args.gextend
# ***************** End of parsing *********************************************
# ******************************************************************************

# ****************************************************************************
# ******** Reading the the FASTA and parsing section *************************

# Reference for parsing fasta in python is the original biopython library documentation
# (https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc11 section 2.4.1) 

seq_name = []
sequence = []

for seq_record in SeqIO.parse(input_path, "fasta"):
  seq_name.append(seq_record.id)
  sequence.append(seq_record.seq)  
print("+--------------++----------------++---------------------------+")
print("The first two sequences as indicated in the read_fasta.R are: ")
print(sequence[1])
print(sequence[2])
print("+--------------++----------------++---------------------------+")
# ***************** End of FASTA preprocessing *************************************
# **********************************************************************************

# **********************************************************************************
# **************** Preprocessing section for the pamx score ************************

# The pam matrix better be indexed by the alphabets of the sequences to calculate the scores easily.
# Hence, score will be returned in a data frame format 
def preprocess_score(score):
  file = open(score, 'r')
  lines = file.readlines()
  file.close()
  columns = []
  mut = []
  for line in lines:
    if line.startswith("#") or line.startswith(" "):
      continue
    columns.append(line[0]) # Amino acid alphabets stored for row & column indexing
    line = line[1:].strip().split() # remove the letters from left and the spaces
    mut.append(line)
  mut = np.array(mut).astype(int) # change the string numbers into type int
  #print(mut.shape) # Uncomment to check the size of the score matrix
  
  return pd.DataFrame(mut, index = columns, columns = columns)

score = preprocess_score(score_path) # score made a data frame to access its values with 
                  # alphabet indices
# ***************** End of score matrix preprocessing *************************************
# *****************************************************************************************

# *****************************************************************************************
# **************** A function to determine pairwise alignment score ***********************

def pairwise(x, y):
  pair_sum = 0
  x_prev = '' # x_prev and y_prev are used to store immediate previous alphabets
  y_prev = '' # so that extension gaps will be identified
  for i in range(len(x)): # Considering that both x and y have the same length
    x_prev = x_prev == '-' and x[i] == '-' # status of extension gap in the first paired sequence
    y_prev = y_prev == '-' and y[i] == '-' # status of extension gap in second paired sequence
    curr = x[i] == '-' or y[i] == '-' # status of open gap in the paired sequence
    if x_prev or y_prev: # extension gap value will be taken if holds
      pair_sum = pair_sum + gex
    elif curr: # open gap value will be taken if holds
      pair_sum = pair_sum + gop
    else: # score matrix value will be taken otherwise
      pair_sum = pair_sum + score.loc[x[i]][y[i]] 
    x_prev = x[i]
    y_prev = y[i]
  return pair_sum

pairwise(sequence[0], sequence[1])
# ***************** End of pairwise alignment scoring function ****************************
# *****************************************************************************************

# pairwise(sequence[0], sequence[1]) # Uncomment to verify the function's validity

# *****************************************************************************************
# **************** A function to determine the total sum of pairs score *******************
def sumOfPairs(seq):
  sop = 0 # to store sum of pairs value
  for i in range(len(seq)-1): # i iterates as long as there is j in seq as its next sequence
    for j in range(i+1, len(seq)): # j will be the next seq index to be compared with seq[i]
      sop = sop + pairwise(seq[i], seq[j])
  return sop
# ***************** End of the total sum of pairs score function **************************
# *****************************************************************************************


# *****************************************************************************************
# ********************************* Displaying Section ************************************
print("+--------------+----------------+")
print("Given:")
print("   + Number of sequences: ", len(sequence))
print("   + Score matrix: pam{}".format(score_path[-10:-4].partition('m')[2])) 
print("+--------------+----------------+")
print("sum-of-pairs (SoP) score: {} ".format(sumOfPairs(sequence)))
print("+--------------+----------------+")
# ************************* End of All Things *********************************************
# *****************************************************************************************

