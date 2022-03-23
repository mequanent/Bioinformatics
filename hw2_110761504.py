# hw2_110761504.py
# ©©© Mequanent Argaw Muluneh ©©©

# **** My Final local paths when writing and testing my code *******
# python_path C:\Users\USER\Desktop\BioInfoHW2\hw2_110761504.py
# Input_path: C:\Users\USER\Desktop\BioInfoHW2\test.fasta
# Score_path: C:\Users\USER\Desktop\BioInfoH2\pam250.txt

# Copy the following command to run the code, with the respective path changes
# python C:\Users\USER\Desktop\BioInfoHW2\hw2_110761504.py --input C:\Users\USER\Desktop\BioInfoHW2\test.fasta --score C:\Users\USER\Desktop\BioInfoHW2\pam250.txt --gap -10
# The arguments (--input test.fasta --score pam250.txt --gap -10) can be shuffled

import sys
import argparse # For parsing arguments when they are shuffled
import numpy as np
import pandas as pd
# pip install Bio   # Library installed to handle parsing with Biostrings in python
from Bio import SeqIO

# *****************************************************************************
# ******** Arguments' parsing section with shuffling possibility **************

if(len(sys.argv) < 4): # There should be four arguments including the .py file.
  print("You entered {} argument(s) which is less than required, please retry!".format(len(sys.argv)))
  sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "test.fasta")
parser.add_argument("--score", help = "pamx.txt")
parser.add_argument("--gap", help = "gap", type = int)

args = parser.parse_args()

input_path = args.input
score_path = args.score
gap = args.gap
# ***************** End of parsing *********************************************
# ******************************************************************************

# ****************************************************************************
# ******** Reading the the FASTA and parsing section *************************

# Reference for parsing the fasta in python is the original biopython library documentation
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

# The pam matrix must be indexed by the alphabets of the sequences to calculate the scores easily.

def preprocess_score(score):
  file = open(score, 'r')
  lines = file.readlines()[10:] # The first 10 lines are strings
  file.close()
  columns = []
  mut = []
  for line in lines:
    columns.append(line[0]) # Amino acid alphabets stored for row & column indexing
    line = line[1:].strip().split() # remove the letters from left and the spaces
    mut.append(line)
  mut = np.array(mut).astype(int) # change the string numbers into type int
  #print(mut.shape) # Uncomment to check the size of the score matrix
  size = len(mut) # save the size of the array for future usages to come below

  return pd.DataFrame(mut, index = columns, columns = columns)

score = preprocess_score(score_path) # score made a data frame to access its values with 
                  # alphabet indices
# ***************** End of score matrix preprocessing *************************************
# *****************************************************************************************


# *****************************************************************************************
# **************** A function to determine pairwise alignment score ***********************
def pairwise(x, y):
  pair_sum = 0
  for i in range(len(x)): # Considering that both x and y have the same length
    if x[i] == '-' or y[i] == '-':
      pair_sum = pair_sum + gap
    else:
      pair_sum = pair_sum + score.loc[x[i]][y[i]]
  return pair_sum
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

