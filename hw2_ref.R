######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
library("Biostrings",verbose=F,quietly=T)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw2_<your student ID>.R --input test.fasta --score pam250.txt --gap -10", call.=FALSE)
}


