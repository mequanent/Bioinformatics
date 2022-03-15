library("Biostrings",verbose=F,quietly=T)

# read fasta file
ff <- readAAStringSet("test.fasta")
seq_name = names(ff)
sequence = paste(ff)
print(sequence[1])
print(sequence[2])



