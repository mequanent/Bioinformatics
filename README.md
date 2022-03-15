# hw2
<your name + student ID>
## Description

* Write R or Python script to calculate the sum-of-pair score (SoP) of the multiple sequence alignment.
* Creating your own script with your student ID, ie. hw2_105753026.R.
* In this program, library Biostrings is only used to parse fasta file.

## File

* hw2_ref.R: You can start from this reference code, and try to write your own comment in English
* pam100.txt
* pam250.txt
* test.fasta

## Parameters

* input: fasta file (ex. test.fasta)
* score: score file (ex. pam250.txt)
* gap: gap score

## Command

Executing your code with the following command.

```R
Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gap -10
```

```Python
python3 hw2_studentID.R --input test.fasta --score pam250.txt --gap -10
```
The answer is 999. You should print it on the screen.

## Evaluation

10 testing data

```R
Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gap -10
Rscript hw2_studentID.R --input test2.fasta --score pam100.txt --gap -8
Rscript hw2_studentID.R --input data/test3.fasta --score pam/pam1.txt --gap -5
```

```Python
python3 hw2_studentID.R --input test.fasta --score pam250.txt --gap -10
python3 hw2_studentID.R --input test2.fasta --score pam100.txt --gap -8
python3 hw2_studentID.R --input data/test3.fasta --score pam/pam1.txt --gap -5
```


Correct answer gets 10 points of each testing data.

### Penalty

* High code similarity to others: YOUR SCORE = 0

