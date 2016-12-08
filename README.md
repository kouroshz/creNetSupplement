# creNetScripts
Example runs and command line scripts for facilitating creNET runs.

After installing the package creNet, use the following command line to run the program.

##USAGE##
For Usage:
```{R}
./runCreNet.R -h
```
For test run:
```{R}
./runCreNet.R -n creNet -m test -e ./KB/STRING_BEL_IREF.ents -r ./KB/STRING_BEL_IREF.rels -d ./test_data/UC1.txt -t ./test_data/UC2.txt -o ./results.txt -z ./significant_genes.txt
```
###Arguments:###
1) -n creNet: run the program in creNet mode.
2) -m test: run the program for independent train and test
3) -i 3: run the program using 3 iterations to get mean and standard deviations for accuracy measures.
4) -e ./KB/STRING_BEL_IREF.ents: Path to the KB entitiy file
5) -r ./KB/STRING_BEL_IREF.rels: Path to the KB relation file
6) -d ./test_data/UC1.txt: Path to the training data
7) -t ./test_data/UC2.txt: Path to the test data
8) -o ./results.txt: output file of accuracy results
9) -z ./significant_genes.txt: Significant predictors picked by the model.

NOTE: The training and testing data are assumed to be tab delimited gene expression data (either normalized microarray or voom normalized RNAseq). The first column should be the binary respose (0 or 1) and labeled as "y" and the rest of the columns should be gene expressions with entrez ids as column names.

