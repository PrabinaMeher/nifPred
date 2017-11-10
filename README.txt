####### nifPred Help manual for executing the source code offline #######
#########################################################################

Here, we are providing the source code of nifPred which has been written in R.

The file "nifpred.r" can be run in a local system to execute nifPred by following the instructions given below:

1) Install R in your work station or Personal computer. 

2) Install the required packages of R i.e., "Biostrings", "BioSeqClass", "protr", "e1071" and "R2HTML". 

3) Unzipp the zipped file in your desired location after downloading the codes from Github.  

4) In the created directory, preapre the FASTA files (.fasta) containing the test sequences (protein sequences containing standard amino acid residues only).

   The name of the sequence file should be "test.fasta" (the name of the file is test and file extension is .fasta). 

5) Set the created directory as the current working directory in R, using the command setwd("PATH OF THE DIRECTORY") 

6) Source the R code (nifpred.r) using the command source("nifpred.r").

7) After the code is executed successfully, you will get two files such as "nifpred.txt" and "nifpred.html" in the directory. 
   The ".txt" file contains the raw result files and ".html" contains the processed result files. 
   For example the .html file will only display the results of predicted categories nif sequences, whereas ".txt" file contains results for all the sequences. 
   However, if no sequence is predicted into nif category, the ".html" file will also show the results for non-nif sequences.
   For detail about the processing of the results, you can refer our manuscript or the help pages of the web server at http://webapp.cabgrid.res.in/nifPred/help.html.
   If you still face any problem or want to give any feedback don't hesitate to contact us http://webapp.cabgrid.res.in/nifPred/contact.html 


########################### Thank YOU ####################################  ##########################################################################
