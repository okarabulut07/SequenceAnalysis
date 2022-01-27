# Load libraries needed for this analysis
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Define a path for the output PDF to be written
outfile="/home/onur/workspace/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf"

# Load phenotype data
pheno_data = read.csv("adv_infection_vs_virus_positive.csv")

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(pheno_data$type)+1,las=2,ylab='log2(FPKM+1)')

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[5790]

# Display the gene name for a single row of data 
ballgown::geneNames(bg)[5790]



boxplot(fpkm[5790,] ~ pheno_data$type, border=c(2,3), main=paste(ballgown::geneNames(bg)[5790],' : ', ballgown::transcriptNames(bg)[5790]),pch=19, xlab="Type", ylab='log2(FPKM+1)')


points(fpkm[5790,] ~ jitter(c(2,2,2,1,1,1)), col=c(2,2,2,1,1,1)+1, pch=16)



plotTranscripts(ballgown::geneIDs(bg)[5790], bg, main=c('RPS8 in all virus_positive samples'), sample=c('SRR12662521_expression', 'SRR12662530_expression', 'SRR12662532_expression'), labelTranscripts=TRUE)
plotTranscripts(ballgown::geneIDs(bg)[5790], bg, main=c('RPS8 in all adv_Infection samples'), sample=c('SRR12662528_expression', 'SRR12662551_expression', 'SRR12662556_expression'), labelTranscripts=TRUE)



# Close the PDF device where we have been saving our plots
dev.off()

# Exit the R session
quit(save="no")