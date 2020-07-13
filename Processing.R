### Load required packages -----------------------------------------------------
library(Biostrings)
library(xlsx)
library(openxlsx)
library(ape)
library(tidyverse)

#### Load data -----------------------------------------------------------------
milnieLTP <- Biostrings::readRNAStringSet("Source_data/Milner2016/arb-silva.de_2019-09-17_id709787_tax_silva.fasta")
accnos.tre<- readxl::read_xlsx("Source_data/Overview_seq_sources.xlsx",sheet="Accessions",
                               col_names = c("Entry","Strain","Accession"))
### Process data ---------------------------------------------------------------
milnieLTP.DNA <- Biostrings::DNAStringSet(milnieLTP) # RNA to DNA

# extract parts
seqids <- names(milnieLTP.DNA)
reads  <- paste(milnieLTP.DNA)  
length <- Biostrings::width(milnieLTP.DNA)


accnos <- sub("(.+?)\\..*$","\\1",seqids)
specn  <- sub(".*;(.*)$","\\1",seqids)

resdf <- data.frame(ID=accnos,Species=specn,Strain=NA,
                    Fasta=reads,Length=length,Source="SILVA LTP")

### Write data  ----------------------------------------------------------------
openxlsx::write.xlsx(resdf,"Source_data/Milner2016/LTPmilner.xlsx")


### Get sequences from Genbank  ------------------------------------------------
reads <- ape::read.GenBank(access.nb=accnos.tre$Accession,species.names=TRUE,
                           as.character = TRUE)

ids <- names(reads)
seqs <- sapply(reads,function(x)toupper(paste0(x,collapse = "")))

resdfgb <- data.frame(queryid=ids,sequence=seqs,accnos.tre$Entry,accnos.tre$Strain,accnos.tre$Accession)

openxlsx::write.xlsx(resdfgb,"Source_data/Milner2016/GenbankMilner.xslx")


### Convert Eddie fasta to excel------------------------------------------------
eddiefasta <- Biostrings::readDNAStringSet("Source_data/Eddie2016/Eddie2018treeDNA.fasta")
eddiedf <- data.frame(ID=sub("(.?)\ .*","\\1",names(eddiefasta)),
                      Sequence=paste(eddiefasta),
                      SeqLen=width(eddiefasta))
# openxlsx::write.xlsx(eddiedf,file = "Source_data/Overview_seq_sources.xlsx",sheetName="Eddie2018tree")

### Read in seq overview, filter >1200 and write to fasta ----------------------
referenceSeqs <- readxl::read_xlsx("Source_data/Overview_seq_sources.xlsx",
                                   sheet="Selected")

hist(referenceSeqs$SeqLen,breaks=25,col="lightblue",
     main="Histogram of sequence lengths before filtering",
     xlab="Sequence length (number of nucleotides)")
box()

referenceslong <- referenceSeqs %>% dplyr::filter(SeqLen>1200L & 
                                                    Classis=="Gammaproteobacteria")
outgroup <- referenceSeqs %>% dplyr::filter(SeqLen>1200L & 
                                              Genus=="Rhizobium")

refseqlong_wOG <- rbind(referenceslong,outgroup)

hist(referenceslong$SeqLen,breaks=25,col="lightblue",
     main="Histogram of sequence lengths after length filtering",
     xlab="Sequence length (number of nucleotides >1200)")
box()

bpseqlen <- data.frame(beforeFilter=referenceSeqs$SeqLen,
                       afterFilter=c(referenceslong$SeqLen,
                                     rep(NA,(length(referenceSeqs$SeqLen)-
                                               length(referenceslong$SeqLen)))))
bpseqlen.m <- reshape2::melt(bpseqlen)
with(bpseqlen.m,boxplot(value~variable,col="lightblue",
                        xlab="",ylab="sequence length"))


dna3 <- DNAStringSet(refseqlong_wOG$Sequence)
names(dna3) <-  gsub(" ","_",refseqlong_wOG$Accession)
Biostrings::writeXStringSet(dna3,"Source_data/Reference_seqs_all_wOG.fasta")


# Programmatically create label file for tree ----------------------------------
flnwk <- read.tree("Source_data/fullLengthWithOG/RAxML_bipartitions.FulllennoOTUwOGMLtreewithBPinfo")
write.csv2(data.frame(tiplabels=flnwk$tip.label),
           "Source_data/fullLengthWithOG/tiplabels.csv")
treefl_labels <- readxl::read_xlsx("Source_data/Overview_seq_sources.xlsx",
                                   sheet = "Tiplabels")
write.table(treefl_labels,file="Source_data/fullLengthWithOG/tiplabels_full_length.tsv",row.names = FALSE,
            col.names = FALSE,quote = FALSE,sep="\t")

# Reduced fasta for SILVA ACT --------------------------------------------------
a <- ape::read.dna("Source_data/fullLengthWithOG/full_lengthwithOG_minimizedDNA.fasta.reduced",
                   format="sequential")
ape::write.dna(a,format="fasta",file = "Source_data/fullLengthWithOG/full_lengthwithOG_minimizedDNA.fasta.reduced.fasta")


# Make sure difference in capitalization does not affect SINA ------------------
flwogOTU <- Biostrings::readDNAStringSet("Source_data/fullLengthWithOG/full_lengthwithOG_withOTU_unaligned.fasta")
writeXStringSet(flwogOTU,"Source_data/fullLengthWithOG/FLwOGwOTU_unaligned_clean.fasta")

# Selecting only minimized OTUs for EPA-NG -------------------------------------
minimwOTU <- Biostrings::readDNAStringSet("Source_data/fullLengthWithOG/FLwOTUwOG_minimizedDNA.fasta")
refwog    <- Biostrings::readDNAStringSet("Source_data/fullLengthWithOG/full_lengthwithOG_minimizedDNA.fasta.reduced.fasta")

names(refwog) %in% names(minimwOTU) # sanity check

OTUonly <- names(minimwOTU)[!(names(minimwOTU) %in% names(refwog))]
otuseqminaln <- minimwOTU[OTUonly]
Biostrings::writeXStringSet(otuseqminaln,"Source_data/fullLengthWithOG/OTUonlyMinimized.fasta")
noOTUs <- names(minimwOTU)[(names(minimwOTU) %in% names(refwog))]
nootuseqminaln <- minimwOTU[noOTUs]
Biostrings::writeXStringSet(nootuseqminaln,"Source_data/fullLengthWithOG/noOTUMinimized.fasta")
