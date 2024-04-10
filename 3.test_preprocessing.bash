#!/bin/bash

# This script records the codes for preprocessing of test data.

# Data Cleaning (R)
R
library(seqinr)
library(readr)
library(dplyr)
library(Biostrings)
library(DECIPHER)
library(ggplot2)
library(stringr)
reads = seqinr::read.fasta(file = 'data/raw_data/rrnDB-5.8_16S_rRNA.fasta', as.string = TRUE,
                           forceDNAtolower = FALSE, whole.header = TRUE)
reads = data.frame(unlist(reads))
reads= cbind(row.names(reads),reads)
row.names(reads) = seq(1,nrow(reads))
colnames(reads) = c('name',"sequence")
reads$accession = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,2)
reads$name = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,1)

reads$name = str_replace_all(reads$name,"'","")
reads$name = str_replace_all(reads$name,"\\]","")
reads$name = str_replace_all(reads$name,"\\[","")


meta = read_tsv("raw_data/rrnDB-5.8.tsv")
meta = select(meta,1:2,12)
colnames(meta) = c("accession","taxid","copy_number")
reads = left_join(reads,meta,by = "accession")
rm(meta)
reads = reads[(duplicated(reads$accession)==F),]
set.seed(114514)
reads = reads[sample(1:nrow(reads)),]

traindata = read.csv("data/cv/datasets/full_length_reads.csv")

testdata = reads[!(reads$accession%in%traindata$accession),]

write.csv(testdata,"data/test/datasets/cleaned_testdata.csv",row.names = F)

fa = character(2*nrow(testdata))
fa[c(TRUE, FALSE)] = sprintf(">%s", testdata$accession)
fa[c(FALSE, TRUE)] = testdata$sequence
writeLines(fa,"data/test/blast/cleaned_testdata.fasta")

q()
n

# BLASTN (Bash)
cd data/test/blast
makeblastdb -dbtype nucl -in marker.fasta -out test
blastn -db test -query cleaned_testdata.fasta -task blastn-short -outfmt "6 qseqid qlen qstart qend stitle slen sstart send pident bitscore evalue" -out output.tsv -max_hsps 1
echo -e "qseqid\tqlen\tqstart\tqend\tstitle\tslen\tsstart\tsend\tpident\tbitscore\tevalue" | cat - output.tsv > results.tsv
cd ../../..

# Unify Orientation (R)
R
library(readr)
library(dplyr)
library(stringr)
library(Biostrings)
reads = read.csv("data/test/datasets/cleaned_testdata.csv")
blast_results = read_tsv("data/test/blast/results.tsv")
colnames(blast_results)[1] = "accession"
reads = left_join(reads,blast_results, by = "accession")

reads = reads%>%
  filter(evalue<0.01)

rev_reads = reads[(reads$sstart>reads$send),1:5]
fwd_reads = reads[(reads$sstart<reads$send),1:5]

rc = c()
for (seq in rev_reads$sequence){
  rc = append(rc,as.character(reverseComplement(DNAString(seq))))
}

rev_reads$sequence = rc

reads = rbind(fwd_reads,rev_reads)

set.seed(114514)
reads = reads[sample(1:nrow(reads)),]

write.csv(reads,"data/test/datasets/oriented_testdata.csv",row.names = F)

fa = character(2*nrow(reads))
fa[c(TRUE, FALSE)] = sprintf(">%s", reads$accession)
fa[c(FALSE, TRUE)] = reads$sequence
writeLines(fa,"data/test/datasets/oriented_testdata.fasta")

q()
n

# Hyperex (Bash)
cd data/test
hyperex -p datasets/tmp_full_length -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V1-V2 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer GCTGCCTCCCGTAGGAGT datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V1-V3 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer ATTACCGCGGCTGCTGG datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V3-V4 -m 2 --forward-primer CCTACGGGAGGCAGCAG --reverse-primer GACTACHVGGGTATCTAATCC datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V4 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer GGACTACHVGGGTWTCTAAT datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V4-V5 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer CCGYCAATTYMTTTRAGTTT datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V6-V8 -m 2 --forward-primer GAATTGACGGGGGCCCGCACAAG --reverse-primer CGGTGTGTACAAGGCCCGGGAACG datasets/oriented_testdata.fasta 
hyperex -p datasets/tmp_V7-V9 -m 2 --forward-primer CAACGAGCGCAACCCT --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_testdata.fasta 
rm datasets/*.gff
awk -F" " '{print $1}' datasets/tmp_full_length.fa > datasets/full_length.fasta
awk -F" " '{print $1}' datasets/tmp_V1-V2.fa > datasets/V1-V2.fasta
awk -F" " '{print $1}' datasets/tmp_V1-V3.fa > datasets/V1-V3.fasta
awk -F" " '{print $1}' datasets/tmp_V3-V4.fa > datasets/V3-V4.fasta
awk -F" " '{print $1}' datasets/tmp_V4.fa > datasets/V4.fasta
awk -F" " '{print $1}' datasets/tmp_V4-V5.fa > datasets/V4-V5.fasta
awk -F" " '{print $1}' datasets/tmp_V6-V8.fa > datasets/V6-V8.fasta
awk -F" " '{print $1}' datasets/tmp_V7-V9.fa > datasets/V7-V9.fasta
rm datasets/*.fa
cd ../..

# Convert Trimmed Full Length from FASTA to CSV (R)
R
reads = read.csv("data/test/datasets/oriented_testdata.csv")
fasta = seqinr::read.fasta(file = 'datasets/full_length_testdata.fasta', as.string = TRUE,
                                   forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
full_length = dplyr::left_join(fasta,reads[,c(1,3:5)],by = "accession")
write.csv(full_length,"data/test/datasets/full_length_testdata.csv",row.names = F)

paprica_bact = read.csv("data/raw_data/paprica.bacteria.csv") # Remove paprica overlap
check_bact = dplyr::anti_join(full_length,paprica_bact,by = "accession")
write.csv(check_bact,"data/test/datasets/full_length_testdata.filtered.csv",row.names = F)
q()
n

# Taxonomic Classification (R)
R
fasta = seqinr::read.fasta(file = 'data/test/datasets/full_length.fasta', as.string = TRUE,
                           forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
fasta = fasta[(duplicated(fasta$sequence==F)),]
input = fasta$sequence

for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
  tag = strsplit(database,"_")[[1]][1]
  taxa <- dada2::assignTaxonomy(input, paste0("data/raw_data/",database))
  taxa = data.frame(taxa)
  taxa$sequence = row.names(taxa)
  reads = read.csv("data/test/datasets/full_length_testdata.csv")
  reads = dplyr::left_join(reads,taxa,by = "sequence")
  write.csv(reads,paste0("data/test/taxa/",tag,"_full_length_testtaxa.csv"),row.names=F)
}

regions = c("V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")
for (region in regions){
  fasta = seqinr::read.fasta(file = paste0('data/test/datasets/',region,".fasta"), as.string = TRUE,
                             forceDNAtolower = FALSE, whole.header = TRUE)
  fasta = data.frame(unlist(fasta))
  colnames(fasta) = "sequence"
  fasta$accession = row.names(fasta)
  fasta = fasta[(duplicated(fasta$sequence==F)),]
  input = fasta$sequence
  for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
    tag = strsplit(database,"_")[[1]][1]
    print(paste(region,tag,sep=" "))
    taxa <- dada2::assignTaxonomy(input, paste0("data/raw_data/",database))
    taxa = data.frame(taxa)
    taxa$accession = fasta$accession
    reads = read.csv("data/test/datasets/full_length_testdata.csv")
    reads = reads[,2:5]
    reads = dplyr::left_join(fasta,reads,by = "accession")
    reads = merge(reads,taxa,by = "accession")
    write.csv(reads,paste0("data/test/taxa/",tag,"_",region,"_testtaxa.csv"),row.names=F)
  }
}
q()
n


# Kmer for Machine Learning (Python)
python
import numpy as np
import pandas as pd
from utils import Preprocessing

def read_region_train(region):
    da = pd.read_csv("data/cv/full_length_reads.csv")
    if region == "full_length":
        return da
    file_handle = open("data/cv/datasets/"+region+".fasta","r")
    seq = []
    seqid = []
    tmp_seq = ""
    for line in file_handle:
        if (line[0] == ">"):
            if tmp_seq != "":
                seq.append(tmp_seq)
            seqid.append(line.split("\n")[0][1:])
            tmp_seq = ""
        else:
            tmp_seq+=line.split("\n")[0]
    seq.append(tmp_seq)
    file_handle.close()
    sub = pd.DataFrame([seq,seqid], index = [region,"accession"])
    sub = sub.transpose()
    da = da[["accession","copy_number"]]
    da = pd.merge(da,sub,on="accession",how="inner")
    return da

def read_region_test(region):
    da = pd.read_csv("data/test/datasets/full_length_testdata.filtered.csv")
    if region == "full_length":
        return da
    file_handle = open("data/test/datasets/"+region+"_testdata.fasta","r")
    seq = []
    seqid = []
    tmp_seq = ""
    for line in file_handle:
        if (line[0] == ">"):
            if tmp_seq != "":
                seq.append(tmp_seq)
            seqid.append(line.split("\n")[0][1:])
            tmp_seq = ""
        else:
            tmp_seq+=line.split("\n")[0]
    seq.append(tmp_seq)
    file_handle.close()
    sub = pd.DataFrame([seq,seqid], index = [region,"accession"])
    sub = sub.transpose()
    da = da[["accession","copy_number"]]
    da = pd.merge(da,sub,on="accession",how="inner")
    return da

pp = Preprocessing()
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4-V5","V4","V6-V8","V7-V9"]:
    da0 = read_region_train(region)
    X0 = da0[region]
    Y0 = da0['copy_number']
    X0 = pp.CountKmers(X0)
    
    da1 = read_region_test(region)
    X1 = da1[region]
    Y1 = da1['copy_number']
    X1 = pp.CountKmers(X1)

    pd.DataFrame(X0).to_pickle(f"data/test/datasets/{region}_X_final_train.gz")
    pd.DataFrame(Y0).to_pickle(f"data/test/datasets/{region}_Y_final_train.gz")
    pd.DataFrame(X1).to_pickle(f"data/test/datasets/{region}_X_final_test.gz")
    pd.DataFrame(Y1).to_pickle(f"data/test/datasets/{region}_Y_final_test.gz")


quit()
