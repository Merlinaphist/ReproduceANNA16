#!/bin/bash

# This script records the codes for preprocessing of cross-validation data.

# Data Cleaning (R)
R
reads = seqinr::read.fasta(file = 'data/raw_data/rrnDB-5.7_16S_rRNA.fasta', as.string = TRUE,
                           forceDNAtolower = FALSE, whole.header = TRUE)
reads = data.frame(unlist(reads))
reads= cbind(row.names(reads),reads)
row.names(reads) = seq(1,nrow(reads))
colnames(reads) = c('name',"sequence")
reads$accession = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,2)
reads$name = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,1)

reads$name = stringr::str_replace_all(reads$name,"'","")
reads$name = stringr::str_replace_all(reads$name,"\\]","")
reads$name = stringr::str_replace_all(reads$name,"\\[","")

meta = read_tsv("raw_data/rrnDB-5.7.tsv")
meta = dplyr::select(meta,1:2,12)
colnames(meta) = c("accession","taxid","copy_number")
reads = dplyr::left_join(reads,meta,by = "accession")
rm(meta)
reads = reads[(duplicated(reads$accession)==F),]
set.seed(114514)
reads = reads[sample(1:nrow(reads)),]

write.csv(reads,"data/cv/datasets/cleaned_reads.csv",row.names = F)

fa = character(2*nrow(reads))
fa[c(TRUE, FALSE)] = sprintf(">%s", reads$accession)
fa[c(FALSE, TRUE)] = reads$sequence
writeLines(fa,"data/cv/blast/cleaned_reads.fasta")
q()
n

# BLASTN (Bash)
cd data/cv/blast
makeblastdb -dbtype nucl -in marker.fasta -out test
blastn -db test -query cleaned_reads.fasta -task blastn-short -outfmt "6 qseqid qlen qstart qend stitle slen sstart send pident bitscore evalue" -out output.tsv -max_hsps 1
echo -e "qseqid\tqlen\tqstart\tqend\tstitle\tslen\tsstart\tsend\tpident\tbitscore\tevalue" | cat - output.tsv > results.tsv
cd ../../..

# Unify Orientation (R)
R
library(readr)
library(dplyr)
library(stringr)
library(Biostrings)
reads = read.csv("data/cv/datasets/cleaned_reads.csv")
blast_results = read_tsv("data/cv/blast/results.tsv")
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

write.csv(reads,"data/cv/datasets/oriented_reads.csv",row.names = F)

fa = character(2*nrow(reads))
fa[c(TRUE, FALSE)] = sprintf(">%s", reads$accession)
fa[c(FALSE, TRUE)] = reads$sequence
writeLines(fa,"data/cv/datasets/oriented_reads.fasta")

q()
n

# Hyperex (Bash)
cd data/cv
hyperex -p datasets/tmp_full_length -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V1-V2 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer GCTGCCTCCCGTAGGAGT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V1-V3 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer ATTACCGCGGCTGCTGG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V3-V4 -m 2 --forward-primer CCTACGGGAGGCAGCAG --reverse-primer GACTACHVGGGTATCTAATCC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V4 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer GGACTACHVGGGTWTCTAAT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V4-V5 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer CCGYCAATTYMTTTRAGTTT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V6-V8 -m 2 --forward-primer GAATTGACGGGGGCCCGCACAAG --reverse-primer CGGTGTGTACAAGGCCCGGGAACG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V7-V9 -m 2 --forward-primer CAACGAGCGCAACCCT --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_reads.fasta 
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
reads = read.csv("data/cv/datasets/oriented_reads.csv")
fasta = seqinr::read.fasta(file = 'data/cv/datasets/full_length.fasta', as.string = TRUE,
                                   forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
full_length = dplyr::left_join(fasta,reads[,c(1,3:5)],by = "accession")
write.csv(full_length,"data/cv/datasets/full_length_reads.csv",row.names = F)
q()
n

# Construct Phylogenetic Trees (Bash)
cd data/cv
## full_length
mafft datasets/full_length.fasta > tree/full_length_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/full_length_aln.fasta > tree/full_length.tre
## V1-V2
mafft datasets/V1-V2.fasta > tree/V1-V2_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V1-V2_aln.fasta > tree/V1-V2.tre
## V1-V3
mafft datasets/V1-V3.fasta > tree/V1-V3_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V1-V3_aln.fasta > tree/V1-V3.tre
## V3-V4
mafft datasets/V3-V4.fasta > tree/V3-V4_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V3-V4_aln.fasta > tree/V3-V4.tre
## V4
mafft datasets/V4.fasta > tree/V4_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V4_aln.fasta > tree/V4.tre
## V4-V5
mafft datasets/V4-V5.fasta > tree/V4-V5_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V4-V5_aln.fasta > tree/V4-V5.tre
## V6-V8
mafft datasets/V6-V8.fasta > tree/V6-V8_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V6-V8_aln.fasta > tree/V6-V8.tre
## V7-V9
mafft datasets/V7-V9.fasta > tree/V7-V9_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V7-V9_aln.fasta > tree/V7-V9.tre
cd ../..

# Split Data for Machine Learning (Python)
python
import numpy as np
import pandas as pd
from utils import Preprocessing

def read_region(region):
    da = pd.read_csv("data/cv/datasets/full_length_reads.csv")
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

def split_test_train(da, multiplicand, region, path="data/cv/datasets/kmer_splits"):
    pp = Preprocessing()
    for i in range(0,5,1):
        X_test = da["sequence"].loc[i*multiplicand:(i+1)*multiplicand]
        Y_test = da["copy_number"].loc[i*multiplicand:(i+1)*multiplicand]
        X_train = pd.concat([da["sequence"].loc[0:i*multiplicand],
                             da["sequence"].loc[(i+1)*multiplicand:]],axis = 0)
        Y_train = pd.concat([da["copy_number"].loc[0:i*multiplicand],
                             da["copy_number"].loc[(i+1)*multiplicand:]],axis = 0)
        X_train = pp.CountKmers(seqs=X_train)
        X_test = pp.CountKmers(seqs=X_test)
        pd.DataFrame(X_train).to_pickle(f"{path}/{region}_X_train_{i}.gz")
        pd.DataFrame(Y_train).to_pickle(f"{path}/{region}_Y_train_{i}.gz")
        pd.DataFrame(X_test).to_pickle(f"{path}/{region}_X_test_{i}.gz")
        pd.DataFrame(Y_test).to_pickle(f"{path}/{region}_Y_test_{i}.gz")
    return

pp = Preprocessing()
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9"]:
    rmse = []
    da = read_region(region)
    multiplicand = int(da.shape[0]*0.2)+1
    split_test_train(da, multiplicand, region)

quit()

# Taxonomic Classification (R)
R
for (region in c("full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")){
  fasta = seqinr::read.fasta(file = paste0('data/cv/datasets/',region,'.fasta'), as.string = TRUE,
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
    reads = read.csv("data/cv/datasets/full_length_reads.csv")
    reads = dplyr::left_join(reads,taxa,by = "sequence")
    write.csv(reads,paste0("data/cv/taxa/",tag,"_",region,"_taxa.csv"),row.names=F)
  }
}
q()
n