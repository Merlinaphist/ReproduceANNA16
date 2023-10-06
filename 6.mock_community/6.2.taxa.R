library(dplyr)
even = read.csv("processed/even.csv")
input = even$sequence
for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
  tag = strsplit(database,"_")[[1]][1]
  taxa <- dada2::assignTaxonomy(input, paste0("/Users/miaojiazheng/Desktop/Projects/Biogeochemistry/ACE/Machine\ Learning/remake/raw_data/",database))
  taxa = data.frame(taxa)
  taxa$GenBank.ID = even$GenBank.ID
  write.csv(taxa,paste0("taxa/",tag,"_full_length_taxa.csv"),row.names=F)
}

regions = c("V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")
for (region in regions){
  fasta = seqinr::read.fasta(file = paste0('genomes/',region,".fasta"), as.string = TRUE,
                             forceDNAtolower = FALSE, whole.header = TRUE)
  fasta = data.frame(unlist(fasta))
  colnames(fasta) = "sequence"
  fasta$accession = row.names(fasta)
  # fasta = fasta[(duplicated(fasta$sequence==F)),]
  input = fasta$sequence
  for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
    tag = strsplit(database,"_")[[1]][1]
    print(paste(region,tag,sep=" "))
    taxa <- dada2::assignTaxonomy(input, paste0("/Users/miaojiazheng/Desktop/Projects/Biogeochemistry/ACE/Machine\ Learning/remake/raw_data/",database))
    taxa = data.frame(taxa)
    taxa$GenBank.ID = fasta$accession
    write.csv(taxa,paste0("taxa/",tag,"_",region,"_taxa.csv"),row.names=F)
  }
}