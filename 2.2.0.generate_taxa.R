for (region in c("full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")){
  fasta = seqinr::read.fasta(file = paste0('datasets/',region,'.fasta'), as.string = TRUE,
                             forceDNAtolower = FALSE, whole.header = TRUE)
  fasta = data.frame(unlist(fasta))
  colnames(fasta) = "sequence"
  fasta$accession = row.names(fasta)
  fasta = fasta[(duplicated(fasta$sequence==F)),]
  input = fasta$sequence
  for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
    tag = strsplit(database,"_")[[1]][1]
    taxa <- dada2::assignTaxonomy(input, paste0("raw_data/",database))
    taxa = data.frame(taxa)
    taxa$sequence = row.names(taxa)
    reads = read.csv("datasets/full_length_reads.csv")
    reads = dplyr::left_join(reads,taxa,by = "sequence")
    write.csv(reads,paste0("taxa/",tag,"_",region,"_taxa.csv"),row.names=F)
  }
}