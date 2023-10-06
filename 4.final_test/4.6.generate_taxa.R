fasta = seqinr::read.fasta(file = 'datasets/full_length.fasta', as.string = TRUE,
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
  reads = read.csv("datasets/full_length_testdata.csv")
  reads = dplyr::left_join(reads,taxa,by = "sequence")
  write.csv(reads,paste0("taxa/",tag,"_full_length_testtaxa.csv"),row.names=F)
}

regions = c("V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")
for (region in regions){
  fasta = seqinr::read.fasta(file = paste0('datasets/',region,".fasta"), as.string = TRUE,
                             forceDNAtolower = FALSE, whole.header = TRUE)
  fasta = data.frame(unlist(fasta))
  colnames(fasta) = "sequence"
  fasta$accession = row.names(fasta)
  fasta = fasta[(duplicated(fasta$sequence==F)),]
  input = fasta$sequence
  for (database in c("rdp_train_set_18.fa.gz","gg_13_8_train_set_97.fa.gz")){
    tag = strsplit(database,"_")[[1]][1]
    print(paste(region,tag,sep=" "))
    taxa <- dada2::assignTaxonomy(input, paste0("raw_data/",database))
    taxa = data.frame(taxa)
    taxa$accession = fasta$accession
    reads = read.csv("datasets/full_length_testdata.csv")
    reads = reads[,2:5]
    reads = dplyr::left_join(fasta,reads,by = "accession")
    reads = merge(reads,taxa,by = "accession")
    write.csv(reads,paste0("taxa/",tag,"_",region,"_testtaxa.csv"),row.names=F)
  }
}