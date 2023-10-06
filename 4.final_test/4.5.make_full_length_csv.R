reads = read.csv("datasets/oriented_testdata.csv")
fasta = seqinr::read.fasta(file = 'datasets/full_length_testdata.fasta', as.string = TRUE,
                                   forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
full_length = dplyr::left_join(fasta,reads[,c(1,3:5)],by = "accession")
write.csv(full_length,"datasets/full_length_testdata.csv",row.names = F)

paprica_bact = read.csv("datasets/paprica.bacteria.csv")
check_bact = dplyr::anti_join(full_length,paprica_bact,by = "accession")
write.csv(check_bact,"datasets/full_length_testdata.filtered.csv",row.names = F)
