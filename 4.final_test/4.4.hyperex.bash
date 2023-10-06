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