cd genomes
hyperex -p tmp_V1-V2 -m 3 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer GCTGCCTCCCGTAGGAGT mock_full_length.fasta 
hyperex -p tmp_V1-V3 -m 3 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer ATTACCGCGGCTGCTGG mock_full_length.fasta 
hyperex -p tmp_V3-V4 -m 3 --forward-primer CCTACGGGAGGCAGCAG --reverse-primer GACTACHVGGGTATCTAATCC mock_full_length.fasta 
hyperex -p tmp_V4 -m 3 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer GGACTACHVGGGTWTCTAAT mock_full_length.fasta 
hyperex -p tmp_V4-V5 -m 3 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer CCGYCAATTYMTTTRAGTTT mock_full_length.fasta 
hyperex -p tmp_V6-V8 -m 3 --forward-primer GAATTGACGGGGGCCCGCACAAG --reverse-primer CGGTGTGTACAAGGCCCGGGAACG mock_full_length.fasta 
hyperex -p tmp_V7-V9 -m 3 --forward-primer CAACGAGCGCAACCCT --reverse-primer TACGGYTACCTTGTTACGACT mock_full_length.fasta 
rm *.gff
awk -F" " '{print $1}' tmp_V1-V2.fa > V1-V2.fasta
awk -F" " '{print $1}' tmp_V1-V3.fa > V1-V3.fasta
awk -F" " '{print $1}' tmp_V3-V4.fa > V3-V4.fasta
awk -F" " '{print $1}' tmp_V4.fa > V4.fasta
awk -F" " '{print $1}' tmp_V4-V5.fa > V4-V5.fasta
awk -F" " '{print $1}' tmp_V6-V8.fa > V6-V8.fasta
awk -F" " '{print $1}' tmp_V7-V9.fa > V7-V9.fasta
rm *.fa