cd genomes
tar -xf strain_genomes.tar.gz
files=`echo *.fa`
for file in $files
do
hyperex -p ${file}_1_tmp -m 3 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer TACGGYTACCTTGTTACGACT $file
hyperex -p ${file}_2_tmp -m 3 --forward-primer TACGGYTACCTTGTTACGACT --reverse-primer AGAGTTTGATCCTGGCTCAG $file
awk -F" " '{print $1}' ${file}_1_tmp.fa > ${file}_1_extracted.fasta
awk -F" " '{print $1}' ${file}_2_tmp.fa > ${file}_2_extracted.fasta
done
rm *.gff
rm *_tmp.fa
rm *.fa
cat *1_extracted.fasta > strain_fwd_16S.fasta
cat *2_extracted.fasta > strain_rev_16S.fasta
rm *extracted.fasta