place_seqs.py -s datasets/full_length_testdata.fasta -o tree/picrust_full_length.tre
hsp.py -i 16S -t tree/picrust_full_length.tre -o predictions/picrust_full_length.tsv.gz

echo "V1-V2 starts"
place_seqs.py -s datasets/V1-V2_testdata.fasta -o tree/picrust_V1-V2.tre
hsp.py -i 16S -t tree/picrust_V1-V2.tre -o predictions/picrust_V1-V2.tsv.gz

echo "V1-V3 starts"
place_seqs.py -s datasets/V1-V3_testdata.fasta -o tree/picrust_V1-V3.tre
hsp.py -i 16S -t tree/picrust_V1-V3.tre -o predictions/picrust_V1-V3.tsv.gz

echo "V3-V4 starts"
place_seqs.py -s datasets/V3-V4_testdata.fasta -o tree/picrust_V3-V4.tre
hsp.py -i 16S -t tree/picrust_V3-V4.tre -o predictions/picrust_V3-V4.tsv.gz

echo "V4 starts"
place_seqs.py -s datasets/V4_testdata.fasta -o tree/picrust_V4.tre
hsp.py -i 16S -t tree/picrust_V4.tre -o predictions/picrust_V4.tsv.gz

echo "V4-V5 starts"
place_seqs.py -s datasets/V4-V5_testdata.fasta -o tree/picrust_V4-V5.tre
hsp.py -i 16S -t tree/picrust_V4-V5.tre -o predictions/picrust_V4-V5.tsv.gz

echo "V6-V8 starts"
place_seqs.py -s datasets/V6-V8_testdata.fasta -o tree/picrust_V6-V8.tre
hsp.py -i 16S -t tree/picrust_V6-V8.tre -o predictions/picrust_V6-V8.tsv.gz

echo "V7-V9 starts"
place_seqs.py -s datasets/V7-V9_testdata.fasta -o tree/picrust_V7-V9.tre
hsp.py -i 16S -t tree/picrust_V7-V9.tre -o predictions/picrust_V7-V9.tsv.gz