# echo "full_length starts"
# ./paprica-run.sh full_length_testdata bacteria
# mkdir results/dustbin/full_length
# mkdir results/full_length
# mv full_length_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/full_length
# mv full_length_testdata.bacteria.edge_data.csv results/full_length
# mv full_length_testdata* results/dustbin/full_length

# echo "V1-V2 starts"
# ./paprica-run.sh V1-V2_testdata bacteria
# mkdir results/dustbin/V1-V2
# mkdir results/V1-V2
# mv V1-V2_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V1-V2
# mv V1-V2_testdata.bacteria.edge_data.csv results/V1-V2
# mv V1-V2_testdata* results/dustbin/V1-V2

echo "V1-V3 starts"
./paprica-run.sh V1-V3_testdata bacteria
# mkdir results/dustbin/V1-V3
# mkdir results/V1-V3
mv V1-V3_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V1-V3
mv V1-V3_testdata.bacteria.edge_data.csv results/V1-V3
mv V1-V3_testdata* results/dustbin/V1-V3

echo "V3-V4 starts"
./paprica-run.sh V3-V4_testdata bacteria
# mkdir results/dustbin/V3-V4
# mkdir results/V3-V4
mv V3-V4_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V3-V4
mv V3-V4_testdata.bacteria.edge_data.csv results/V3-V4
mv V3-V4_testdata* results/dustbin/V3-V4

echo "V4 starts"
./paprica-run.sh V4_testdata bacteria
# mkdir results/dustbin/V4
# mkdir results/V4
mv V4_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V4
mv V4_testdata.bacteria.edge_data.csv results/V4
mv V4_testdata* results/dustbin/V4

echo "V4-V5 starts"
./paprica-run.sh V4-V5_testdata bacteria
# mkdir results/dustbin/V4-V5
# mkdir results/V4-V5
mv V4-V5_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V4-V5
mv V4-V5_testdata.bacteria.edge_data.csv results/V4-V5
mv V4-V5_testdata* results/dustbin/V4-V5

echo "V6-V8 starts"
./paprica-run.sh V6-V8_testdata bacteria
# mkdir results/dustbin/V6-V8
# mkdir results/V6-V8
mv V6-V8_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V6-V8
mv V6-V8_testdata.bacteria.edge_data.csv results/V6-V8
mv V6-V8_testdata* results/dustbin/V6-V8

echo "V7-V9 starts"
./paprica-run.sh V7-V9_testdata bacteria
# mkdir results/dustbin/V7-V9
# mkdir results/V7-V9
mv V7-V9_testdata.bacteria.combined_16S.bacteria.tax.placements.csv results/V7-V9
mv V7-V9_testdata.bacteria.edge_data.csv results/V7-V9
mv V7-V9_testdata* results/dustbin/V7-V9