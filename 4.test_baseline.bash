#!/bin/bash

# This script implements the baseline methods for the final test.

# rrnDB & CopyRighter (Python)
python
import pandas as pd
import numpy as np
from math import sqrt
from sklearn.metrics import mean_squared_error
from statistics import mean, stdev
from utils import CopyRighterSimulator, rrnDBSimulator

def test_rmse(model,X_test,Y_test):
    test_preds = model.predict(X_test)
    mse = mean_squared_error(Y_test, test_preds)
    rmse = sqrt(mse)
    return rmse

filtered = pd.read_csv("data/test/datasets/full_length_testdata.filtered.csv")

data = pd.read_csv("data/test/datasets/copyrighterdata.txt",sep="\t")
raw_taxa = data["taxonomy"].values.tolist()
copy_number = data["16S rRNA count"].values.tolist()
CopyRighterData = []
for i in range(0,len(raw_taxa)):
    line = raw_taxa[i]
    lineage = []
    lineage.append(line.replace("; ",";"))
    lineage.append(copy_number[i])
    CopyRighterData.append(lineage)

crs = CopyRighterSimulator(CopyRighterData)

crs_performance = {}
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9"]:
    greengenetaxa = pd.read_csv("data/test/taxa/gg_"+region+"_testtaxa.csv")
    greengenetaxa = greengenetaxa.merge(filtered[["accession"]],how="inner",on = ["accession"])
    greengenetaxa = greengenetaxa.dropna(subset=["copy_number"])
    ggX = greengenetaxa[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', "Species"]]
    ggY = greengenetaxa['copy_number']
    newggX = []
    for line in ggX.values.tolist():
        lineage=""
        for element in line:
            lineage+=str(element)
            lineage+=";"
        lineage=lineage[:-1]
        newggX.append(lineage)
    crs_performance[region] = test_rmse(crs,newggX,ggY)

rrnDBData = pd.read_csv("data/test/datasets/rrnDB_pan-taxa stats_RDP.csv")
rrnDBData = rrnDBData.values.tolist()
rds = rrnDBSimulator(rrnDBData)
rds_performance = {}
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9"]:
    rdptaxa = pd.read_csv("data/test/taxa/rdp_"+region+"_testtaxa.csv")
    rdptaxa = rdptaxa.merge(filtered[["accession"]],how="inner",on = ["accession"])
    rdptaxa = rdptaxa.dropna(subset=["copy_number"])
    ggX = rdptaxa[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']]
    ggY = rdptaxa['copy_number']
    rds_performance[region] = test_rmse(rds,ggX.values.tolist(),ggY.tolist())

pd.concat([pd.DataFrame(crs_performance,index=["CopyRighter"]).transpose(),
pd.DataFrame(rds_performance,index=["rrnDB"]).transpose()],axis = 1).to_csv("performance/test/cprt_rrndb.csv")

quit()

# PICRUSt2 (Bash)
cd data/test
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

cd ../..

# PAPRICA (Bash)
path="data/test/datasets"
./paprica-run.sh $path/full_length_testdata bacteria
./paprica-run.sh $path/V1-V2_testdata bacteria
./paprica-run.sh $path/V1-V3_testdata bacteria
./paprica-run.sh $path/V3-V4_testdata bacteria
./paprica-run.sh $path/V4_testdata bacteria
./paprica-run.sh $path/V4-V5_testdata bacteria
./paprica-run.sh $path/V6-V8_testdata bacteria
./paprica-run.sh $path/V7-V9_testdata bacteria

mv *_testdata.bacteria.edge_data.csv performance/test/predictions/paprica
rm *_testdata*
