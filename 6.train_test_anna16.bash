#!/bin/bash

# This script trains and tests ANNA16 (Python)
python
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from math import sqrt
from anna16 import CopyNumberPredictor
import warnings
warnings.filterwarnings("ignore")

path = "data/cv/datasets/kmer_splits"
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9"]:
    rmse = []
    for i in range(0,5,1):
        X_train = pd.read_pickle(f"{path}/{region}_X_train_{i}.gz")
        Y_train = pd.read_pickle(f"{path}/{region}_Y_train_{i}.gz")
        X_test = pd.read_pickle(f"{path}/{region}_X_test_{i}.gz")
        Y_test = pd.read_pickle(f"{path}/{region}_X_test_{i}.gz")
        model = CopyNumberPredictor(region)
        model.fit(X_train,Y_train,verbose=False)
        pred = model.predict(X_test)
        rmse.append(sqrt(mean_squared_error(Y_test,pred)))
        print(f"{region}: {rmse[i]}")

        pred_records = pd.DataFrame(np.array([Y_test,pred]).transpose(), columns=["Y_test","Y_pred"])
        pred_records["group"] = i
        pred_records.to_csv("performance/cv/predictions/anna16/{region}_{i}.csv",index=False)
    performance[region] = rmse

pd.DataFrame(performance).to_csv("performance/cv/anna16.csv",index=False)

path = "data/test/datasets/kmer_splits"
performance = {}
for region in ["full_length","V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9"]:
    X_train = pd.read_pickle(f"{path}/{region}_X_train.gz")
    Y_train = pd.read_pickle(f"{path}/{region}_Y_train.gz")
    X_test = pd.read_pickle(f"{path}/{region}_X_test.gz")
    Y_test = pd.read_pickle(f"{path}/{region}_X_test.gz")
    model = CopyNumberPredictor(region)
    model.fit(X_train,Y_train,verbose=False)
    pred = model.predict(X_test)
    rmse = sqrt(mean_squared_error(Y_test,pred))
    print(f"{region}: {rmse}")

    pred_records = pd.DataFrame(np.array([Y_test,pred]).transpose(), columns=["Y_test","Y_pred"])
    pred_records["group"] = i
    pred_records.to_csv("performance/test/predictions/anna16/{region}.csv",index=False)
    performance[region] = rmse

pd.DataFrame(performance).to_csv("performance/test/anna16.csv")

quit()