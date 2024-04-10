#!/bin/bash

# This script calculates the Shapley value
python
import pandas as pd
import numpy as np
from utils import Preprocessing
from anna16 import CopyNumberPredictor
import shap

da0 = pd.read_csv("data/cv/datasets/full_length_reads.csv")
da1 = pd.read_csv("data/test/datasets/full_length_testdata.csv")
da = pd.concat([da0,da1],axis=0)
model = CopyNumberPredictor(region = "full_length")
da.index = [i for i in range(da.shape[0])]
X = da["sequence"]
Y = da['copy_number']
pp = Preprocessing()
X = pp.CountKmers(seqs=X)
X_with_indices = pd.concat([pd.DataFrame(X, columns = pp.vectorizer.get_feature_names_out()),
           pd.DataFrame([i for i in range(27579)],columns = ["index"])],axis = 1)
X_train,X_test,Y_train,Y_test = train_test_split(X_with_indices, Y, test_size=0.1, random_state=0)
X_test.iloc[:,4096:4097].to_csv("explainability/selected_X_test.csv",index=False)
X_test = X_test.iloc[:,0:4096]
X_train = X_train.iloc[:,0:4096]
X_train_summary = shap.kmeans(X_train, 10)
# X_summary = shap.kmeans(X, 10)
model.fit(X,Y)
# X = pd.DataFrame(X,columns = pp.vectorizer.get_feature_names_out())
explainer = shap.KernelExplainer(model.predict, X_train_summary)
shap_values = explainer.shap_values(X_test)
pd.DataFrame(shap_values,columns = pp.vectorizer.get_feature_names_out()).to_csv("explainability/shap_values.csv",index = False)
quit()