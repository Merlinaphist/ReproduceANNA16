#!/bin/bash

# This script searches the hyperparameters for ANNA16 (Python)
python
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from math import sqrt
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import matplotlib.pyplot as plt
from tensorflow.keras.optimizers import Adam
from tensorflow import cast,float32
from keras import backend as kb
from statistics import mean, stdev
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge, Lasso
from smac.facade.smac_bb_facade import SMAC4BB
from smac.scenario.scenario import Scenario
from smac.configspace import ConfigurationSpace
from ConfigSpace.hyperparameters import (
    CategoricalHyperparameter,
    UniformFloatHyperparameter,
    UniformIntegerHyperparameter,
)

import warnings
warnings.filterwarnings("ignore")

class CopyNumberPredictor():
    def __init__(self,config):
        self.state =[config["epochs"], config["learning_rate"], config["validation_split"],
                     config["batch_size"], config["n_neurons_0"], config["activation_0"], 
                     config["n_neurons_1"], config["activation_1"], config["n_neurons_2"], 
                     config["activation_2"], config["n_neurons_3"], config["activation_3"],
                     config["n_neurons_4"], config["activation_4"], config["n_neurons_5"], 
                     config["activation_5"]]
        # self.state = [55, 0.0001, 0.1, 100, 469, 1, 956, 0, 649, 0, 941, 3, 597, 1, 65, 4]
        self.activation_indices={0:"relu",1:"gelu",2:"selu",3:"elu",4:"linear"}
        self.mlp = self.create_mlp()
        self.ridge = Ridge(alpha = config["alpha"])
        self.pca = PCA(n_components=100)
        self.svr = SVR(kernel='rbf',C=config["C"],gamma=config["gamma"])
        # self.svr = thundersvm.SVR(kernel='rbf',C=config["C"],gamma=config["gamma"])
        
    def save(self,filename):
        if filename[-4:] != ".zip":
            raise ValueError('Invalid filename. Expect a zip file.')
            
        path = filename[:-4]
        prefix = path + "/" + self.region
            
        if not os.path.exists(path):
            os.makedirs(path)
        
        with open(prefix+'_pca.pkl','wb') as file:
            pickle.dump(self.pca,file)
        with open(prefix+'_ridge.pkl','wb') as file:
            pickle.dump(self.ridge,file)
        with open(prefix+'_svr.pkl','wb') as file:
            pickle.dump(self.svr,file)
        self.mlp.save(prefix+"_mlp.h5")
        
        shutil.make_archive(path, 'zip', path)
        shutil.rmtree(path)
            
    def load(self,filename):
        if filename[-4:] != ".zip":
            raise ValueError('Invalid input file. Expect a zip file.')
            
        path = filename[:-4]
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        with ZipFile(filename,'r') as zObject:
            zObject.extractall(path=path)
        
        prefix = path + "/" + self.region
        with open(prefix+'_pca.pkl', 'rb') as file:
            self.pca = pickle.load(file) 
        with open(prefix+'_ridge.pkl', 'rb') as file:
            self.ridge = pickle.load(file)
        with open(prefix+'_svr.pkl', 'rb') as file:
            self.svr = pickle.load(file)
        self.mlp = load_model(prefix + "_mlp.h5",
                              custom_objects={"root_mean_squared_error": self.root_mean_squared_error})
        shutil.rmtree(path)
        
    def fit(self,X_train,Y_train,verbose=True):
        self.echo(text = "------Training Starts------", verbose = verbose)
        X_train_pca = self.pca.fit_transform(X_train)
        self.fit_mlp(X_train,Y_train)
        mlp_pred = self.mlp.predict(X_train,verbose=0)
        self.echo(text = "Model 1: MLP done.", verbose = verbose)
        
        reshaped_Y_train = Y_train.values.reshape(Y_train.shape[0])
        new_X_train = pd.concat([pd.DataFrame(X_train_pca),pd.DataFrame(mlp_pred)], axis = 1)
        
        signals = ["Model 2: SVR done."]
        for model in [self.svr]:
            model.fit(X_train_pca,reshaped_Y_train)
            pred = model.predict(X_train_pca)
            new_X_train = pd.concat([new_X_train, pd.DataFrame(pred)], axis = 1)
            self.echo(text = signals[0], verbose = verbose)
            del signals[0]
        
        self.ridge.fit(new_X_train,reshaped_Y_train)
        self.echo(text = "Meta-Model: Ridge done.", verbose = verbose)
    
    def echo(self,text,verbose):
        if verbose not in [True,False]:
            raise ValueError('verbose must be True or False')
        if verbose:
            print(text)
        
    def predict(self, X_test):
        X_test_pca = self.pca.transform(X_test)
        mlp_pred = self.mlp.predict(X_test,verbose=0)
        new_X_test = pd.concat([pd.DataFrame(X_test_pca),pd.DataFrame(mlp_pred)], axis = 1)
        
        for model in [self.svr]:
            pred = model.predict(X_test_pca)
            new_X_test = pd.concat([new_X_test, pd.DataFrame(pred)], axis = 1)
            
        final_pred = self.ridge.predict(new_X_test)
        return final_pred
    
    def create_mlp(self):
        learning_rate = self.state[1]
        model = Sequential()
        for i in range(4,len(self.state),2):
            n_neurons = self.state[i]
            activation = self.activation_indices[self.state[i+1]]
            if n_neurons != 0:
                model.add(Dense(n_neurons,activation=activation))
        model.add(Dense(1, activation="linear"))
        model.compile(loss=self.root_mean_squared_error,optimizer=Adam(learning_rate))
        return model
    
    def root_mean_squared_error(self, y_true, y_pred):
        y_true = cast(y_true,float32)
        return kb.sqrt(kb.mean(kb.square(y_pred - y_true)))
    
    def fit_mlp(self,X_train,Y_train):
        epochs = self.state[0]
        validation_split = self.state[2]
        batch_size = self.state[3]
        self.mlp.fit(X_train,Y_train,validation_split=validation_split,
                     batch_size=batch_size,epochs=epochs,verbose=0)
        
def train(config, seed: int = 0) -> float:
    rmse = []
    rmse1 = []
    cv_path = "data/cv/kmer_splits"
    test_path = "data/test"
    for i in range(5):
        X_train = pd.read_csv(f"{cv_path}/X_train_{i}.csv")
        Y_train = pd.read_csv(f"{cv_path}/Y_train_{i}.csv")
        X_test = pd.read_csv(f"{cv_path}/X_test_{i}.csv")
        Y_test = pd.read_csv(f"{cv_path}/Y_test_{i}.csv")
        model = CopyNumberPredictor(config)
        model.fit(X_train,Y_train,verbose=False)
        pred = model.predict(X_test)
        rmse.append(sqrt(mean_squared_error(Y_test,pred)))
    X_train = pd.read_csv(f"{test_path}/X_final_train.csv")
    Y_train = pd.read_csv(f"{test_path}/Y_final_train.csv")
    X_test = pd.read_csv(f"{test_path}/X_final_test.csv")
    Y_test = pd.read_csv(f"{test_path}/Y_final_test.csv")
    for j in range(3):
        model = CopyNumberPredictor(config)
        model.fit(X_train,Y_train,verbose=False)
        pred = model.predict(X_test)
        rmse1.append(sqrt(mean_squared_error(Y_test,pred)))
    return {"Train": mean(rmse), "Test": mean(rmse1)}


cs = ConfigurationSpace()

n_neurons_0 = UniformIntegerHyperparameter("n_neurons_0", 0, 1024, default_value=490)
activation_0 = UniformIntegerHyperparameter("activation_0", 0, 4, default_value=1)
n_neurons_1 = UniformIntegerHyperparameter("n_neurons_1", 0, 1024, default_value=940)
activation_1 = UniformIntegerHyperparameter("activation_1", 0, 4, default_value=0)
n_neurons_2 = UniformIntegerHyperparameter("n_neurons_2", 0, 1024, default_value=648)
activation_2 = UniformIntegerHyperparameter("activation_2", 0, 4, default_value=0)
n_neurons_3 = UniformIntegerHyperparameter("n_neurons_3", 0, 1024, default_value=909)
activation_3 = UniformIntegerHyperparameter("activation_3", 0, 4, default_value=3)
n_neurons_4 = UniformIntegerHyperparameter("n_neurons_4", 0, 1024, default_value=590)
activation_4 = UniformIntegerHyperparameter("activation_4", 0, 4, default_value=1)
n_neurons_5 = UniformIntegerHyperparameter("n_neurons_5", 0, 1024, default_value=83)
activation_5 = UniformIntegerHyperparameter("activation_5", 0, 4, default_value=4)
batch_size = UniformIntegerHyperparameter("batch_size", 64, 128, default_value=100)

validation_split = UniformFloatHyperparameter("validation_split", 0.1, 0.33, default_value=0.11)
learning_rate = UniformFloatHyperparameter("learning_rate", 0.00005, 0.01, default_value=0.000492238598052977)

epochs = UniformIntegerHyperparameter("epochs", 20, 60, default_value=59)

alpha = UniformIntegerHyperparameter("alpha", 1, 100,default_value=22)
C = UniformIntegerHyperparameter("C", 1, 100, default_value=13)
gamma = CategoricalHyperparameter("gamma", ["auto", "scale"], default_value="auto")

cs.add_hyperparameters([epochs,learning_rate,validation_split,batch_size,n_neurons_0,activation_0,
                       n_neurons_1,activation_1,n_neurons_2,activation_2,n_neurons_3,activation_3,
                       n_neurons_4,activation_4,n_neurons_5,activation_5,alpha,C,gamma])
max_epochs = 10
scenario = Scenario(
        {
            "run_obj": "quality",  # we optimize quality (alternatively runtime)
            "runcount-limit": max_epochs,  # max. number of function evaluations
            "cs": cs,  # configuration space
            "multi_objectives": "Train, Test",
            "limit_resources": False,
        }
    )

smac = SMAC4BB(
        scenario=scenario,
        rng=np.random.RandomState(42),
        tae_runner=train
    )

# smac.solver.intensifier.tae_runner.use_pynisher = False

incumbent = smac.optimize()

print("Optimized configuration %s" % str(incumbent))

quit()