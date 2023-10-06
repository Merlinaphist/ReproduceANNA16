import pandas as pd
import numpy as np
import pickle, os, shutil, shap
from math import sqrt
from zipfile import ZipFile
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow import cast,float32
from keras import backend as kb
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Ridge, Lasso
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

class Preprocessing():
    def __init__(self,k_size=6):
        self.k_size = k_size
        kmers = self.generate_kmers("",self.k_size)
        self.vectorizer = CountVectorizer(vocabulary = kmers)
        self.seqs = []
    
    def generate_kmers(self,current_kmer,current_depth):
        if current_depth == 1:
            return [current_kmer+"a",current_kmer+"t",current_kmer+"c",current_kmer+"g"]
        else:
            ret = self.generate_kmers(current_kmer+"a",current_depth-1)
            for nt in ['t','c','g']:
                ret += self.generate_kmers(current_kmer+nt,current_depth-1)
            return ret
    
    def generate_kmer_multiple(self,seqlist,k):
        kmer_list = []
        n = -1
        for seq in seqlist:
            kmer_list.append(self.generate_kmer_single(str(seq),k))
        return kmer_list
    
    def generate_kmer_single(self,seq,k):
        kmer = ""
        for i in range(0,len(seq)-k,1):
            kmer += seq[i:i+k]+" "
        return kmer[:-1]
    
    def CountKmers(self,seqs):
        if type(seqs) in [type([]),type(pd.core.series.Series([1]))]:
            kmer = self.generate_kmer_multiple(seqs, self.k_size)
            transformed_X = self.vectorizer.transform(kmer).toarray()
            return transformed_X
        else:
            raise ValueError("""Invalid 'seqs' format.
            Expected formats are 'list' or 'pandas.core.series.Series'.""")
            
    def ReadFASTA(self,filename,as_pd=True):
        if filename.split(".")[-1] not in ["fasta","fna","fa"]:
            raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
        file_handle = open(filename,"r")
        seqs = []
        seqid = []
        tmp_seq = ""
        for line in file_handle:
            if (line[0] == ">"):
                if tmp_seq != "":
                    seqs.append(tmp_seq)
                seqid.append(line.split("\n")[0][1:])
                tmp_seq = ""
            else:
                tmp_seq+=line.split("\n")[0]
        seqs.append(tmp_seq)
        file_handle.close()
        if as_pd:
            fasta = {}
            for i in range(len(seqs)):
                fasta[seqid[i]] = seqs[i]
            return pd.DataFrame(fasta,index=["sequence"]).transpose()["sequence"]
        else:
            return seqs, seqid

class CopyNumberPredictor():
    def __init__(self,region):
        self.region = region
        self.state = [59, 0.00016770313599, 0.11, 100, 489, 1, 926, 0, 645, 0, 929, 3, 582, 1, 82, 4]
        self.activation_indices={0:"relu",1:"gelu",2:"selu",3:"elu",4:"linear"}
        self.mlp = self.create_mlp()
        self.ridge = Ridge(alpha = 49)
        self.pca = PCA(n_components=100)
        self.svr = SVR(kernel='rbf',C=11,gamma='auto')
        
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
    

da0 = pd.read_csv("datasets/full_length_reads.csv")
# da1 = pd.read_csv("4.final_test/datasets/full_length_testdata.csv")
da1 = pd.read_csv("datasets/full_length_testdata.csv")
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