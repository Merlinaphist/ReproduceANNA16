import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer

class Preprocessing():
    def __init__(self, k_size=6):
        self.k_size = k_size
        kmers = self.ref_kmers("", self.k_size)
        self.vectorizer = CountVectorizer(vocabulary = kmers)
        self.seqs = []

    def ref_kmers(self, current_kmer, current_depth):
        if current_depth == 1:
            return [current_kmer+"a",current_kmer+"u",current_kmer+"c",current_kmer+"g"]
        else:
            ret = self.ref_kmers(current_kmer+"a",current_depth-1)
            for nt in ['u','c','g']:
                ret += self.ref_kmers(current_kmer+nt,current_depth-1)
            return ret

    def seq2kmer(self, seq, k):
        kmer = ""
        for i in range(0,len(seq)-k,1):
            kmer += seq[i:i+k]+" "
        return kmer[:-1]

    def CountKmers(self,seqs):
        if type(seqs) in [type([]),type(pd.core.series.Series([1]))]:
            kmer = pd.Series(seqs).apply(lambda x: self.seq2kmer(x, self.k_size))
            transformed_X = self.vectorizer.transform(kmer).toarray()
            return transformed_X
        else:
            raise ValueError("Invalid 'seqs' format. Expected formats are 'list' or 'pandas.core.series.Series'.")

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
        
class Tree():
    def __init__(self,name,rank):
        self.subtree = []
        self.mean = 0
        self.name = name
        self.rank = rank
        self.pan_stats = []

    def enter(self,lineage):
        if len(lineage) == 1:
            self.mean = int(lineage[0])
            return
        
        result = self.find(lineage[0])
        if result == None:
            new_node = Tree(lineage[0],self.rank+1)
            lineage = lineage[1:]
            new_node.enter(lineage)
            self.subtree.append(new_node)
        else:
            lineage = lineage[1:]
            self.subtree[result].enter(lineage)
        self.averaging()
    
    def find(self,target_name):
        for i in range(0,len(self.subtree),1):
            if self.subtree[i].name == target_name:
                return i
        return None
    
    def averaging(self):
        add = 0
        for i in range(0,len(self.subtree),1):
            add += self.subtree[i].mean
        self.mean = add/len(self.subtree)

    def output(self):
        pan_stats = []
        pan_stats = self.pull(pan_stats)
        return pan_stats
        
    def pull(self,pan_stats):
        if str(self.name) != "nan":
            pan_stats.append([self.rank,str(self.name),self.mean])
        if self.subtree == []:
            return pan_stats
        for i in range(0,len(self.subtree),1):
            pan_stats = self.subtree[i].pull(pan_stats)
        return pan_stats
    
class GCNlibrary():
    def __init__(self):
        self.lib = {}

    def record(self, rank, name, mean):
        if rank in list(self.lib.keys()):
            self.lib[rank][name] = mean
        else:
            rank_dict = {}
            rank_dict[name] = mean
            self.lib[rank] = rank_dict
    
    def search(self, index, name):
        return self.lib.get(index).get(name,"N/A")

    def fit(self,pan_stats):
        for i in range(0,len(pan_stats),1):
            rank = pan_stats[i][0]
            name = pan_stats[i][1]
            mean = float(pan_stats[i][2])
            self.record(rank,name,mean)
            
    def predict(self,taxa):
        ultimean = self.lib.get(0).get('prokaryotes')
        output = []
        for lineage in taxa:
            for index in range(len(lineage),0,-1):
                result = self.search(index,lineage[index-1])
                if result == 'N/A':
                    if index == 1:
                        output.append(ultimean)
                        break
                    else:
                        continue
                else:
                    output.append(result)
                    # lineage.append(search)
                    break
        return output
    
class TaxAvg():
    def __init__(self):
        self.prokaryotes = Tree("prokaryotes",0)
        self.model = GCNlibrary()
        self.pan_stats = []
        
    def fit(self,da):
        for lineage in da:
            self.prokaryotes.enter(lineage)
        self.pan_stats = self.prokaryotes.output()
        self.model.fit(self.pan_stats)
        
    def predict(self,X):
        pred = self.model.predict(X)
        return pred

class CopyRighterSimulator():
    def __init__(self,CopyRighterData):
        self.lib = {}
        for line in CopyRighterData:
            self.lib[line[0]] = line[1]

    def predict(self,X):
        pred = []
        for line in X:
            key=line
            value = self.lib.get(key,"N/A")
            taxalist = line.split(";")
            while value == "N/A":
                taxalist = taxalist[:-1]
                if len(taxalist) == 0:
                    value = 2.3539721350613916
                else:
                    key = self.combine(taxalist)
                    value = self.lib.get(key,"N/A")
            pred.append(value)
        return pred
                
    def combine(self,taxalist):
        string = ""
        for element in taxalist:
            string+=element
            string+=";"
        return string[:-1]

class rrnDBSimulator():
    def __init__(self,rrnDBData):
        self.lib = {"domain":{},"phylum":{},"class":{},"order":{},"family":{},"genus":{}}
        self.ranks = ["domain","phylum","class","order","family","genus"]
        for line in rrnDBData:
            if line[0] in self.ranks:
                self.lib[line[0]][line[1]] = line[2]

    def predict(self,X):
        pred = []
        for line in X:
            key=line
            ranks = self.ranks
            value = self.lib.get(ranks[-1],{}).get(key[-1],"N/A")
            while value == "N/A":
                ranks = ranks[:-1]
                key = key[:-1]
                if len(key) == 0:
                    value = 2.3539721350613916
                else:
                    value = self.lib.get(ranks[-1],{}).get(key[-1],"N/A")
            pred.append(value)
        return pred