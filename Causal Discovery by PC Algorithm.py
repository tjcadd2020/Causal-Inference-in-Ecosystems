#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 5 20:06:07 2023

@author: gaosheng
"""


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes

import networkx as nx
from cdt.causality.graph import PC
from cdt.data import load_dataset


# graph
CD_edge = pd.read_csv('/home/data/CD_edge.csv')
CD_node = pd.read_csv('/home/data/CD_node.csv')

edges = [(list(CD_edge['Source'])[i],list(CD_edge['Target'])[i],{'weight':list(CD_edge['Weight'])[i],'pos':list(CD_edge['Pos_Neg'])[i],'association':list(CD_e
dge['Association'])[i],'FDR':list(CD_edge['FDR'])[i]}) for i in range(CD_edge.shape[0])]

#生成图
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes

import networkx as nx
from cdt.causality.graph import PC
from cdt.data import load_dataset


# graph
CD_edge = pd.read_csv('/home/data/CD_edge.csv')
CD_node = pd.read_csv('/home/data/CD_node.csv')

edges = [(list(CD_edge['Source'])[i],list(CD_edge['Target'])[i],{'weight':list(CD_edge['Weight'])[i],'pos':list(CD_edge['Pos_Neg'])[i],'association':list(CD_e
dge['Association'])[i],'FDR':list(CD_edge['FDR'])[i]}) for i in range(CD_edge.shape[0])]

#生成图
G_CD = nx.Graph()                    #无向图
G_edges = edges
G_CD.add_edges_from(G_edges)
labels={}
for node in G_CD.nodes():
    labels[node]=str(node)

node_dic = {}
for i in CD_node.index:
    name = list(CD_node['Id'])[i]
    ty = list(CD_node['Type'])[i]
    label = list(CD_node['Label'])[i]
    node_dic[name] = {'name':name,'type':ty,'label':label}

for i in G_CD._node.keys():
    if i in node_dic.keys():
        G_CD._node[i] = node_dic[i]



# data
data_CD = pd.read_csv('/home/data/sp_marker_CD.tsv',sep = '\t',index_col=0)
data_CD = data_CD.fillna(0)
data_CD = data_CD.loc[data_CD.index.isin(list(G_CD.nodes)),:]
data_CD = data_CD.loc[:,list(data_CD.sum(0)!=0)]


#定义对象
obj = PC()
output_CD = obj.predict(data_CD.T, G_CD)  #With an undirected graph


import pickle
f=open('/home/outputs/output_CD.pkl','wb')
pickle.dump(output_CD,f)
f.close()

nx.write_gexf(output_CD,'/home/outputs/output_CD.gexf')
