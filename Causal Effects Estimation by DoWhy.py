#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 5 20:06:07 2023

@author: gaosheng
"""


import numpy as np
import pandas as pd

import dowhy
from dowhy import CausalModel
import dowhy.datasets



# =============================================================================
# Causal Inference
# =============================================================================
def write_dot(graph_edge): ### convert edge csv into .dot
    dot_str = 'digraph {'
    for i in range(graph_edge.shape[0]):
        edge = ' ' + graph_edge.iloc[i][0] + '->' + graph_edge.iloc[i][1] + ';'
        dot_str = dot_str + edge
    #dot_str = dot_str + ' Unobserved Confounders->Condition;}'
    dot_str = dot_str + '}'
    return dot_str

def write_gml(graph_edge): ### convert edge csv into .dot
     dot_str = 'graph[directed 1node[ id "Condition" label "Condition"]node[ id "Unobserved Confounders" label "Unobserved Confounders"]edge[source "Unobserved Confounders" target "Condition"]'
     for i in range(graph_edge.shape[0]):
         edge = ' edge[ source "' + graph_edge.iloc[i][0] + '" target "' + graph_edge.iloc[i][1] + '"]'
         dot_str = dot_str + edge
         if graph_edge['Target'].iloc[i] == 'Condition':
                 node = 'node[ id "' + graph_edge.iloc[i][0]  + '" label "' + graph_edge.iloc[i][0]  + '"]'
                 edge = ' edge[source "Unobserved Confounders"' + ' target "' + graph_edge.iloc[i][0] + '"]'
                 dot_str = dot_str + node + edge
     dot_str = dot_str + ']'
     return dot_str




graph = pd.read_csv('dowhy_graph_input_CD.csv',index_col = 0)
graph_dot = write_dot(graph)
graph_gml = write_gml(graph)

abundance_df = pd.read_csv('dowhy_abundance_input_CD.csv',index_col = 0)
nodes = list(abundance_df.columns)

# print(graph_dot)


# =============================================================================
# Significance Evaluation of Causal Inference
# =============================================================================

ci_df_temp = abundance_df.copy()


import random

class Causal_Inference_Modeling:
    def __init__(self,causal_inference_df,dot,Source,Target):
        self.causal_inference_df = causal_inference_df
        self.dot = dot
        self.index = Source
        self.outcome = self.causal_inference_df[Target]


    def Causal_Inference_class(self):   ######## causal inference process
        model = CausalModel(
            data = self.causal_inference_df,
            treatment = Source,
            outcome = Target,
            graph = self.dot
            )

        identified_estimand = model.identify_effect() #proceed_when_unidentifiable=True ; method_name="id-algorithm"

        causal_estimate = model.estimate_effect(identified_estimand,method_name="backdoor.linear_regression") ##  LogisitcRegression

        value = causal_estimate.value
        # print(value)

        expr = causal_estimate.realized_estimand_expr   ##### extract 构造表达式，用于置换检验

        ######## Permutation
        random_cev = self.permutation_func(expr)
        random_cev.append(value)
        random_cev = pd.Series(random_cev)

        return random_cev


    def permutation_func(self,expr):  ### 蒙特卡洛方法随机999次

        random_cev = []

        for i in range(999):
            outcome = np.random.permutation(self.outcome) #### 随机最后一列disease状态值
   
            # 解析表达式，抽提source features
            temp_feature = [item.split('*') for item in expr.split('~')[1].split('+')]
            temp = []
            for item in temp_feature:
                temp += item
            temp = list(set(temp))

            features = self.causal_inference_df[temp] # 构建source feature的表达矩阵

            from sklearn import linear_model
            model = linear_model.LinearRegression()
            model.fit(features, outcome) # 线性拟合
            coefficients = model.coef_ 

            value = coefficients[0] # 拟合系数
            random_cev.append(value)
        return random_cev

class Significance_Evaluation:
    def __init__(self,CEV_df,effect_relation):
        self.CEV_df = CEV_df
        self.Effect_Relation = effect_relation

    def compute_significance(self,CEV_series): #### 根据每个series计算pvalue
        CEV_list = list(CEV_series)
        bias = [item for item in CEV_list if item >= CEV_list[-1]]
        sig = len(bias)/1000
        return sig

    def applyfunc(self):
        sig_df = self.CEV_df.apply(self.compute_significance,axis = 'rows') ### 根据每列,apply函数计算所有gene的pvalue
        combined_df = pd.DataFrame({'Effect_Relation':self.Effect_Relation, 'Causal_Estimated_Value':list(self.CEV_df.iloc[-1,:]), 'Pvalue':list(sig_df)} )

        return combined_df

import time
time1 = time.time()


######### permutation in CD
import time
time1 = time.time()


CEV_list = [] # save causal effect values
effect_relation = [] # causal relations



# compute causal effects
for i in range(graph.shape[0]):
    # for i in range(1):
    print(i)
    Source = graph['Source'][i]
    Target = graph['Target'][i]

    effect_relation.append(str(Source)+'->'+str(Target))
    ci_obj = Causal_Inference_Modeling(ci_df_temp,graph_dot,Source,Target)
    cev = ci_obj.Causal_Inference_class()
    CEV_list.append(cev)
CEV_df = pd.DataFrame(CEV_list,index = effect_relation).T

# compute pvalues
sig_obj = Significance_Evaluation(CEV_df,effect_relation)
combined_df = sig_obj.applyfunc()


time2 = time.time()
time = (time2 - time1)/3600
print('\n\n\n it has cost %s Hrs \n\n\n' %time)
combined_df.to_csv('causal_inference_df.csv',header = True, index = True)
