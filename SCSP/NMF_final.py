# -*- coding: utf-8 -*-
# @Time    : 2023/4/10  20:29
# @Author  : Gou Yujie
# @File    : NMF_final.py
import pandas as pd
import numpy as np
def matrix_factorisation_getQ(R, P, Q, alpha=0.00000018, beta=0.00002, steps=20):
    R = np.array(R)
    Q = np.array(Q)
    P = np.array(P)
    print(R.shape,P.shape,Q.shape)
    for step in range(steps):
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                eij = R[i][j] - np.dot(P[i, :], Q[:, j])
                for k in range(P.shape[1]):
                    plus = alpha * (2 * eij * P[i][k] - beta * Q[k][j])
                    Q[k][j] = Q[k][j] + plus
    Q[Q < 0] = 0
    return Q

R=pd.read_csv("R.txt",sep='\t',index_col=0)
Q=pd.read_csv("Q.txt",sep='\t',index_col=0)
P=pd.read_csv("P.txt",sep='\t',index_col=0)
single_frame=pd.read_csv("single_plus.txt",sep='\t',index_col=0)
start=pd.read_csv("random_start_next.txt",sep='\t',index_col=0)
newQ=matrix_factorisation_getQ(R, P, start)
newQ=pd.DataFrame(newQ,index=Q.index,columns=Q.columns)
out=pd.concat([single_frame,newQ.T])
out.to_csv("getQ_SCC.txt",sep='\t')