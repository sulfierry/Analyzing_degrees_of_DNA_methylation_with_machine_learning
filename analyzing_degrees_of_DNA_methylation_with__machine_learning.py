#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import statistics as st


# In[2]:


df_original = pd.read_csv('./data_3_methylation_medina_class.tsv', '\t')
df_original2 = pd.read_csv('./data_3_methylation_medina_v2.txt', '\t', usecols=list(df_original.columns))
df_original.info()


# In[3]:


#df = df_original.drop(['chr','start', 'end', 'gene_name'], axis=1).copy()

df_original2 = df_original2.dropna()
df = df_original2.drop(['chr','start', 'end', 'gene_name'], axis=1).copy()


# In[4]:


#l_col = df_original2.columns
#for i in l_col[4:]:
#    df_original2[i][df_original2[df_original2[i] >= 0.65][i].index] = 2
#    df_original2[i][df_original2[(df_original2[i] > 0.35) & (df_original2[i] < 0.65)][i].index] = 1
#    df_original2[i][df_original2[df_original2[i] <= 0.35][i].index] = 0


# In[5]:


rets3 = df[:].pct_change()
corr = rets3.corr()  
plt.figure(figsize=(12,12))
plt.imshow(corr, cmap='hot', interpolation='none')  
plt.colorbar()  
plt.xticks(range(len(corr)), corr.columns, rotation=90, size=15)  
plt.yticks(range(len(corr)), corr.columns, size=15);
plt.show()


# In[6]:


#used = []
#corrs = []
#for i, j in enumerate(df.corr().columns):
#    for k in range(len(df.corr())):
#        if ((j != df.corr().index[k]) & (abs(df.corr().iloc[k, i]) > .7)):
#            used.append(j)
#            corrs.append((j, df.corr().index[k],
#                          np.round(df.corr().iloc[k, i], 2)))
#corrs


# In[7]:


from sklearn.cluster import KMeans

numClusters = [1,2,3,4,5,6,7,8,9,10]
SSE = []
for k in numClusters:
    k_means = KMeans(n_clusters=k)
    k_means.fit(df)
    SSE.append(k_means.inertia_)

plt.plot(numClusters, SSE)
plt.xlabel('Number of Clusters')
plt.ylabel('SSE')


# In[8]:


kmeans = KMeans(n_clusters=3).fit(df)
centroids = kmeans.cluster_centers_


# In[9]:


df_original2['GRUPO'] = kmeans.labels_

def f(x):
    if x == 0:
        return 'A'
    if x == 1:
        return 'B'
    if x == 2:
        return 'C'
    if x == 3:
        return 'D'
    if x == 4:
        return 'E'

df_original2['GRUPO'] = df_original2['GRUPO'].map(lambda x: f(x))

A = df_original2[df_original2['GRUPO'] == 'A'].copy()
B = df_original2[df_original2['GRUPO'] == 'B'].copy()
C = df_original2[df_original2['GRUPO'] == 'C'].copy()
D = df_original2[df_original2['GRUPO'] == 'D'].copy()
E = df_original2[df_original2['GRUPO'] == 'E'].copy()


# In[10]:


df_original2['GRUPO'].value_counts()


# In[11]:


#for i in A.columns:
#    print('A')
#    print(A[i].value_counts())
#    print('B')
#    print(B[i].value_counts())
#    print('C')
#    print(C[i].value_counts())
#    print('D')
#    print(D[i].value_counts())
#    print('E')
#    print(E[i].value_counts())


# In[12]:


l_hist = list(A['HSPC'])
for i in A.columns[5:-1]:
    l_hist = l_hist + list(A[i])

print(st.mean(l_hist))
print(st.stdev(l_hist))

plt.figure(figsize=(10,6))
plt.hist(l_hist, bins=80)
plt.show()


# In[13]:


l_hist = list(B['HSPC'])
for i in B.columns[5:-1]:
    l_hist = l_hist + list(B[i])

print(st.mean(l_hist))
print(st.stdev(l_hist))

plt.figure(figsize=(10,6))
plt.hist(l_hist, bins=80)
plt.show()


# In[14]:


l_hist = list(C['HSPC'])
for i in C.columns[5:-1]:
    l_hist = l_hist + list(C[i])

print(st.mean(l_hist))
print(st.stdev(l_hist))

plt.figure(figsize=(10,6))
plt.hist(l_hist, bins=80)
plt.show()


# In[ ]:


l_hist = list(D['HSPC'])
for i in D.columns[5:-1]:
    l_hist = l_hist + list(D[i])

print(st.mean(l_hist))
print(st.stdev(l_hist))

plt.figure(figsize=(10,6))
plt.hist(l_hist, bins=80)
plt.show()


# In[ ]:


l_hist = list(E['HSPC'])
for i in E.columns[5:-1]:
    l_hist = l_hist + list(E[i])

print(st.mean(l_hist))
print(st.stdev(l_hist))

plt.figure(figsize=(10,6))
plt.hist(l_hist, bins=80)
plt.show()


# In[23]:


for i in C['gene_name']:
    print(i)


# In[26]:


C[C['gene_name'] == 'PARD6G-AS1']


# In[28]:


A[A['gene_name'] == 'WRB']


# In[32]:


A[A['gene_name'] == 'LINC00319']


# In[ ]:


A.to_csv('alto.csv', index=False)
B.to_csv('baixo.csv', index=False)
C.to_csv('medio.csv', index=False)

