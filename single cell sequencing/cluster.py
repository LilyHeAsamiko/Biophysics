# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 14:57:42 2019

@author: Lily\he\asamiko
"""
import pandas as pd
from sklearn import preprocessing
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

path = r'D:\TUT\Medical\biophysics\experiment'
#importing the dataset
dataset = pd.read_csv(f'{path}/Data_Cortex_Nuclear.csv')
# peeking at the dataset
dataset.head().T
#Descriptive stats of the variables in data
dataset.describe()

#analysis
col = range(1,78)
dataset1 = pd.read_csv(f'{path}/Data_Cortex_Nuclear.csv', usecols = col, header =1, skiprows = 1)

#standardize the data to normal distribution
dataset1_std = preprocessing.scale(dataset1)
dt1_df = pd.DataFrame(dataset1_std.astype(np.float16))
#replace blank cells with mean
dt1_df.fillna(np.mean(np.mean(dt1_df)), inplace=True)
dt1_df.dropna(how = 'all')
# find the appropriate cluster number
plt.figure()
wcss = []
for i in range(1, 31):
    kmeans = KMeans(n_clusters = i, init = 'k-means++', random_state = 42)
    kmeans.fit(dt1_df)
    wcss.append(kmeans.inertia_)
plt.plot(range(1, 31), wcss)
plt.title('The Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.show()
plt.savefig('best_cluster_number.pdf', dpi = 200)

# Fitting K-Means to the dataset
kmeans = KMeans(n_clusters = 4, init = 'k-means++', random_state = 42)
y_kmeans = kmeans.fit_predict(dt1_df)
#beginning of  the cluster numbering with 1 instead of 0
y_kmeans1=y_kmeans+1
# New Dataframe called cluster
cluster = pd.DataFrame(y_kmeans1)
# Adding cluster to the Dataset1
dt1_df['cluster_KNN'] = cluster


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x,y = np.meshgrid(dt1_df.iloc[:,0].astype('float16'), np.linspace(min(dt1_df.iloc[0,:].astype('float16')),max(dt1_df.iloc[0,:].astype('float16')),1078))
#Mean of clusters
#kmeans_mean_cluster = pd.DataFrame(round(dt1_df.groupby('cluster_KNN').mean(),1))
#kmeans_mean_cluster
z = pd.DataFrame(round(dt1_df.groupby(['cluster_KNN']).mean()))
for cluster, zlow, zhigh in [(1, min(z.iloc[0,:]),max(z.iloc[0,:])),(2, min(z.iloc[1,:]),max(z.iloc[1,:])),(3, min(z.iloc[2,:]),max(z.iloc[2,:])),(4, min(dt1_df.iloc[3,:]),max(dt1_df.iloc[3,:]))]:
    ax.scatter(x[0,:].astype('float16'), y[:,0],  cluster, cmap='prism')  # plot points with cluster dependent colors

#plt.title('Protein Data - KNN Clustering')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plt.savefig('KNN Clustering Scatter_3D.pdf', dpi = 200)



# Hierarchical clustering for the same dataset
# creating a dataset for hierarchical clustering
#dt2_standardized = dataset1_standardized

from scipy.cluster.hierarchy import dendrogram, linkage
np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation
#creating the linkage matrix
H_cluster = linkage(dt1_df.iloc[:,0:77],'ward')
plt.title('Hierarchical Clustering Dendrogram (truncated)')
plt.xlabel('sample index or (cluster size)')
plt.ylabel('distance')
dendrogram(
    H_cluster,
    truncate_mode='lastp',  # show only the last p merged clusters
    p=4,  # show only the last p merged clusters
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,  # to get a distribution impression in truncated branches
)

plt.show()
plt.savefig('Hierarchical Clustering.pdf', dpi = 200)
# Assigning the clusters and plotting the observations as per hierarchical clustering
from scipy.cluster.hierarchy import fcluster
k=4
cluster_2 = fcluster(H_cluster, k, criterion='maxclust')
cluster_2[0:70:,]
plt.figure()
plt.scatter(dt1_df.iloc[:,0], dt1_df.iloc[:,1],c=cluster_2, cmap='prism')  # plot points with cluster dependent colors
plt.title('Protein Data - Hierarchical Clutering')
plt.show()
plt.savefig('Hierarchical Clustering Scatter_1_2.pdf', dpi = 200)
# New Dataframe called cluster
cluster_Hierarchical = pd.DataFrame(cluster_2)
# Adding the hierarchical clustering to dataset
dt1_df['cluster_HC'] = cluster_Hierarchical

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x,y = np.meshgrid(dt1_df.iloc[:,0].astype('float16'), np.linspace(min(dt1_df.iloc[0,:].astype('float16')),max(dt1_df.iloc[0,:].astype('float16')),1078))
#Mean of clusters
#kmeans_mean_cluster = pd.DataFrame(round(dt1_df.groupby('cluster_KNN').mean(),1))
#kmeans_mean_cluster
z = pd.DataFrame(round(dt1_df.groupby(['cluster_HC']).mean()))
for cluster_2, zlow, zhigh in [(1, min(z.iloc[0,:]),max(z.iloc[0,:])),(2, min(z.iloc[1,:]),max(z.iloc[1,:])),(3, min(z.iloc[2,:]),max(z.iloc[2,:])),(4, min(dt1_df.iloc[3,:]),max(dt1_df.iloc[3,:]))]:
    ax.scatter(x[0,:].astype('float16'), y[:,0],  cluster_2, cmap='prism')  # plot points with cluster dependent colors

#plt.title('Protein Data - KNN Clustering')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plt.savefig('Hierarchical Clustering Scatter_3D.pdf', dpi = 200)
