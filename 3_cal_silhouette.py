import os
import netCDF4
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.cluster import KMeans

f= Dataset('dTS_all_ocean.nc','r')
lats = f.variables['lat'][:]
lons = f.variables['lon'][:]
time = f.variables['time'][:]
temp = f.variables['TS'][:,:,:]

SSE = []
sil_coef = []

for i in range(2,9):
    kmeans = KMeans(n_clusters=i,init='random', random_state=0, n_init="auto")
    cluster_labels = kmeans.fit_predict(temp[i,:,:])
    score = sklearn.metrics.silhouette_score(temp[i,:,:],kmeans.labels_,metric='euclidean')
    sil_coef.append(score)

X = range(2,9)
plt.plot(X, sil_coef, marker='o')
plt.title('Silhouette Score')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Silhouette Coefficient')
plt.show()
