import os
import netCDF4
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from kneed import KneeLocator
from netCDF4 import Dataset
from sklearn.cluster import KMeans

f= Dataset('dTS_all_ocean.nc','r')
lats = f.variables['lat'][:]
lons = f.variables['lon'][:]
time = f.variables['time'][:]
temp = f.variables['TS'][:,:,:]
X = temp.reshape(temp.shape[0],(temp.shape[1]*temp.shape[2]))
sse = []

for i in range(1,9):
    kmeans = KMeans(n_clusters=i,init='random', random_state=0, n_init="auto")
    kmeans.fit(X)
    sse.append(kmeans.inertia_)
kl = KneeLocator(range(1,9),sse, curve="convex", direction="decreasing")
print(kl.elbow)
plt.plot(range(1,9), sse, marker='o')
plt.title('Elbow Method')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Inertia')
plt.show()
