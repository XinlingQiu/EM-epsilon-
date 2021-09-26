import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt

mat = loadmat('D:/Python/ex7faces.mat')
X = mat['X']
print(X.shape)  #(5000, 1024)


def displayData(X, row, col):
    fig, axs = plt.subplots(row, col, figsize=(8,8))
    for r in range(row):
        for c in range(col):
            axs[r][c].imshow(X[r*col + c].reshape(32,32).T, cmap = 'Greys_r')
            axs[r][c].set_xticks([])
            axs[r][c].set_yticks([])
            
displayData(X, 10, 10)

def featureNormalize(X):
    means = X.mean(axis=0)
    stds = X.std(axis=0, ddof=1)
    X_norm = (X - means) / stds
    return X_norm, means, stds


def pca(X):
    sigma = (X.T @ X) / len(X)
    U, S, V = np.linalg.svd(sigma)
    
    return U, S, V    

X_norm, means, stds = featureNormalize(X)
U, S, V = pca(X_norm)    
#print(U.shape, S.shape)  #(1024, 1024) (1024,)

displayData(U[:,:36].T, 6, 6)

#Dimensionality Reduction

def projectData(X, U, K):
    Z = X @ U[:,:K]
    
    return Z

z = projectData(X_norm, U, K=36)
def recoverData(Z, U, K):
    X_rec = Z @ U[:,:K].T
    
    return X_rec
X_rec = recoverData(z, U, K=36)
displayData(X_rec, 10, 10)