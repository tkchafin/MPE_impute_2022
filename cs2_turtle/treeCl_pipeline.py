
import sys
import os
import treeCl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import random
from tqdm import tqdm
import numpy as np
from sklearn.cluster import SpectralClustering, AffinityPropagation


# read trees and alignments
c = treeCl.Collection(input_dir='fasta', file_format='fasta', trees_dir="trees")

# geodesic distances for tree set
dm = c.get_inter_tree_distances('geo',show_progress=False)
dm.to_csv("pairwise_geodesic.tsv", index=False, sep="/t", quoting=None)

# Spectral Clustering for K 1-10
spclust = treeCl.Spectral(dm)

partitions = list()
with open("cluster_assignments.tsv", "w") as of:
    for n in range(1, 11):
        partition = spclust.cluster(n)
        partitions.append(partition)
        ol=str(n)
        ol = ol + "\t"
        o = [str(p) for p in partition.partition_vector]
        ol = ol + '\t'.join(o) + "\n"
        of.write(ol)
        print(partition)

#
# Getting transformed coordinates
embedding = spclust.spectral_embedding(2) # spectral embedding in 2 dimensions
embedding.df.to_csv("embedding.tsv", index=False, sep="\t", quoting=None)

# Estimate best K using eigen gap heuristic 
# Code from github.com/ciortanmadalina/high_noise_clustering
from scipy.spatial.distance import pdist, squareform
def getAffinityMatrix(coordinates, k = 7):
    """
    Calculate affinity matrix based on input coordinates matrix and the numeber
    of nearest neighbours.
    
    Apply local scaling based on the k nearest neighbour
        References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    """
    # calculate euclidian distance matrix
    dists = squareform(pdist(coordinates)) 
    
    # for each row, sort the distances ascendingly and take the index of the 
    #k-th position (nearest neighbour)
    knn_distances = np.sort(dists, axis=0)[k]
    knn_distances = knn_distances[np.newaxis].T
    
    # calculate sigma_i * sigma_j
    local_scale = knn_distances.dot(knn_distances.T)

    affinity_matrix = dists * dists
    affinity_matrix = -affinity_matrix / local_scale
    # divide square distance matrix by local scale
    affinity_matrix[np.where(np.isnan(affinity_matrix))] = 0.0
    # apply exponential
    affinity_matrix = np.exp(affinity_matrix)
    np.fill_diagonal(affinity_matrix, 0)
    return affinity_matrix
import scipy
from scipy.sparse import csgraph
# from scipy.sparse.linalg import eigsh
from numpy import linalg as LA
def eigenDecomposition(A, plot = True, topK = 5):
    """
    :param A: Affinity matrix
    :param plot: plots the sorted eigen values for visual inspection
    :return A tuple containing:
    - the optimal number of clusters by eigengap heuristic
    - all eigen values
    - all eigen vectors
    
    This method performs the eigen decomposition on a given affinity matrix,
    following the steps recommended in the paper:
    1. Construct the normalized affinity matrix: L = D−1/2ADˆ −1/2.
    2. Find the eigenvalues and their associated eigen vectors
    3. Identify the maximum gap which corresponds to the number of clusters
    by eigengap heuristic
    
    References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf
    """
    L = csgraph.laplacian(A, normed=True)
    n_components = A.shape[0]
    
    # LM parameter : Eigenvalues with largest magnitude (eigs, eigsh), that is, largest eigenvalues in 
    # the euclidean norm of complex numbers.
#     eigenvalues, eigenvectors = eigsh(L, k=n_components, which="LM", sigma=1.0, maxiter=5000)
    eigenvalues, eigenvectors = LA.eig(L)
    
    if plot:
        plt.title('Largest eigen values of input matrix')
        plt.scatter(np.arange(len(eigenvalues)), eigenvalues)
        plt.grid()
        
    # Identify the optimal number of clusters as the index corresponding
    # to the larger gap between eigen values
    index_largest_gap = np.argsort(np.diff(eigenvalues))[::-1][:topK]
    nb_clusters = index_largest_gap + 1
        
    return nb_clusters, eigenvalues, eigenvectors
 
#affinity_matrix = getAffinityMatrix(embedding.df, k = 10)
k, eigenvals,  eigenvecs = eigenDecomposition(spclust._affinity, topK=10)
print(f'Optimal number of clusters {k}')
#print(eigenvecs)

bestk=pd.DataFrame(k, columns=["K"])
bestk.to_csv("bestK.tsv", sep="\t", index=False, quoting=None)

eig=pd.DataFrame({'K':list(range(len(eigenvals))), 'Eig':eigenvals})
eig.to_csv("eigengap.tsv", index=False, sep="\t", quoting=None)

raxml = treeCl.tasks.RaxmlTaskInterface()
part_names=dict()
k=1
for p in partitions:
    sc=treeCl.Scorer(c, cache_dir='scorer', task_interface=raxml)
    loci = sc.get_partition_members(p)
    if "Name" not in part_names:
        names=[i for sublist in loci for i in sublist]
        part_names["Name"]=names
    locdict=dict()
    clust=1
    for sublist in loci:
        for i in sublist:
            locdict[i]=str(clust)
        clust=clust+1
    labels=list()
    for name in part_names["Name"]:
        labels.append(locdict[name])
    part_names[str(k)]=labels
    k=k+1
labdf=pd.DataFrame(part_names)
labdf.to_csv("labeled_partitions.tsv", index=False, sep="\t", quoting=None)

print("Done!")
