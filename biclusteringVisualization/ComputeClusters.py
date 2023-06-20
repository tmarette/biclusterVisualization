import os
import sys
import tempfile
import scipy.io
from scipy.sparse import csr_matrix
from datetime import datetime

sys.path.append("includes/pcv")
import PCV


def cluster(A, k, algo):
    leftClusters = []
    rightClusters = []
    if algo == "PCV":
        [leftClusters, rightClusters] = pcv(
            A, k, 0.5, minClusterSize=1, useHeuristic=True
        )
    elif algo == "basso":
        [leftClusters, rightClusters] = basso(A, k)
    else:
        print("Unknown clustering algorithm method.")

    return [leftClusters, rightClusters]


def pcv(A, k, threshold, minClusterSize, useHeuristic):
    [leftClusters, rightClusters] = PCV.PCV(
        A, k, 0.5, minClusterSize=1, useHeuristic=True
    )
    return [leftClusters, rightClusters]


def basso(A, k):
    pathToBasso = "includes/basso-0.5/src/cmdline/basso"

    tmpDir = tempfile.gettempdir()
    now = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")

    Asparse = csr_matrix(A)
    mmFileA = f"{tmpDir}/matrix-{now}.mtx"
    scipy.io.mmwrite(mmFileA, Asparse)

    clustersFilePrefix = f"{tmpDir}/factor-{now}"

    command = f"{pathToBasso} -k{k} -w1 -z8 -t0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 -o {clustersFilePrefix} {mmFileA}"

    os.system(command)

    leftClustersFile = f"{clustersFilePrefix}_L.mtx"
    rightClustersFile = f"{clustersFilePrefix}_R.mtx"
    leftClusters, rightClusters = MMtoClusters(leftClustersFile, rightClustersFile)

    os.remove(mmFileA)
    os.remove(leftClustersFile)
    os.remove(rightClustersFile)

    return leftClusters, rightClusters


# convert MM files from BMF into clusters
def MMtoClusters(leftClustersFile, rightClustersFile):
    leftClustersMatrix = scipy.io.mmread(leftClustersFile).tocsr()
    rightClustersMatrix = scipy.io.mmread(rightClustersFile).tocsr()
    n, k = leftClustersMatrix.shape
    m = rightClustersMatrix.shape[1]
    leftClusters = []
    rightClusters = []
    for i in range(k):
        cluster = []
        for j in range(n):
            if leftClustersMatrix[j, i] != 0:
                cluster.append(j)
        leftClusters.append(cluster)
        cluster = []
        for j in range(m):
            if rightClustersMatrix[i, j] != 0:
                cluster.append(j)
        rightClusters.append(cluster)
    return leftClusters, rightClusters
