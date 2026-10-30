from biclusteringVisualization.plotCreator import visualizeClustering

"""
Basic demonstration, easy to use

"""


inputFile = "demo/basic.csv"
outputFile = "plots/demo/basic.pdf"

visualizeClustering(
    inputFile=inputFile,
    outputFile=outputFile,
)

"""
The above command is equivalent to

visualizeClustering(
    inputFile=inputFile,
    outputFile=outputFile,
    clusteringAlgorithm="PCV",
    nbClusters=5,
)

"""
