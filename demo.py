from biclusteringVisualization.plotCreator import visualizeClustering

"""
Basic demonstration, easy to use

"""


inputFile = "demo/example1.csv"
outputFile = "example1.pdf"

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
