from biclusteringVisualization.plotCreator import visualizeClustering

"""
Advanced demonstration, the clustering is given in csv files

"""

inputFile = "demo/advanced.csv"
rowCluster = "demo/advanced_rowClusters.csv"
columnCluster = "demo/advanced_columnClusters.csv"

outputFile = "plots/demo/advanced.pdf"
orderingMethod = "TSPheuristic"
heuristicTime = 4
plotUnorderedMatrix = True
showObjectiveFunctionValues = True

visualizeClustering(
    inputFile=inputFile,
    outputFile=outputFile,
    rowClusters=rowCluster,
    columnClusters=columnCluster,
    orderingMethod=orderingMethod,
    heuristicTime=heuristicTime,
    clusteringAlgorithm="PCV",
    plotUnorderedMatrix=plotUnorderedMatrix,
    showObjectiveFunctionValues=showObjectiveFunctionValues,
)
