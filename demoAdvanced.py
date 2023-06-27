from biclusteringVisualization.plotCreator import visualizeClustering

"""
Advanced demonstration

"""

inputFile = "demo/example2.csv"
rowCluster = "demo/example2_rowClusters.csv"
columnCluster = "demo/example2_columnClusters.csv"

outputFile = "example2.pdf"
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
    plotUnorderedMatrix=plotUnorderedMatrix,
    showObjectiveFunctionValues=showObjectiveFunctionValues,
)
