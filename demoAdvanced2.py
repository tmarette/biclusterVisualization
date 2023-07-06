from biclusteringVisualization.plotCreator import visualizeClustering

"""
Advanced demonstration, hte clustering is given as a python list

"""

inputFile = "demo/advanced.csv"
rowCluster = [[0, 1, 2, 4], [1, 2, 3]]
columnCluster = [[0, 1, 5, 6, 8], [2, 3, 4, 5, 6, 9]]
inputIsList = True


outputFile = "advanced2.pdf"
orderingMethod = "TSPheuristic"
heuristicTime = 4
plotUnorderedMatrix = True

visualizeClustering(
    inputFile=inputFile,
    outputFile=outputFile,
    rowClusters=rowCluster,
    columnClusters=columnCluster,
    inputIsList=inputIsList,
    orderingMethod=orderingMethod,
    heuristicTime=heuristicTime,
    plotUnorderedMatrix=plotUnorderedMatrix,
)
