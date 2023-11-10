from biclusteringVisualization.plotCreator import visualizeClustering

"""
Advanced demonstration, hte clustering is given as a python list

"""

inputMatrix = [
    [1, 1, 0, 0, 0, 1, 1, 0, 0, 1],
    [1, 1, 0, 0, 0, 1, 0, 0, 0, 1],
    [0, 0, 1, 1, 1, 1, 1, 0, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 0, 1, 0],
    [0, 1, 0, 0, 0, 1, 1, 0, 1, 1],
]
rowCluster = [[0, 1, 2, 4], [2, 3]]
columnCluster = [[0, 1, 5, 6, 9], [2, 3, 4, 5, 6, 8]]

outputFile = "advanced2.pdf"
orderingMethod = "TSPheuristic"
heuristicTime = 4
plotUnorderedMatrix = True

visualizeClustering(
    inputFile=inputMatrix,
    outputFile=outputFile,
    rowClusters=rowCluster,
    columnClusters=columnCluster,
    orderingMethod=orderingMethod,
    heuristicTime=heuristicTime,
    plotUnorderedMatrix=plotUnorderedMatrix,
    printPermutations=False,
    customColorScheme=True,
)
