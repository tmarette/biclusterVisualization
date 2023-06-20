# Visualizing Overlapping Biclusterings and Boolean Matrix Factorizations

README updated on the 20/06/2023

This archive contains the code to visualize biclusters from the paper (to be published)


If you want to use basso clustering algorithm, extract `include/basso-0.5.tar.gz` and follow the setup instructions

## Execute the code

The code is available in the folder `biclusteringVisualization/`
Code examples are available in `demo.py` and `demoAdvanced.py`

## Documentation

```python
visualizeClustering.visualizeClustering(
    inputFile,
    outputFile,
    clusteringAlgorithm="PCV",
    nbClusters="5",
    orderingMethod="TSPheuristic", 
    heuristicTime="20",
    plotUnorderedMatrix=False,
    showObjectiveFunctionValues="False"
    )
```
### Parameters:
- `inputFile: string`
  
   Location of the input file
- `outputFile: string`
  
  Location of the output file 
- `clusteringAlgorithm: string`
 
  Biclustering algorithm to use. `PCV` (default) and `basso` are supported.
- `nbClusters: int`
 
   Number of clusters for the clustering algorithm. Default value is 5.
- `orderingMethod: string`
 
    Algorithm to use for the ordering of the matrix in the visualization. Default value is TSPHeuristic. ADVISER, as well as the greedy algorithms used for the experiments in the paper are also available.

- `heuristicTime: int`
 
    Execution time (in seconds) for TSPHeuristic. Default value is 20 seconds.
- `plotUnorderedMatrix: boolean`
 
   Export the unordered matrix as pdf. Default value is `False`.

- `saveObjectiveFunctionValues: boolean`
 
   Export the values of the objective functions as csv. Default value is `False`.


### Notes
If there is a clustering passed as an argument through `rowClustersFile` and `columnClustersFile`, then `clusteringAlgorithm` and `nbClusters` will not be used.