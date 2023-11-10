# Visualizing Overlapping Biclusterings and Boolean Matrix Factorizations

README updated on the 27/06/2023

This archive contains the code to visualize biclusters from the paper (to be published)
Thibault Marette, Pauli Miettinen, Stefan Neumann:
Visualizing Overlapping Biclusterings and Boolean Matrix Factorizations. (to be published)


If you want to use basso clustering algorithm, extract `includes/basso-0.5.tar.gz`. Once the archive is extracted, follow the installation instructions at `includes/basso-0.5/README`.


For Python packages requirements, see `requirements.txt`.

## Execute the code

The code is available in the folder `biclusteringVisualization/`
Code examples are available in `demo.py`, `demoAdvanced.py` and `demoAdvanced2.py`

## Documentation

```python
visualizeClustering.visualizeClustering(
    inputFile,
    outputFile,
    clusteringAlgorithm="PCV",
    nbClusters=5,
    rowClusters="",
    columnClusters="",
    inputIsList=False,
    orderingMethod="TSPheuristic", 
    heuristicTime=20,
    plotUnorderedMatrix=False,
    showObjectiveFunctionValues=False
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

- `rowClusters: string`

   location of the file for the clustering of the rows. Optionnal (only use if you do not want to use clustering algorithm)

- `columnClusters: string | list`

   location of the file for the clustering of the columns. Optionnal (only use if you do not want to use clustering algorithm)

- `orderingMethod: string | list`
 
    Algorithm to use for the ordering of the matrix in the visualization. Default value is TSPHeuristic. ADVISER, as well as the greedy algorithms used for the experiments in the paper are also available.

- `heuristicTime: int`
 
    Execution time (in seconds) for TSPHeuristic. Default value is 20 seconds.
- `plotUnorderedMatrix: boolean`
 
   Export the unordered matrix as pdf. Default value is `False`.

- `saveObjectiveFunctionValues: boolean`
 
   Export the values of the objective functions as csv. Default value is `False`.

- `printPermutations: boolean`

   Print in the standard output the row and columns permutations. Default value is `False`.

- `customColorSchemes: boolean`

   Use the custom colour scheme. Default value is `True`.



### Notes
If there is a clustering passed as an argument through `rowClusters` and `columnClusters`, then `clusteringAlgorithm` and `nbClusters` will not be used.

Hence, there are two ways to use the code:
- using `clusteringAlgorithm` and `nbClusters` to visualize the clustering obtained through the given clustering algorithm (see `demo.py`).
- using `rowClusters` and `columnClusters` to visualize the clustering given as an input files (see `demoAdvanced.py`).
- using `rowClusters`, `columnCluster` and `inputIsList` to visualize the clustering given as python lists (see `demoAdvanced2.py`).


Contact address: Thibault Marette: marette@kth.se
