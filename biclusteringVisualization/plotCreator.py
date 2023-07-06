import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import numpy as np
import time
from .ComputeClusters import cluster
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

"""
This file contains the implementation of the TSPHeuristic algorithm from our paper

Thibault Marette, Pauli Miettinen, Stefan Neumann:
Visualizing Overlapping Biclusterings and Boolean Matrix Factorizations. (to be published)


It also contains an implementation by Thibault Marette and Stefan Neumann of the
ADVISER algorithm from the following paper:

Alessandro Colantonio, Roberto Di Pietro, Alberto Ocello, Nino Vincenzo Verde:
Visual Role Mining: A Picture Is Worth a Thousand Roles. IEEE Trans. Knowl. Data Eng. 24(6): 1120-1133 (2012)

Please do not share this code without our expressed permission.
"""


"""
	
	Parameters:
		A:
			the original matrix as numpy matrix
		rowClusters:
			the row clusters
		columnClusters:
			the column clusters
		method:
			the method that shall be used for clustering
			should be either TSPheuristic or ADVISER
		plotUnorderedMatrix:
			if True, plots the original matrix A (without reordering)
		transpose:
			if True, plots the transpose of the matrices
		outputfile:
			if not None, writes the plot to the file given by the string in
			outputfile (e.g., outputfile='~/Desktop/myGreatPlot.pdf')
"""

ROW = "row"
COL = "column"


def readClusterFromFile(inputFile):
    f = open(inputFile, "r")
    clusters = []
    for line in f.readlines():
        clusters.append([int(l) for l in line.split(",")])
    return clusters


def visualizeClustering(
    inputFile,
    outputFile,
    clusteringAlgorithm="PCV",
    nbClusters=5,
    orderingMethod="TSPheuristic",
    rowClusters="",
    columnClusters="",
    inputIsList=False,
    heuristicTime=20,
    plotUnorderedMatrix=False,
    showObjectiveFunctionValues=False,
):
    A = np.loadtxt(inputFile)
    if inputIsList:
        pass
    elif rowClusters != "" and columnClusters != "":
        rowClusters = readClusterFromFile(rowClusters)
        columnClusters = readClusterFromFile(columnClusters)
    
    else:
        [rowClusters, columnClusters] = cluster(A, nbClusters, algo=clusteringAlgorithm)
    " visualize the matrix "
    plotClusteredMatrix(
        A,
        rowClusters,
        columnClusters,
        method=orderingMethod,
        heuristicTime=heuristicTime,
        transpose=False,
        outputfile=outputFile,
        plotUnorderedMatrix=plotUnorderedMatrix,
        showObjectiveFunctionValues=showObjectiveFunctionValues,
    )


def plotClusteredMatrix(
    A,
    rowClusters,
    columnClusters,
    method="TSPheuristic",
    heuristicTime=5,
    plotUnorderedMatrix=False,
    transpose=False,
    outputfile=None,
    showObjectiveFunctionValues=False,
):
    if transpose:
        plotClusteredMatrix(
            A.transpose(),
            columnClusters,
            rowClusters,
            method=method,
            plotUnorderedMatrix=plotUnorderedMatrix,
            transpose=False,
            outputfile=outputfile,
        )
        return

    " this is where the re-ordering and post-processing happens "
    blockBlockMatrix = NewBlockBlockMatrix(rowClusters, columnClusters, A)
    timeStart = time.time()
    blockBlockMatrix.order(method=method, heuristicTime=heuristicTime)

    runningTime = time.time() - timeStart
    clusterAreaScore = blockBlockMatrix.clusterAreaScore()
    proximityScore = blockBlockMatrix.totalProximityScore()
    demeritScore = blockBlockMatrix.globalDemerit(ROW) + blockBlockMatrix.globalDemerit(
        COL
    )
    visualisationCost = blockBlockMatrix.globalVisualizationCost()
    uninterScore = blockBlockMatrix.blockSliceAreaScore()

    " the rest from here is for creating the plots "

    if plotUnorderedMatrix:
        fig1, ax1 = plt.subplots()
        fig1.tight_layout()
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            labelbottom=False,
            right=False,
            left=False,
            labelleft=False,
        )
        ax1.matshow(blockBlockMatrix.originalMatrix)
        plt.axis("off")

        fig1.savefig(
            outputfile.replace(".pdf", "") + "_unordered.pdf",
            bbox_inches="tight",
            pad_inches=0,
        )
    fig2, ax2 = plt.subplots()
    fig2.tight_layout()
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        labelbottom=False,
        right=False,
        left=False,
        labelleft=False,
    )
    blockBlockMatrix.plotOriginalMatrix(
        axs=ax2,
        permutationRow=blockBlockMatrix.getFinalPermutation(ROW),
        permutationCol=blockBlockMatrix.getFinalPermutation(COL),
    )
    plt.axis("off")
    fig2.savefig(
        outputfile,
        bbox_inches="tight",
        pad_inches=0,
    )

    plt.close("all")
    if showObjectiveFunctionValues:
        f = open(outputfile.replace(".pdf", "") + "_metrics.csv", "w")
        f.write(
            "runningTime, clusterAreaScore, uninterScore, proximityScore, demeritScore, visualizationCost\n"
        )
        f.write(
            f"{runningTime}, {clusterAreaScore}, {uninterScore}, {proximityScore}, {demeritScore}, {visualisationCost}"
        )


def createBlocks(clusters):
    indexToClusters = {}
    for i in range(len(clusters)):
        cluster = clusters[i]
        for index in cluster:
            if index not in indexToClusters:
                indexToClusters[index] = []
            indexToClusters[index].append(i)
    blocks = {}
    allIndices = indexToClusters.keys()

    for index in allIndices:
        blockName = str(indexToClusters[index])
        if blockName not in blocks:
            blocks[blockName] = []
        blocks[blockName].append(index)
    blocks = list(blocks.values())
    return blocks


def createOrphanBlocks(fosterFamilies):
    indexToClusters = {}
    indexToBlueprint = {}
    for i in range(len(fosterFamilies)):
        orphans = fosterFamilies[i]["orphans"]
        for index in list(orphans):
            if index not in indexToClusters:
                indexToClusters[index] = []
            if index not in indexToBlueprint:
                indexToBlueprint[index] = None
            indexToClusters[index].append(fosterFamilies[i]["cluster"])
            indexToBlueprint[index] = fosterFamilies[i]["blueprint"]
    blocks = {}
    allIndices = indexToClusters.keys()
    for index in allIndices:
        blockName = str(indexToClusters[index])
        if blockName not in blocks:
            blocks[blockName] = {
                "block": blockName,
                "orphans": [],
                "blueprint": indexToBlueprint[index],
            }
        blocks[blockName]["orphans"].append(index)
    orphanBlocks = []
    for block in blocks.values():
        orphanBlocks.append(block)
    return orphanBlocks


class NewBlockBlockMatrix:
    def __init__(self, rowClusters, columnClusters, originalMatrix=[]):
        self.originalMatrix = np.array(originalMatrix)
        self.rows, self.cols = originalMatrix.shape
        self.rowClusters = []
        self.columnClusters = []
        for i in range(len(rowClusters)):
            # We remove empty clusters (if the row or column side of the cluster is empty)
            if (rowClusters[i] != []) and (columnClusters[i] != []):
                self.rowClusters.append(rowClusters[i])
                self.columnClusters.append(columnClusters[i])
        self.nbClusters = len(self.rowClusters)
        self.clustersSize = [
            len(self.rowClusters[i]) * len(self.columnClusters[i])
            for i in range(self.nbClusters)
        ]
        # blocks: block view of the matrix, the rows and columns are grouped together in a block when they belong to the same clusters
        self.rowBlocks = createBlocks(self.rowClusters)
        self.columnBlocks = createBlocks(self.columnClusters)
        self.m = len(self.rowBlocks)
        self.n = len(self.columnBlocks)
        self.rowBlocksSize = [len(block) for block in self.rowBlocks]
        self.columnBlocksSize = [len(block) for block in self.columnBlocks]

        # block clusters: clusters, but with block index as reference instead of rows/columns
        self.rowClustersBlock = self.createBlockedClusters(
            self.rowClusters, self.rowBlocks
        )
        self.columnClustersBlock = self.createBlockedClusters(
            self.columnClusters, self.columnBlocks
        )
        # clusters block: blocks, but use clusters index as reference instead of rows/columns
        self.rowBlocksCluster = self.createClusteredBlocks(
            self.rowBlocks, self.rowClustersBlock
        )
        self.columnBlocksCluster = self.createClusteredBlocks(
            self.columnBlocks, self.columnClustersBlock
        )
        # Collapsed block to speedup computation of lambda score
        self.rowCollapsedBlocks = self.createCollapsedBlocks(ROW)
        self.columnCollapsedBlocks = self.createCollapsedBlocks(COL)
        # blocks permutation: the permutation used in the plot
        self.idealizedRows = self.createIdealizedItems(self.rowClusters, self.rows)
        self.idealizedColumns = self.createIdealizedItems(
            self.columnClusters, self.cols
        )

        self.rowBlocksPermutation = [i for i in range(self.m)]
        self.columnBlocksPermutation = [i for i in range(self.n)]

        # ORPHAN ORDERING
        self.rowClustersBlueprint = self.getClustersBlueprint(ROW)
        self.columnClustersBlueprint = self.getClustersBlueprint(COL)
        self.rowClustersSparcity = self.getClustersSparcity(ROW)
        self.columnClustersSparcity = self.getClustersSparcity(COL)
        self.rowOrphans = self.getOrphans(ROW)
        self.columnOrphans = self.getOrphans(COL)

        self.rowFosterFamilies = self.getFosterFamilies(ROW)
        self.columnFosterFamilies = self.getFosterFamilies(COL)
        self.rowOrphanBlocks = createOrphanBlocks(self.rowFosterFamilies)
        self.columnOrphanBlocks = createOrphanBlocks(self.columnFosterFamilies)

        # block adjacency matrix
        self.C = (
            self.createBlockBlockAdj()
        )  # C[i,j] = 1 means that block i and j belong to the same cluster, and thus all indices in these block belong to the same cluster.

    def getClustersSparcity(self, mode):
        if mode == ROW:
            clusters = self.rowClusters
            items = [
                [np.array(self.originalMatrix[row]) for row in cluster]
                for cluster in clusters
            ]
            idealizedImages = [
                [1 if i in cluster else 0 for i in range(self.cols)]
                for cluster in self.columnClusters
            ]
            size = self.cols

        if mode == COL:
            clusters = self.columnClusters
            items = [
                [np.array(self.originalMatrix[:, column]) for column in cluster]
                for cluster in clusters
            ]
            idealizedImages = [
                [1 if i in cluster else 0 for i in range(self.rows)]
                for cluster in self.rowClusters
            ]
            size = self.rows
        sparcities = []
        for i in range(len(items)):
            item = items[i]
            idealizedImage = idealizedImages[i]
            sparcity = 0
            for cluster in item:
                sparcity += np.sum(cluster + idealizedImage == 2)
            sparcities.append(sparcity / (len(item) * size))
        return sparcities

    def createIdealizedItems(self, xClusters, itemsLength):
        idealizedItems = [[] for _ in range(itemsLength)]
        for i in range(len(xClusters)):
            for item in xClusters[i]:
                idealizedItems[item].append(i)
        return idealizedItems

    def createBlockBlockAdj(self):
        C = np.zeros([self.m, self.n])
        for i in range(self.m):
            for j in range(self.nbClusters):
                rowBlock = set(self.rowBlocks[i])
                rowCluster = set(self.rowClusters[j])
                if rowBlock.issubset(rowCluster):
                    for k in range(self.n):
                        columnBlock = set(self.columnBlocks[k])
                        columnCluster = set(self.columnClusters[j])
                        if columnBlock.issubset(columnCluster):
                            C[i, k] = len(rowBlock) * len(columnBlock)
        return C

    def createCollapsedBlocks(self, mode):
        if mode == ROW:
            blockLength = self.cols
            xBlocks = self.rowBlocks
        if mode == COL:
            blockLength = self.rows
            xBlocks = self.columnBlocks
        xCollapsedBlocks = []
        for block in xBlocks:
            collapsedBlock = np.zeros(blockLength)
            if mode == ROW:
                collapsedBlock = sum(self.originalMatrix[block])
            if mode == COL:
                collapsedBlock = np.transpose(
                    sum(np.transpose(self.originalMatrix[:, block]))
                )
            collapsedBlock = np.where(collapsedBlock / len(block) > 0.4, 1, 0)
            xCollapsedBlocks.append(collapsedBlock)
        return xCollapsedBlocks

    def createBlockedClusters(self, clusters, blocks):
        blockedClusters = []
        for cluster in clusters:
            blockedCluster = []
            for block in blocks:
                if set(block).issubset(set(cluster)):
                    blockedCluster.append(blocks.index(block))
            blockedClusters.append(blockedCluster)
        return blockedClusters

    def createClusteredBlocks(self, xBlocks, xClustersBlock):
        xBlocksCluster = [set() for _ in xBlocks]
        for i in range(len(xClustersBlock)):
            blocks = xClustersBlock[i]  # blocks that belong to cluster i
            for block in list(blocks):
                xBlocksCluster[block].add(i)

        return xBlocksCluster

    ################################################################################################
    ########## Translating blockPermutations into a permutation for the original matrix ############
    ################################################################################################

    def getRowPermutation(self):  # , sn=False):
        finalPermutation = []
        for blockIndex in self.rowBlocksPermutation:
            newBlock = self.rowBlocks[blockIndex]
            finalPermutation += newBlock
        orphans = []
        for node in range(self.rows):
            if node not in finalPermutation:
                orphans.append(node)
        return finalPermutation + orphans

    def getOrphans(self, mode):
        if mode == ROW:
            items = self.rows
            blocks = self.rowBlocks
        if mode == COL:
            items = self.cols
            blocks = self.columnBlocks
        allBlocks = []
        for block in blocks:
            allBlocks += block
        orphans = []
        for node in range(items):
            if node not in allBlocks:
                orphans.append(node)
        return orphans

    def getClusterPermutation(self, mode):
        finalPermutation = []
        if mode == ROW:
            xBlocksPermutation = self.rowBlocksPermutation
            xBlocks = self.rowBlocks
        if mode == COL:
            xBlocksPermutation = self.columnBlocksPermutation
            xBlocks = self.columnBlocks

        for blockIndex in xBlocksPermutation:
            newBlock = xBlocks[blockIndex]
            finalPermutation += newBlock
        return finalPermutation

    def getOrphanPermutation(self, mode):
        if mode == ROW:
            adoptedOrphans = self.rowOrphanBlocks
            blockPermutation = self.rowBlocksPermutation
            blocksCluster = self.rowBlocksCluster

        if mode == COL:
            adoptedOrphans = self.columnOrphanBlocks
            blockPermutation = self.columnBlocksPermutation
            blocksCluster = self.columnBlocksCluster

        adoptedOrphanPermutation = []
        for b in blockPermutation:
            bc = str(sorted(blocksCluster[b]))
            for orphan in adoptedOrphans:
                if orphan["block"] == bc:
                    adoptedOrphanPermutation.append(orphan)
        for orphans in adoptedOrphans:
            if orphans not in adoptedOrphanPermutation:
                adoptedOrphanPermutation.append(orphans)
        return adoptedOrphanPermutation

    def getFinalPermutation(self, mode):
        if mode == ROW:
            size = self.rows
        if mode == COL:
            size = self.cols
        clusterPermutation = self.getClusterPermutation(mode)
        orphanPermutation = self.getOrphanPermutation(mode)
        leftOrphanPermutation = orphanPermutation[::-1][len(orphanPermutation) // 2 :]
        rightOrphanPermutation = orphanPermutation[::-1][: len(orphanPermutation) // 2]
        orphansLeft = 0
        finalPermutation = []
        for orphan in leftOrphanPermutation:
            finalPermutation += orphan["orphans"]
            orphansLeft += len(orphan["orphans"])
        finalPermutation += clusterPermutation
        for orphan in rightOrphanPermutation:
            finalPermutation += orphan["orphans"]
        noise = []
        for node in range(size):
            if node not in finalPermutation:
                noise.append(node)
        if mode == ROW:
            self.startOfRowCluster = len(noise[: len(noise) // 2]) + orphansLeft
        if mode == COL:
            self.startOfColumnCluster = len(noise[: len(noise) // 2]) + orphansLeft
        self.endOfCluster = len(noise[len(noise) // 2 :])
        return noise[: len(noise) // 2] + finalPermutation + noise[len(noise) // 2 :]

    def getColumnPermutation(self):
        finalPermutation = []
        for blockIndex in self.columnBlocksPermutation:
            newBlock = self.columnBlocks[blockIndex]
            finalPermutation += newBlock
        orphans = []
        for node in range(self.cols):
            if node not in finalPermutation:
                orphans.append(node)
        return finalPermutation + orphans

    def itemSimilarityScore(self, i, j, mode):
        if mode == ROW:
            itemi = self.originalMatrix[i]
            itemj = self.originalMatrix[j]
        if mode == COL:
            itemi = self.originalMatrix[:, i]
            itemj = self.originalMatrix[:, j]
        fusionCollapsedBlocks = itemi + itemj
        interScore = np.count_nonzero(fusionCollapsedBlocks == 2) + np.count_nonzero(
            fusionCollapsedBlocks == 0
        )
        unionScore = np.count_nonzero(fusionCollapsedBlocks == 1) + interScore
        if unionScore == 0:
            return 0
        return interScore / unionScore

    def getClustersBlueprint(self, mode):
        if mode == ROW:
            xCluster = self.rowClusters
            items = self.rows
        if mode == COL:
            xCluster = self.columnClusters
            items = self.cols
        bases = []
        for i in range(len(xCluster)):
            clusters = xCluster[i]
            blueprint = np.array([i in clusters for i in range(items)])
            bases.append({"cluster": i, "blueprint": blueprint})
        return bases

    def getFosterFamilies(self, mode):
        if mode == ROW:
            orphansIndexes = self.rowOrphans
            orphansPattern = [np.array(self.originalMatrix[o]) for o in orphansIndexes]
            # print(self.columnClustersBlueprint)
            blueprint = self.columnClustersBlueprint
            clustersSparcity = self.rowClustersSparcity
        if mode == COL:
            orphansIndexes = self.columnOrphans
            orphansPattern = [
                np.array(self.originalMatrix[:, o]) for o in orphansIndexes
            ]
            blueprint = self.rowClustersBlueprint
            clustersSparcity = self.columnClustersSparcity
        if orphansIndexes == []:
            return []
        fosterFamilies = []
        for j in range(len(blueprint)):
            fosterFamily = {
                "blueprint": blueprint[j]["blueprint"],
                "cluster": blueprint[j]["cluster"],
                "orphans": set(),
            }
            b = blueprint[j]["blueprint"]
            for i in range(len(orphansPattern)):
                o = orphansPattern[i]
                similarity = sum((b + o == 2) / len(o))
                if similarity > clustersSparcity[j] / 2:
                    fosterFamily["orphans"].add(orphansIndexes[i])
            fosterFamilies.append(fosterFamily)
        return fosterFamilies

    #########################################################################################
    ############## Dispatching which ordering technique to what function. ###################
    #########################################################################################

    def metaTSP(self, heuristicTime=20):
        sigmaRows = self.TSPheuristic(ROW, heuristicTime // 2)
        sigmaCols = self.TSPheuristic(COL, heuristicTime // 2)
        bestClusterArea = 0
        bestSliceBlockArea = 0
        bestSigmaRows = sigmaRows
        bestSigmaCols = sigmaCols
        self.columnBlocksPermutation = sigmaCols
        for i in range(len(sigmaRows)):
            currentSigmaRows = sigmaRows[i:] + sigmaRows[:i]
            self.rowBlocksPermutation = currentSigmaRows
            blockSliceAreaScore = self.blockSliceAreaScore()
            clusterAreaScore = self.clusterAreaScore()
            if clusterAreaScore > bestClusterArea:
                bestSigmaRows = currentSigmaRows
                bestClusterArea = clusterAreaScore
                bestSliceBlockArea = bestSliceBlockArea
            if (
                clusterAreaScore == bestClusterArea
                and blockSliceAreaScore > bestSliceBlockArea
            ):
                bestSigmaRows = currentSigmaRows
                bestClusterArea = clusterAreaScore
                bestSliceBlockArea = bestSliceBlockArea
        for j in range(len(sigmaCols)):
            currentSigmaCols = sigmaCols[j:] + sigmaCols[:j]
            self.columnBlocksPermutation = currentSigmaCols
            blockSliceAreaScore = self.blockSliceAreaScore()
            clusterAreaScore = self.clusterAreaScore()
            if clusterAreaScore > bestClusterArea:
                bestSigmaCols = currentSigmaCols
                bestClusterArea = clusterAreaScore
                bestSliceBlockArea = bestSliceBlockArea
            if (
                clusterAreaScore == bestClusterArea
                and blockSliceAreaScore > bestSliceBlockArea
            ):
                bestSigmaCols = currentSigmaCols
                bestClusterArea = clusterAreaScore
                bestSliceBlockArea = bestSliceBlockArea
        self.rowBlocksPermutation = bestSigmaRows
        self.columnBlocksPermutation = bestSigmaCols

    def order(self, method="TSPheuristic", heuristicTime=20):
        if method == "TSPheuristic":
            self.metaTSP(heuristicTime)
        elif method == "ADVISER":
            self.adviser()
        elif method == "greedyProximityScore":
            self.greedyProximityScore()
        elif method == "greedyClusterArea":
            self.greedyClusterArea()
        elif method == "greedyUninter":
            self.greedyBlockSliceArea(ROW)
            self.greedyBlockSliceArea(COL)
        elif method == "greedyDemerit":
            self.greedyDemerit2D()
        elif "random" in method:
            np.random.shuffle(self.rowBlocksPermutation)
            np.random.shuffle(self.columnBlocksPermutation)
        else:
            print(f"Unknown ordering method {method}. No changes to ordering.")

    ### greedy algorithms
    def greedyProximityScore(self):
        print("Running greedy greedyProximityScore")
        rowToAddOrder = self.sortForGreedy(ROW)

        columnToAddOrder = self.sortForGreedy(COL)

        sigmaRow = [rowToAddOrder[0]]
        sigmaColumn = [columnToAddOrder[0]]

        # we sort the blocks by the criterion in the paper.
        while len(sigmaRow) != len(self.rowBlocks) or len(sigmaColumn) != len(
            self.columnBlocks
        ):
            if len(sigmaRow) != len(self.rowBlocks):
                sigmaRow = self.addProximityScore(
                    sigmaRow, sigmaColumn, rowToAddOrder[len(sigmaRow)], ROW
                )
            if len(sigmaColumn) != len(self.columnBlocks):
                sigmaColumn = self.addProximityScore(
                    sigmaRow, sigmaColumn, columnToAddOrder[len(sigmaColumn)], COL
                )
        self.rowBlocksPermutation = sigmaRow
        self.columnBlocksPermutation = sigmaColumn

    def addProximityScore(self, sigmaRow, sigmaColumn, i, mode):
        if mode == ROW:
            sigma = sigmaRow
        if mode == COL:
            sigma = sigmaColumn
        bestScore = 1e100
        place = 0
        for j in range(len(sigma)):
            if mode == COL:
                score = self.partialproximityScore(
                    sigmaRow, sigma[j:] + [i] + sigma[:j]
                )
            if mode == ROW:
                score = self.partialproximityScore(
                    sigma[j:] + [i] + sigma[:j], sigmaColumn
                )
            if score < bestScore:
                bestScore = score
                place = j
        sigma = sigma[:place] + [i] + sigma[place:]
        if mode == ROW:
            self.rowBlocksPermutation = sigma
        if mode == COL:
            self.columnBlocksPermutation = sigma
        return sigma

    def greedyClusterArea(self):
        print("Running greedy clusterArea")
        rowToAddOrder = self.sortForGreedy(ROW)

        columnToAddOrder = self.sortForGreedy(COL)

        sigmaRow = [rowToAddOrder[0]]
        sigmaColumn = [columnToAddOrder[0]]

        # we sort the blocks by the criterion in the paper.
        while len(sigmaRow) != len(self.rowBlocks) or len(sigmaColumn) != len(
            self.columnBlocks
        ):
            if len(sigmaRow) != len(self.rowBlocks):
                sigmaRow = self.addGreedyClusterArea(
                    sigmaRow, sigmaColumn, rowToAddOrder[len(sigmaRow)], ROW
                )
            if len(sigmaColumn) != len(self.columnBlocks):
                sigmaColumn = self.addGreedyClusterArea(
                    sigmaRow, sigmaColumn, columnToAddOrder[len(sigmaColumn)], COL
                )
        self.rowBlocksPermutation = sigmaRow
        self.columnBlocksPermutation = sigmaColumn

    def addGreedyClusterArea(self, sigmaRow, sigmaColumn, i, mode):
        if mode == ROW:
            sigma = sigmaRow
        if mode == COL:
            sigma = sigmaColumn
        bestScore = 0
        place = 0
        for j in range(len(sigma)):
            if mode == COL:
                score = self.partialClusterArea(sigmaRow, sigma[j:] + [i] + sigma[:j])
            if mode == ROW:
                score = self.partialClusterArea(
                    sigma[j:] + [i] + sigma[:j], sigmaColumn
                )
            if score > bestScore:
                bestScore = score
                place = j
        sigma = sigma[:place] + [i] + sigma[place:]
        if mode == ROW:
            self.rowBlocksPermutation = sigma
        if mode == COL:
            self.columnBlocksPermutation = sigma
        return sigma

    def sortForGreedy(self, mode):
        if mode == ROW:
            tiles = self.rowBlocks
            xClustersTiles = self.rowClustersBlock
            clusters = self.rowClusters
        if mode == COL:
            tiles = self.columnBlocks
            xClustersTiles = self.columnClustersBlock
            clusters = self.columnClusters
        IA = []
        for i in range(len(tiles)):
            block = tiles[i]
            item = block[0]
            for r in range(len(xClustersTiles)):
                c = clusters[r]  # c is cluster number r
                if item in c:
                    IA.append((i, block, r))
        return sorted(
            list(range(len(tiles))), key=lambda x: -self.importanceScore(x, IA)
        )

    def greedyBlockSliceArea(self, mode):
        print("Running greedy blockSliceArea")
        if mode == ROW:
            tiles = self.rowBlocks
        if mode == COL:
            tiles = self.columnBlocks
        nextTiles = self.sortForGreedy(mode)

        # here items is itemsBar already, since the items have been grouped by block.
        sigma = [nextTiles[0]]
        # we sort the blocks by the criterion in the paper.
        while len(sigma) != len(tiles):
            maxScore = 0
            i = nextTiles[len(sigma)]
            bestPlace = 0
            for j in range(len(sigma)):
                score = self.partialblockSliceArea(
                    0, len(sigma) + 1, sigma[j:] + [i] + sigma[:j], mode
                )
                if score > maxScore:
                    maxScore = score
                    bestPlace = j
            sigma = sigma[:bestPlace] + [i] + sigma[bestPlace:]
        if mode == ROW:
            self.rowBlocksPermutation = sigma
        if mode == COL:
            self.columnBlocksPermutation = sigma

    def blockSliceAreaScore(self):
        return self.blockSliceAreaScore1D(COL) + self.blockSliceAreaScore1D(ROW)

    def partialClusterArea(self, partialRowPermutation, partialColumnPermutation):
        nbRows = len(partialRowPermutation)
        nbColumns = len(partialColumnPermutation)
        rowBlocksCluster = [
            set(self.rowBlocksCluster[partialRowPermutation[i]]) for i in range(nbRows)
        ]
        rowBlockSize = [
            len(self.rowBlocks[partialRowPermutation[i]]) for i in range(nbRows)
        ]
        columnBlocksCluster = [
            set(self.columnBlocksCluster[partialColumnPermutation[i]])
            for i in range(nbColumns)
        ]
        columnBlockSize = [
            len(self.columnBlocks[partialColumnPermutation[i]])
            for i in range(nbColumns)
        ]

        score = 0
        for c in range(self.nbClusters):
            compressedRow = [0]
            for i in range(len(rowBlocksCluster)):
                rowBlockCluster = rowBlocksCluster[i]
                if c in rowBlockCluster:
                    compressedRow[-1] += rowBlockSize[i]
                else:
                    compressedRow.append(0)

            compressedColumn = [0]
            for i in range(len(columnBlocksCluster)):
                columnBlockCluster = columnBlocksCluster[i]
                if c in columnBlockCluster:
                    compressedColumn[-1] += columnBlockSize[i]
                else:
                    compressedColumn.append(0)
            for r in compressedRow:
                for c in compressedColumn:
                    score += (r * c) ** 2
        return score

    def partialblockSliceArea(self, firstTile, lastTile, partialTilePermutation, mode):
        if mode == ROW:
            line = self.columnBlocks
            clusters = self.columnBlocksCluster
            tLines = self.rowBlocks
            tClusters = self.rowBlocksCluster
        if mode == COL:
            line = self.rowBlocks
            clusters = self.rowBlocksCluster
            tLines = self.columnBlocks
            tClusters = self.columnBlocksCluster

        interval = range(firstTile, lastTile)
        tClusters = [set(tClusters[partialTilePermutation[i]]) for i in interval]
        tSize = [len(tLines[partialTilePermutation[i]]) for i in interval]
        score = 0
        for i in range(len(line)):
            ci = set(clusters[i])
            if ci.intersection(set(tClusters[firstTile])) != set():
                compressedLine = [tSize[firstTile]]
            else:
                compressedLine = [0]
            for j in interval[1:]:
                cj = set(tClusters[j])
                if cj.intersection(ci) == set():
                    compressedLine.append(tSize[j])
                else:
                    compressedLine[-1] += tSize[j]
            for consecutiveElement in compressedLine:
                score += (len(line[i]) * consecutiveElement) ** 2
        return score

    def partialRoleSize(self, c, partialSigmaRow, partialSigmaColumn):
        blocksInRowCluster = self.rowClustersBlock[c]  # indices of blocks in cluster c
        blocksInColumnCluster = self.columnClustersBlock[c]
        blocksInRowCluster = sorted(
            [
                partialSigmaRow.index(i)
                for i in blocksInRowCluster
                if i in partialSigmaRow
            ]
        )
        blocksInColumnCluster = sorted(
            [
                partialSigmaColumn.index(i)
                for i in blocksInColumnCluster
                if i in partialSigmaColumn
            ]
        )
        if blocksInRowCluster == [] or blocksInColumnCluster == []:
            return 0
        firstIndex = blocksInRowCluster[0]
        lastIndex = blocksInRowCluster[-1]
        sizeRow = 0
        for index in range(firstIndex, lastIndex + 1):
            currentBlock = partialSigmaRow[index]
            sizeRow += len(self.rowBlocks[currentBlock])
        firstIndex = blocksInRowCluster[0]
        lastIndex = blocksInRowCluster[-1]
        sizeColumn = 0
        for index in range(firstIndex, lastIndex + 1):
            currentBlock = partialSigmaRow[index]
            sizeColumn += len(self.rowBlocks[currentBlock])
        return sizeRow * sizeColumn

    def partialproximityScore(self, sigmaRows, sigmaColumn):
        score = 0
        for c in range(self.nbClusters):
            score += self.partialRoleSize(c, sigmaRows, sigmaColumn)
        return score

    def constantBlockSliceAreaScore1D(self, mode):
        if mode == COL:
            line = self.columnBlocks
            clusters = self.columnBlocksCluster
            tLine = self.rowBlocks
            tClusters = self.rowBlocksCluster
        if mode == ROW:
            line = self.rowBlocks
            clusters = self.rowBlocksCluster
            tLine = self.columnBlocks
            tClusters = self.columnBlocksCluster

        tSize = [len(tLine[i]) for i in range(len(tLine))]
        score = 0
        for i in range(len(line)):
            ci = set(clusters[i])
            for j in range(0, len(tLine)):
                cj = set(tClusters[j])
                if cj.intersection(ci) == set():
                    score += tSize[j] * len(line[i]) ** 2
        return score

    def blockSliceAreaScore1D(self, mode):
        if mode == COL:
            line = self.columnBlocks
            clusters = self.columnBlocksCluster
            tLine = self.rowBlocks
            tClusters = self.rowBlocksCluster
            tPermutation = self.rowBlocksPermutation
        if mode == ROW:
            line = self.rowBlocks
            clusters = self.rowBlocksCluster
            tLine = self.columnBlocks
            tClusters = self.columnBlocksCluster
            tPermutation = self.columnBlocksPermutation

        tLine = [tLine[tPermutation[i]] for i in range(len(tLine))]
        tClusters = [set(tClusters[tPermutation[i]]) for i in range(len(tLine))]
        tSize = [len(tLine[i]) for i in range(len(tLine))]
        score = 0
        for i in range(len(line)):
            ci = set(clusters[i])
            if ci.intersection(set(tClusters[0])) != set():
                compressedLine = [tSize[0]]
            else:
                compressedLine = [0]
            for j in range(1, len(tLine)):
                cj = set(tClusters[j])
                if cj.intersection(ci) == set():
                    compressedLine.append(tSize[j])
                else:
                    compressedLine[-1] += tSize[j]
            for consecutiveElement in compressedLine:
                score += (len(line[i]) * consecutiveElement) ** 2
        return score

    def demerit(self, s1, s2):
        if s1 == set() or s2 == set():
            return 2 * len(s1.union(s2))
        if s1.intersection(s2) == set():
            return len(s1.union(s2))  # the common cluster
        return len(s1.union(s2)) - len(s1.intersection(s2))

    def globalDemerit(self, mode):
        if mode == ROW:
            line = self.rowBlocks
            linePermutation = self.rowBlocksPermutation
        if mode == COL:
            line = self.columnBlocks
            linePermutation = self.columnBlocksPermutation
        return sum(
            [
                self.totalDemerit(linePermutation[i], linePermutation[i + 1], mode)
                for i in range(len(line) - 1)
            ]
        )

    def totalDemerit(self, i, j, mode):
        if mode == ROW:
            # we order the rows, so we look at the columns
            line = self.columnBlocks
            clusters = self.columnBlocksCluster
            tiles = self.rowBlocks
            tilesCluster = self.rowBlocksCluster
        if mode == COL:
            line = self.rowBlocks
            clusters = self.rowBlocksCluster
            tiles = self.columnBlocks
            tilesCluster = self.columnBlocksCluster

        ci = set(tilesCluster[i])
        cj = set(tilesCluster[j])
        totalDemerit = 0
        for k in range(len(line)):
            lc = set(clusters[k])
            d = self.demerit(lc.intersection(ci), lc.intersection(cj))
            totalDemerit += d * len(line[k])  # * len(tiles[j])
        return totalDemerit

    def greedyDemerit2D(self):
        print("running greedu demerit")
        self.greedyDemerit(ROW)
        self.greedyDemerit(COL)

    def TSPheuristic(self, mode, heuristicTime):
        print("running TSP heuristic")
        if mode == ROW:
            # we order the rows, so we look at the columns
            tiles = self.rowBlocks
        if mode == COL:
            tiles = self.columnBlocks
        distanceMatrix = np.array(
            [
                np.array([self.totalDemerit(i, j, mode) for i in range(len(tiles))])
                for j in range(len(tiles))
            ]
        )
        manager = pywrapcp.RoutingIndexManager(len(tiles), 1, 0)

        # Create Routing Model.
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            """Returns the distance between the two nodes."""
            # Convert from routing variable Index to distance matrix NodeIndex.
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return distanceMatrix[from_node][to_node]

        transit_callback_index = routing.RegisterTransitCallback(distance_callback)

        # Define cost of each arc.
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        # Setting first solution heuristic.
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (
            routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        )
        search_parameters.local_search_metaheuristic = (
            routing_enums_pb2.LocalSearchMetaheuristic.SIMULATED_ANNEALING
        )
        search_parameters.time_limit.seconds = max(1, heuristicTime)

        # Solve the problem.
        solution = routing.SolveWithParameters(search_parameters)
        index = routing.Start(0)
        sigma = []
        while not routing.IsEnd(index):
            sigma.append(manager.IndexToNode(index))
            index = solution.Value(routing.NextVar(index))
        return sigma

    def greedyDemerit(self, mode):
        nextTiles = self.sortForGreedy(mode)
        permutation = [nextTiles[0]]
        nextTiles = nextTiles[1:]
        for nextTile in nextTiles:
            if self.totalDemerit(nextTile, permutation[0], mode) < self.totalDemerit(
                nextTile, permutation[-1], mode
            ):
                p = 0
                j = self.totalDemerit(nextTile, permutation[0], mode)
            else:
                p = len(permutation)
                j = self.totalDemerit(nextTile, permutation[-1], mode)
            for i in range(2, len(permutation)):
                jprec = self.totalDemerit(nextTile, permutation[i - 1], mode)
                jsucc = self.totalDemerit(nextTile, permutation[i], mode)
                jcurr = self.totalDemerit(permutation[i], permutation[i - 1], mode)
                if min(jprec, jsucc) < j and max(jprec, jsucc) <= jcurr:
                    p = i
                    j = min(jprec, jsucc)
            permutation.insert(p, nextTile)
        if mode == ROW:
            self.rowBlocksPermutation = permutation
        if mode == COL:
            self.columnBlocksPermutation = permutation

    #######################################################################
    ######################### ADVISER ALGORITHM ###########################
    #######################################################################

    def adviser(self):
        print("Running ADVISER...")
        # We keep the same notations as in the paper
        # row = users
        # column = perms
        # clusters = roles
        # UA = all role-user associations

        users = self.rowBlocks
        perms = self.columnBlocks

        UA = []
        for b in range(len(users)):
            block = users[b]
            user = block[0]  # all users in the block share the same rank
            for r in range(len(self.rowClustersBlock)):
                c = self.rowClusters[r]  # c is cluster number r
                if user in c:
                    UA.append((b, block, r))

        PA = []
        for b in range(len(perms)):
            block = perms[b]
            user = block[0]  # all users in the block share the same rank
            for r in range(len(self.columnClustersBlock)):
                c = self.columnClusters[r]  # c is cluster number r
                if user in c:
                    PA.append((b, block, r))
        userOrder = self.sortSet(range(len(self.rowBlocks)), UA, PA)
        permsOrder = self.sortSet(range(len(self.columnBlocks)), PA, UA)
        self.rowBlocksPermutation = userOrder
        self.columnBlocksPermutation = permsOrder

    def sortSet(self, items, IA, IAbar):
        # here items is itemsBar already, since the items have been grouped by block.
        sigma = []
        # we sort the blocks by the criterion in the paper.
        items = sorted(items, key=lambda x: -self.importanceScore(x, IA))
        for I in items:
            if len(sigma) < 2:
                sigma.append(I)
            else:
                if self.jacc(I, sigma[0], IA, IAbar) > self.jacc(
                    I, sigma[-1], IA, IAbar
                ):
                    p = 0
                    j = self.jacc(I, sigma[0], IA, IAbar)
                else:
                    p = len(sigma)
                    j = self.jacc(I, sigma[-1], IA, IAbar)
                for i in range(2, len(sigma)):
                    jprec = self.jacc(I, sigma[i - 1], IA, IAbar)
                    jsucc = self.jacc(I, sigma[i], IA, IAbar)
                    jcurr = self.jacc(sigma[i - 1], sigma[i], IA, IAbar)
                    if max(jprec, jsucc) > j and min(jprec, jsucc) >= jcurr:
                        p = i
                        j = max(jprec, jsucc)
                sigma.insert(p, I)
        return sigma

    def roleArea(self, r):
        return len(self.rowClusters[r]) * len(self.columnClusters[r])

    def importanceScore(self, item, IA):
        R = []
        for i, _, r in IA:
            if i == item:
                R.append(r)
        return sum([self.roleArea(r) for r in R])

    def jacc(self, i, j, IA, IAbar):
        ri, rj = [], []
        for item, items, role in IA:
            if item == i:
                ri.append(role)
            if item == j:
                rj.append(role)
        riUrj = []
        riNrj = []
        for r in ri:
            if r in rj:
                riNrj.append(r)
            riUrj.append(r)
        for r in rj:
            if r not in riUrj:
                riUrj.append(r)

        rU = 0
        rN = 0
        for u, items, r in IAbar:
            if r in riUrj:
                rU += len(items)
            if r in riNrj:
                rN += len(items)
        if rU == 0:
            return 0
        return rN / rU

    def blockClusterSimilarityScore(
        self, i, j, yBlocks, xBlocksCluster, yClustersBlock
    ):
        blocki = xBlocksCluster[i]
        blockj = xBlocksCluster[j]
        interCluster = blocki.intersection(blockj)
        unionCluster = blocki.union(blockj)
        interScore = 0
        unionScore = 0
        for c in unionCluster:
            for block in yClustersBlock[c]:
                unionScore += len(yBlocks[block])
        if unionScore == 0:
            return 0
        for c in interCluster:
            for block in yClustersBlock[c]:
                interScore += len(yBlocks[block])
        return interScore / unionScore

    def getGlobalBlockSimilarityScore(self):
        totalScore = 0
        for i in range(1, self.m):
            totalScore += self.blockElementSimilarityScore(
                self.rowBlocksPermutation[i - 1], self.rowBlocksPermutation[i], ROW
            )
        for i in range(1, self.n):
            totalScore += self.blockElementSimilarityScore(
                self.columnBlocksPermutation[i - 1],
                self.columnBlocksPermutation[i],
                COL,
            )
        return totalScore

    def getGlobalBlockClusterSimilarityScore(self):
        totalScore = 0
        for i in range(1, self.m):
            totalScore += self.blockClusterSimilarityScore(
                self.rowBlocksPermutation[i - 1],
                self.rowBlocksPermutation[i],
                self.columnBlocks,
                self.rowBlocksCluster,
                self.columnClustersBlock,
            )
        for i in range(1, self.n):
            totalScore += self.blockClusterSimilarityScore(
                self.columnBlocksPermutation[i - 1],
                self.columnBlocksPermutation[i],
                self.rowBlocks,
                self.columnBlocksCluster,
                self.rowClustersBlock,
            )
        return totalScore

    #######################################################################################################
    #######################        Adviser visualisation cost metric      #################################
    #######################################################################################################

    def totalProximityScore(self):
        return sum([self.ProximityScore(c) for c in range(self.nbClusters)])

    def ProximityScore(self, c):
        return self.roleSize(ROW, c) * self.roleSize(COL, c)

    def globalVisualizationCost(self):
        return sum(
            [self.roleVisualizationCost(role) for role in range(self.nbClusters)]
        )

    def roleVisualizationCost(self, role):
        piu = self.numberRoleFragments(ROW, role)
        pip = self.numberRoleFragments(COL, role)
        wu = self.roleSize(ROW, role)
        wp = self.roleSize(COL, role)
        return piu * pip * (wu * wp - self.clustersSize[role])

    def numberRoleFragments(self, mode, role):
        if mode == ROW:
            roleInBlock = sorted(
                [
                    list(self.rowBlocksPermutation).index(b)
                    for b in self.rowClustersBlock[role]
                ]
            )
        if mode == COL:
            roleInBlock = sorted(
                [
                    list(self.columnBlocksPermutation).index(b)
                    for b in self.columnClustersBlock[role]
                ]
            )
        nbFragments = 1
        for i in range(1, len(roleInBlock)):
            if roleInBlock[i] != roleInBlock[i - 1] + 1:
                nbFragments += 1
        return nbFragments

    def roleSize(self, mode, role):
        if mode == ROW:
            xBlocks = self.rowBlocks
            xBlocksPermutation = self.rowBlocksPermutation
            xClustersBlock = self.rowClustersBlock
        if mode == COL:
            xBlocks = self.columnBlocks
            xBlocksPermutation = self.columnBlocksPermutation
            xClustersBlock = self.columnClustersBlock

        blocksInCluster = xClustersBlock[role]
        indexOfBlocksInPermutation = sorted(
            [list(xBlocksPermutation).index(b) for b in blocksInCluster]
        )
        firstIndex = indexOfBlocksInPermutation[0]
        lastIndex = indexOfBlocksInPermutation[-1]
        w = 0
        for index in range(firstIndex, lastIndex + 1):
            currentBlock = xBlocksPermutation[index]
            w += len(xBlocks[currentBlock])
        return w

    ##################################################################################################
    ##############################	 Cluster Area Score	 #############################################
    ##################################################################################################

    def clusterAreaScore(self):
        totalScore = 0
        for i in range(self.nbClusters):
            totalScore += self.clusterScore(i)
        return totalScore  # maxTheoreticalScore - totalScore

    def clusterScore(self, i):
        score = 0
        compactRowCluster = self.compactCluster(
            self.rowClustersBlock[i], self.rowBlocksPermutation, self.rowBlocksSize
        )
        compactColumnCluster = self.compactCluster(
            self.columnClustersBlock[i],
            self.columnBlocksPermutation,
            self.columnBlocksSize,
        )
        for row in compactRowCluster:
            for col in compactColumnCluster:
                # print(f"[{i}] ({row},{col})")
                score += (row * col) ** 2
        return score

    def compactCluster(self, xBlocksCluster, xBlocksPermutation, xBlocksSize):
        xPermutedBlocksCluster = [
            list(xBlocksPermutation).index(blockCluster)
            for blockCluster in xBlocksCluster
        ]
        xPermutedBlocksCluster.sort()
        xPermutedBlocksSize = [xBlocksSize[i] for i in xBlocksPermutation]
        compactedCluster = [xPermutedBlocksSize[xPermutedBlocksCluster[0]]]
        for j in range(1, len(xPermutedBlocksCluster)):
            # if the blocks are contiguous in the permutation
            if xPermutedBlocksCluster[j] == xPermutedBlocksCluster[j - 1] + 1:
                compactedCluster[-1] += xPermutedBlocksSize[xPermutedBlocksCluster[j]]
            else:
                compactedCluster.append(xPermutedBlocksSize[xPermutedBlocksCluster[j]])
        return compactedCluster

    #######################################################################################################
    ##############################   	 Similarity score  	 #############################################
    #######################################################################################################

    def similarity(
        self, i1, i2, mode="3"
    ):  # returns the sBlocks() imilarity score between B_{pi(i1)} and B_{pi(i2)}
        if i1 == i2:
            if mode == "TSP":
                return np.inf
            return 0
        similarityValue = 0
        columnBlock1 = self.C[:, self.columnBlocksPermutation[i1]]
        columnBlock2 = self.C[:, self.columnBlocksPermutation[i2]]
        only_zeroes_1 = True
        only_zeroes_2 = True
        for j in range(self.m):
            if (columnBlock1[self.rowBlocksPermutation[j]]) != 0:
                only_zeroes_1 = False
            if (columnBlock2[self.rowBlocksPermutation[j]]) != 0:
                only_zeroes_2 = False
            if (
                columnBlock1[self.rowBlocksPermutation[j]]
                * columnBlock2[self.rowBlocksPermutation[j]]
            ):
                similarityValue += self.closeness(i1, i2, j, mode)
        if only_zeroes_2 or only_zeroes_1:
            return 0
        if mode == "TSP":
            if similarityValue == 0:
                return 0
            return 1 / similarityValue

        # print(similarityValue)
        return similarityValue

    # row similarity score, used for the disjoint row case
    def similarityRows(
        self, i1, i2, mode="3"
    ):  # returns the sBlocks() imilarity score between B_{pi(i1)} and B_{pi(i2)}
        if i1 == i2:
            if mode == "TSP":
                return np.inf
            return 0
        similarityValue = 0
        rowBlock1 = self.C[self.rowBlocksPermutation[i1], :]
        rowBlock2 = self.C[self.rowBlocksPermutation[i2], :]
        for j in range(self.n):
            if (
                rowBlock1[self.columnBlocksPermutation[j]]
                * rowBlock2[self.columnBlocksPermutation[j]]
            ):
                similarityValue += self.closenessRows(i1, i2, j)

        # print(similarityValue)
        if mode == "TSP":
            if similarityValue == 0:
                return 0
            return 1 / similarityValue
        return similarityValue

    " compute the weights of having two sepecific 1s aligned from B_{pi(i1)} and B_{pi(i2)} "

    def closenessRows(self, i1, i2, j):
        block1 = self.rowBlocks[self.rowBlocksPermutation[i1]]
        block2 = self.rowBlocks[self.rowBlocksPermutation[i2]]
        return (len(block1) + len(block2)) * len(
            self.columnBlocks[self.columnBlocksPermutation[j]]
        )

    " compute the weights of having two sepecific 1s aligned from B_{pi(i1)} and B_{pi(i2)} "

    def closeness(self, i1, i2, j, mode):
        block1 = self.columnBlocks[self.columnBlocksPermutation[i1]]
        block2 = self.columnBlocks[self.columnBlocksPermutation[i2]]
        if mode == 1:
            return (len(block1) * len(block2)) ** 5
        if mode == 2:
            return len(block1) + len(block2)
        if mode == 3 or mode == "TSP":
            return (len(block1) + len(block2)) * len(
                self.rowBlocks[self.rowBlocksPermutation[j]]
            )

    #######################################################################################################
    ##############################   	  Plot matrix   	 #############################################
    #######################################################################################################

    def plotOriginalMatrix(self, axs=None, permutationRow=None, permutationCol=None):
        mat = self.returnOriginalMatrix(permutationRow, permutationCol)
        for i in range(self.rows):
            for j in range(self.cols):
                mat[i, j] = 1 if mat[i, j] != 0 else 0
        for i in range(self.rows):
            for j in range(self.cols):
                if mat[i, j] != 0 and mat[i, j] != 1:
                    print(mat[i, j], i, j)
        rp = self.getFinalPermutation(ROW)
        cp = self.getFinalPermutation(COL)
        for i in range(len(self.rowClusters)):
            cr, cl = self.rowClusters[i], self.columnClusters[i]
            for i in cr:
                i = rp.index(i)
                for j in cl:
                    j = cp.index(j)
                    if mat[i, j] % 2 == 0:
                        mat[i, j] = 2
                    else:
                        mat[i, j] = 3
        for o in self.rowOrphanBlocks:
            orphans = o["orphans"]
            blueprint = o["blueprint"]
            indexes = [cp.index(index[0]) for index in np.argwhere(blueprint == True)]
            # indexes of the opposite side of the orphans
            for orphan in orphans:
                orphan = rp.index(orphan)
                for j in indexes:
                    mat[orphan, j] = 5 if mat[orphan, j] == 1 else 4
        for o in self.columnOrphanBlocks:
            orphans = o["orphans"]
            blueprint = o["blueprint"]
            indexes = [rp.index(index[0]) for index in np.argwhere(blueprint == True)]
            for orphan in orphans:
                orphan = cp.index(orphan)
                for i in indexes:
                    mat[i, orphan] = 5 if mat[i, orphan] == 1 else 4
        hexcolors = []
        if self.rowOrphanBlocks == [] and self.columnOrphanBlocks == []:
            hexcolors += [
                "#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c",
            ]
        else:
            hexcolors += [
                "#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c",
                "#fb9a99",
                "#e31a1c",
            ]

        colors = [mplcolors.to_rgb(color) for color in hexcolors]

        cmap = mplcolors.ListedColormap(colors)
        if axs == None:
            print(mat)
        else:
            axs.matshow(mat, cmap=cmap)

    def returnOriginalMatrix(self, permutationRow=None, permutationCol=None):
        if permutationRow == None:
            rowPermut = self.getFinalPermutation(ROW)
        else:
            rowPermut = permutationRow
        if permutationCol == None:
            columnPermut = self.getFinalPermutation(COL)
        else:
            columnPermut = permutationCol
        mat = self.originalMatrix[rowPermut, :]
        mat[:] = mat[:, columnPermut]
        return mat
