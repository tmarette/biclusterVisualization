import argparse

from biclusteringVisualization.plotCreator import visualizeClustering


def build_parser():
    parser = argparse.ArgumentParser(
        description="Cluster and visualize a matrix, or visualize clusters from files."
    )
    parser.add_argument("--input", required=True, help="Input matrix file")
    parser.add_argument(
        "--output", required=True, help="Output plot file, e.g. plots/data/result.pdf"
    )
    parser.add_argument(
        "--clustering-algorithm",
        choices=("PCV", "basso"),
        default="PCV",
        help="Clustering algorithm (default: PCV)",
    )
    parser.add_argument(
        "--nb-clusters", type=int, default=5, help="Number of clusters (default: 5)"
    )
    parser.add_argument(
        "--ordering-method",
        choices=(
            "TSPheuristic",
            "ADVISER",
            "greedyProximityScore",
            "greedyClusterArea",
            "greedyUninter",
            "greedyDemerit",
            "random",
        ),
        default="TSPheuristic",
        help="Block ordering algorithm (default: TSPheuristic)",
    )
    parser.add_argument(
        "--row-clusters",
        help="File containing row clusters; must be paired with --column-clusters",
    )
    parser.add_argument(
        "--column-clusters",
        help="File containing column clusters; must be paired with --row-clusters",
    )
    parser.add_argument(
        "--heuristic-time",
        type=int,
        default=20,
        help="Time in seconds for TSPheuristic (default: 20)",
    )
    parser.add_argument(
        "--plot-unordered-matrix",
        action="store_true",
        help="Also save the matrix without reordering",
    )
    parser.add_argument(
        "--show-objective-function-values",
        action="store_true",
        help="Save objective function values as a CSV alongside the plot",
    )
    parser.add_argument(
        "--print-permutations",
        action="store_true",
        help="Print final row and column permutations",
    )
    parser.add_argument(
        "--no-custom-color-schemes",
        action="store_false",
        dest="custom_color_schemes",
        help="Use the default matrix colormap",
    )
    parser.set_defaults(custom_color_schemes=True)
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if (args.row_clusters is None) != (args.column_clusters is None):
        parser.error("--row-clusters and --column-clusters must be provided together")

    visualizeClustering(
        inputFile=args.input,
        outputFile=args.output,
        clusteringAlgorithm=args.clustering_algorithm,
        nbClusters=args.nb_clusters,
        orderingMethod=args.ordering_method,
        rowClusters=args.row_clusters or "",
        columnClusters=args.column_clusters or "",
        heuristicTime=args.heuristic_time,
        plotUnorderedMatrix=args.plot_unordered_matrix,
        showObjectiveFunctionValues=args.show_objective_function_values,
        printPermutations=args.print_permutations,
        customColorSchemes=args.custom_color_schemes,
    )


if __name__ == "__main__":
    main()
