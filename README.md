# BRIC

BRIC is a novel biclustering method for the detection of the repertoire of active gene regulatory signals(GRSs) within each single cell, based on which, we annotate the type and/or physiological state of each cell. 

BRIC consists of two major steps: (i) detecting all gene co-regulation modules (GCM) each corresponding to a GRS, and (ii) inferring the active GRS of each cell, and further its cell type or physiological state.

The program of BRIC consists of two part, one is biclustering and the other one is cell type prediction. For biclustering, C++ code is provided. The input in this part is expression matrix and the output is .blocks file, in which each block represent a bicluster.For cell type prediction part, two R scripts are provided. Users should first run the constructGraph.R to generate weighted graph. The input here is the .blocks files output in the biclustering step, and the output is weighted graphs. Then users could run MCL_clustering.R to group cells based on the weighted graph. Note that MCL is just one clustering choice,users may use other clustering methods.
