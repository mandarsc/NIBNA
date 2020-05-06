## Community Detection-based Node Importance for Identifying Coding Driver Genes

This notebook contains code for identifying important nodes in the cancer network. The cancer network was built in [1]. In this work, a community detection algorithm is used to partition the cancer network into communities and a centrality-based metric is implemented to compute the importance of each node in the network [2]. Specifically, the well-known Louvain algorithm is used for detecting communities, and a centrality-based metric is used to compute the importance of each node with respect to its community. 

The following steps are performed to compute the node importance and validate the important nodes,
1. Build an undirected graph G using the edge list of the cancer network.
2. Partition the graph into communities using the Louvain algorithm.
3. Compute eigenvectors of the adjancency matrix of graph G.
4. Compute node importance of each node in graph G using a centrality-based metric.
5. Sort the coding genes in the cancer network in descending order of their importance score.
6. Validate the coding genes with gold standard CGC.

## Steps to run experiments
In order to run experiments, there are two command line arguments that must be specified. The first argument specifies the number of times community detection must be performed and node importance to be computed. This step has been added since community detection is non-deterministic and therefore performing iterations of this step will provide some variance of the results. The second argument indicates whether an unweighted or weighted graph should be used for detecting communities.

The following command runs community detection and node importance 10 times on a weighted graph.

`python3 cd_script.py -n 10 -weighted 1`

To run experiments on an unweighted graph, you can specify the `weighted` argument to 0.

`python3 cd_script.py -n 10 -weighted 0`

The results are saved in a csv file saved in `Output` directory where each row indicates the number of top-k coding genes found by this approach.

### References:
1. [Pham, Vu VH, Lin Liu, Cameron P. Bracken, Gregory J. Goodall, Qi Long, Jiuyong Li, and Thuc D. Le. "CBNA: A control theory based method for identifying coding and non-coding cancer drivers." PLoS Computational Biology 15, no. 12 (2019).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007538#sec009)
2. [Wang, Y., Di, Z., & Fan, Y. (2011). Identifying and characterizing nodes important to community structure using the spectrum of the graph. PloS one, 6(11).](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027418)
