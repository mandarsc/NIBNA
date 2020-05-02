## Community Detection-based Node Importance for Identifying Coding Genes

This notebook contains code for identifying important nodes in the cancer network. The cancer network was built in [1]. In this work, a community detection algorithm is used to partition the cancer network into communities and a centrality-based metric is implemented to compute the importance of each node in the network [2]. Specifically, the well-known Louvain algorithm is used for detecting communities, and a centrality-based metric is used to compute the importance of each node with respect to its community. 

The following steps are performed to compute the node importance and validate the important nodes,
1. Build an undirected graph G using the edge list of the cancer network.
2. Partition the graph into communities using the Louvain algorithm.
3. Compute eigenvectors of the adjancency matrix of graph G.
4. Compute node importance for each node in G using a centrality-based metric.
5. Sort the coding genes in the cancer network in descending order.
6. Validate the coding genes with gold standard CGC.

### References:
1. [Pham, Vu VH, Lin Liu, Cameron P. Bracken, Gregory J. Goodall, Qi Long, Jiuyong Li, and Thuc D. Le. "CBNA: A control theory based method for identifying coding and non-coding cancer drivers." PLoS Computational Biology 15, no. 12 (2019).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007538#sec009)
2. [Wang, Y., Di, Z., & Fan, Y. (2011). Identifying and characterizing nodes important to community structure using the spectrum of the graph. PloS one, 6(11).](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027418)
