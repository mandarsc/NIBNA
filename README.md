## Network-based Node Importance Approach for Identifying Cancer Driver Genes

This repository contains code for identifying driver genes in the cancer network. In this work, a community detection algorithm is used to partition the cancer network into communities and a centrality-based metric is implemented to compute the importance of each node in the network [2]. Specifically, the well-known Louvain algorithm is used for detecting communities, and a centrality-based metric is used to compute node importance by measuring the distortion in the community structure upon the node's removal from the network. 

The following are the main steps involved in computing node importance,
1. Build a undirected weighted graph G using the edge list of the cancer network.
2. Partition the graph into communities using the Louvain algorithm.
3. Compute eigenvectors of the adjancency matrix of graph G.
4. Compute node importance of each node in graph G using a centrality-based metric.
5. Sort the coding genes in the cancer network in descending order of their importance score.

## Steps to run experiments
1. ### Identify critical nodes in the cancer network: Run the script `nibna_cancer_driver_script.py` to obtain the list of predicted coding drivers with mutations, coding drivers without mutations and non-coding drivers in the cancer network.
```
python3 nibna_cancer_driver_script.py
```
The script will load all the input data files and output results under the `NIBNA/Output/CancerDriver/Cancer/` directory. The list of files created by the script are as follows,
1. `critical_nodes.csv` contains list of all predicted cancer drivers.
2. cancer_node_importance.jpg contains a plot showing the distribution of node importance scores.
3. top_k_50_validated_genes.csv contains top-50 predicted coding cancer drivers. Similarly, the remaining file names with same name convention contain predicted cancer drivers for different values of threshold.
4. top_k_validated_genes_weighted.csv contains the number of predicted coding drivers validated using CGC.
5. coding_candidate_drivers_mutations.csv contains list of predicted coding drivers with mutations.
6. coding_candidate_drivers_no_mutations.csv contains list of predicted coding drivers without mutations.
7. noncoding_candidate_drivers.csv contains list of predicted non-coding drivers.
8. performance_metrics.csv contains precision, recall and f1-score of the predicted coding cancer drivers.

The results are saved in a csv file saved in `Output` directory where each row indicates the number of top-k coding genes found by this approach.

### References:
1. [Pham, Vu VH, Lin Liu, Cameron P. Bracken, Gregory J. Goodall, Qi Long, Jiuyong Li, and Thuc D. Le. "CBNA: A control theory based method for identifying coding and non-coding cancer drivers." PLoS Computational Biology 15, no. 12 (2019).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007538#sec009)
2. [Wang, Y., Di, Z., & Fan, Y. (2011). Identifying and characterizing nodes important to community structure using the spectrum of the graph. PloS one, 6(11).](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027418)
