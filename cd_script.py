# Standard libraries
import argparse
import datetime
import logging
from os.path import join
import pickle
from typing import List

# Libraries for graphs
import community
import networkx as nx

# Libraries for matrix computations
import numpy as np
import pandas as pd

# Libraries for sparse matrices and eigenvectors
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.stats import pearsonr

DATA_DIR = "/home/mandar/Data/NCSU/CBNA/cbna-community-detection/Data/"
OUT_DIR = "/home/mandar/Data/NCSU/CBNA/cbna-community-detection/Output/"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def configure_logging():
    # create console handler and set level to info
    ch = logging.StreamHandler() 
    ch.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)


def compute_node_importance(n_nodes, n_communities, eig_vectors):
    """
    Computes node importance using a centrality-based metric which uses the eigenvectors of the adjacency matrix.
    Args:
        n_nodes: Number of nodes in the graph.
        n_communities: Number of communities estimated.
        eig_vectors: np.ndarray of size n x c where n is number of nodes and c is the number of communities.
        The eigenvectors are computed from the adjacency matrix and are stored as columns.
    Returns: 
        np.array: A numpy array containing the node importance score normalized by the
        number of communities.
    """
    p_k = [0] * n_nodes
    I = [0] * n_nodes
    for k in range(n_nodes):
        for i in range(n_communities):
            p_k[k] += eig_vectors[k, i]**2/np.sum(eig_vectors[:, i]**2)
        I[k] = p_k[k]/n_communities
    return np.array(I)


def computePearsonCorrelation(cancer_network: pd.DataFrame, cancer_data: pd.DataFrame) -> pd.DataFrame:
    """
    This function computes pearson correlation coefficient for every pair of adjacent nodes in the cancer
    network. The correlation coefficient is stored in the `correlation` column and its absolute value is
    stored in the `weight` column. 
    Args:
        cancer_network: Pandas dataframe containing edge list of the cancer network.
        cancer_data: Pandas dataframe containing the expression level of nodes in the cancer network.
    Returns:
        pd.DataFrame containing 3 more columns added to cancer_network dataframe. 
    """
    correlation_coeffs = []
    correlation_pvals = []
    x = cancer_network.cause.values.tolist()
    y = cancer_network.effect.values.tolist()

    for x_i, y_i in zip(x, y):
        corr, pval = pearsonr(cancer_data.loc[:, x_i], cancer_data.loc[:, y_i])
        correlation_coeffs.append(corr)
        correlation_pvals.append(pval)
    
    cancer_network["correlation"] = pd.Series(correlation_coeffs)
    cancer_network["pval"] = pd.Series(correlation_pvals)
    cancer_network["weight"] = cancer_network["correlation"].abs()
    return cancer_network


def buildGraphFromEdgeList(cancer_network: pd.DataFrame, is_weighted=False) -> nx.Graph:
    """
    This function builds a networkx graph object from the edge list. A boolean option is provided to build
    a weighted undirected graph.
    Args:
        cancer_network: Pandas dataframe containing edge list of the cancer network.
        is_weighted: boolean specifying whether a weighted graph should be created.
    Returns:
        nx.Graph object containing undirected edges between pairs of nodes specified in cancer_network dataframe.
    """
    if is_weighted:
        try:
            # Build graph by reading the edge list of the cancer network
            G = nx.from_pandas_edgelist(cancer_network, source="cause", target="effect", edge_attr=['weight'])
        except:
            raise ValueError("Dataframe does not contain a column containing edge weights")
    else:
        # Build graph by reading the edge list of the cancer network
        G = nx.from_pandas_edgelist(cancer_network, source="cause", target="effect")
    logger.info("Number of nodes: {0}, and number of edges: {1}".format(len(G.nodes), len(G.edges)))
    logger.info("Graph weights: {0}".format(G.size(weight='weight')))
    return G


def performCommunityDetection(G: nx.Graph) -> int:
    """
    This function performs community detection on graph G using the Louvain algorithm.
    Args:
        G: nx.Graph object representing an undirected graph
    Returns:
        int specifying the number of communities detected in the graph.
    """
    # Run community detection algorithm on graph G
    partition = community.best_partition(G)
    n_communities = len(set(partition.values()))
    
    logger.info(f"Number of communities in G: {n_communities}")
    
    return n_communities


def computeEigenvectors(G: nx.Graph, n_communities: int) -> np.ndarray:
    """
    This function computes the eigenvectors of the adjacency matrix. The number of eigenvectors returned 
    is specified by the n_communities integer. 
    """
    # Create adjacency matrix of graph G
    A = nx.adjacency_matrix(G)
    
    # Compute eigenvalues and eigenvectors of the adjacency matrix A
    _, adj_eig_vectors = eigs(A.toarray().astype(float), k=n_communities, which='LR')
    
    return adj_eig_vectors
          
                  
def buildNodeImportanceDataFrame(G: nx.Graph, I: np.array) -> pd.DataFrame:
    """
    This function builds a dataframe containing two columns. First column contains the names of the nodes in
    graph G and second column contains their importance score. The nodes are sorted in descending order
    of their importance score.
    Args:
        G: networkx graph object
        I: numpy array containing importance score.
    Returns:
        pd.DataFrame containing node names and their importance score.
    """
    # Get list of all nodes in the cancer network
    nodes = list(G.nodes)
    
    # Sort the nodes by their importance score in descending order
    important_nodes = [nodes[i] for i in np.argsort(I)[::-1]]
    node_importance = np.sort(I)[::-1]
    
    # Build a dataframe with two columns, one containing node name and second containing importance score
    node_importance_df = pd.concat([pd.Series(important_nodes), pd.Series(node_importance)], axis=1)
    node_importance_df.columns = ['node', 'importance']
    return node_importance_df
          

def validateTopKCodingGenes(node_importance_df: pd.DataFrame, mRNAs_data: pd.DataFrame, gold_standard_cgc: pd.Series) -> List[str]:
    """
    This function performs the following steps,
        1. Find the coding genes in the cancer network by matching the node names with those in mRNAs data
        2. Select the top-K coding genes found in the cancer network and validate them against the gold standard cgc.
    Args:
        node_importance_df: Pandas dataframe containig nodes sorted in descending order of their importance
        scores.
        mRNAs_data: Pandas dataframe with columns containing names of coding genes.
        gold_standard_cgc: Pandas series containing gold standard coding genes
    """
    # Step 1: Find coding genes in the cancer network
    # Perform intersection of the nodes in the cancer network and the mRNAs
    coding_genes = node_importance_df.loc[node_importance_df.node.isin(mRNAs_data.columns), ]
    # Select top-K coding genes and save the results to csv files
    top_k = [50, 100, 150, 200]
    top_k_validated_coding_genes = []
    for K in top_k:
        top_k_coding_genes = coding_genes.iloc[:K, ].copy()
        # Step 2: Validate top-k coding genes against the gold standard cgc
        coding_genes_gold_standard = top_k_coding_genes.loc[top_k_coding_genes.node.isin(gold_standard_cgc)]
#         logger.info("Coding genes in top-{0}: {1}".format(K, coding_genes_gold_standard.shape[0]))
        
        # Store top-k validated coding genes in a list
        top_k_validated_coding_genes.append(coding_genes_gold_standard.node.values)
    return top_k_validated_coding_genes
          

if __name__ =="__main__":
    configure_logging()
    
    logger.info("Initialize arg parser")
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-n', '--num_iterations', type=int, required=True, default=10)
    arg_parser.add_argument('-weighted', '--weighted_graph', type=int, required=True, default=0)
    args = vars(arg_parser.parse_args())
    n_iter = args['num_iterations']
    weighted_graph = args['weighted_graph']
    try:
        isinstance(n_iter, int)
    except:
        raise TypeError("num_itrations argument must be int type")
    
    start_time = datetime.datetime.now()
    
    logger.info("Reading cancer network edge list")
    cancer_network = pd.read_csv(join(DATA_DIR, "pVal_cancer_network.csv"))

    logger.info("Reading cancer genes expression data")
    cancer_data = pd.read_csv(join(DATA_DIR, "cancer_data.csv"))

    logger.info("Reading gold standard CGC data")
    gold_standard_cgc_df = pd.read_csv(join(DATA_DIR, "Census_allFri Sep 28 07_39_37 2018.tsv"), sep="\t")
    gold_standard_cgc = gold_standard_cgc_df["Gene Symbol"]

    logger.info("Reading mRNAs genes data")
    mRNAs_df = pd.read_csv(join(DATA_DIR, "mRNAs_data.csv"))
          
    logger.info("Computing pearson correlation coefficient")
    cancer_network = computePearsonCorrelation(cancer_network, cancer_data)
     
    logger.info("Building graph from cancer network edge list")
    if weighted_graph:
        logger.info("Building graph with pearson correlation as edge weights")
        G = buildGraphFromEdgeList(cancer_network, is_weighted=True)
    else:
        G = buildGraphFromEdgeList(cancer_network)
    
    logger.info("Computing node importance {0} times".format(n_iter))
    validated_coding_genes = dict()
    num_top_k_validated_genes = np.zeros((n_iter, 4))
    avg_num_top_k_validated_genes = np.zeros(4)
    std_num_top_k_validated_genes = np.zeros(4)
    
    np.random.seed(10)
    for i in range(n_iter):    
        logger.info("Iteration: {0}".format(i))
        logger.info("Partitioning graph into communities")
        n_communities = performCommunityDetection(G)

        logger.info("Computing eigenvectors from the adjacency matrix of graph G")
        adj_eigen_vectors = computeEigenvectors(G, n_communities)

        logger.info("Computing importance of nodes in the cancer network")
        I = compute_node_importance(G.number_of_nodes(), n_communities, np.real(adj_eigen_vectors))
        node_importance_df = buildNodeImportanceDataFrame(G, I)

        logger.info("Validating node importance scores with gold standard cgc")
        validated_coding_genes[i] = validateTopKCodingGenes(node_importance_df, mRNAs_df, gold_standard_cgc)
       
        for j in range(len(validated_coding_genes[i])):
            num_top_k_validated_genes[i, j] = len(validated_coding_genes[i][j])
            avg_num_top_k_validated_genes[j] += len(validated_coding_genes[i][j])
    
    num_top_k_validated_genes = pd.DataFrame(num_top_k_validated_genes, columns=['top-50', 'top-100', 'top-150', 'top-200'])
    if weighted_graph:
        num_top_k_validated_genes.to_csv(join(OUT_DIR, "top_k_validated_genes_weighted.csv"))
    else:
        num_top_k_validated_genes.to_csv(join(OUT_DIR, "top_k_validated_genes_unweighted.csv"))
    logger.info("Average top-K coding driver genes found")
    logger.info("{0}".format(np.mean(num_top_k_validated_genes, axis=0).ravel()))
    logger.info("Std top-K coding driver genes found")
    logger.info("{0}".format(np.std(num_top_k_validated_genes, axis=0).ravel()))
    logger.info("Experiment took {0} seconds".format((datetime.datetime.now() - start_time).seconds))