# Standard libraries
import datetime
import logging
from os.path import join
from typing import List

# Libraries for graphs
import community
import networkx as nx

# Libraries for matrix computations
import numpy as np
import pandas as pd

# Import functions from nibna module
from nibna.network_node_importance import build_graph_from_edge_list, compute_pearson_correlation, compute_top_k_precision_recall, nibna, validate_top_k_coding_genes
from nibna.utils import configure_logging, find_coding_noncoding_drivers, plot_node_importance, plot_precision_recall_curves, NUM_miR

# DATA_DIR = "/home/mandar/Data/NCSU/CBNA/cbna-community-detection/Data/"
# OUT_DIR = "/home/mandar/Data/NCSU/CBNA/cbna-community-detection/Output/"

DATA_DIR = "/Users/mandar.chaudhary/Research Wednesday/NIBNA/Data/"
OUT_DIR = "/Users/mandar.chaudhary/Research Wednesday/NIBNA/Output/"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


if __name__ =="__main__":
    configure_logging(logger)
    
    start_time = datetime.datetime.now()
    top_k = [50, 100, 150, 200]
    top_k_validated_coding_genes = []
    
    logger.info("Reading cancer network edge list")
    cancer_network = pd.read_csv(join(DATA_DIR, "pVal_cancer_network.csv"))

    logger.info("Reading cancer gene expression data")
    cancer_df = pd.read_csv(join(DATA_DIR, "cancer_data.csv"))

    logger.info("Reading gold standard CGC data")
    gold_standard_cgc_df = pd.read_csv(join(DATA_DIR, "Census_allFri Sep 28 07_39_37 2018.tsv"), sep="\t")
    gold_standard_cgc = gold_standard_cgc_df["Gene Symbol"]

    logger.info("Reading mRNAs data")
    mRNAs_df = pd.read_csv(join(DATA_DIR, "mRNAs_data.csv"))
    
    protein_mutations = pd.read_csv(join(DATA_DIR, 'protein_affecting_mutations.csv'))
    mutation_subtypes = pd.read_csv(join(DATA_DIR, 'mutation_subtypes.csv'))
          
    logger.info("Computing pearson correlation coefficient")
    cancer_network = compute_pearson_correlation(cancer_network, cancer_df)
     
    logger.info("Building graph from cancer network edge list with pearson correlation as edge weights")
    G = build_graph_from_edge_list(cancer_network, is_weighted=True)
    
    num_top_k_validated_genes = []    
    top_k_validated_coding_genes = list()
    np.random.seed(10)

    node_importance_df = nibna(G)
    threshold_value = 1/len(G.nodes())
    node_threshold = (node_importance_df.importance > threshold_value).sum() # select nodes with importance score > 1/n
    logger.info(f'Node importance threshold: {threshold_value}, threshold: {node_threshold}')

    critical_nodes = node_importance_df.iloc[:node_threshold].copy()
    critical_nodes.to_csv(join(OUT_DIR, 'critical_nodes.csv'))
    plot_node_importance(node_importance_df.importance.to_numpy(), node_threshold, 'cancer', OUT_DIR)

    logger.info('Identifying coding drivers with mutations, coding drivers without mutations and non-coding drivers')
    cancer_drivers_dict = find_coding_noncoding_drivers(critical_nodes, protein_mutations,
                                                            cancer_df.columns[:NUM_miR])
    
    cancer_drivers_dict['coding_drivers_mutations'].to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'coding_candidate_drivers_mutations.csv'))
    cancer_drivers_dict['coding_drivers_no_mutations'].to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'coding_candidate_drivers_no_mutations.csv'))
    cancer_drivers_dict['noncoding_drivers'].to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'noncoding_candidate_drivers.csv'))    
    logger.info("Validating node importance scores with gold standard cgc")
    for k in top_k:
        top_k_coding_genes = validate_top_k_coding_genes(cancer_drivers_dict['coding_drivers_mutations'], mRNAs_df,
                                                         gold_standard_cgc, k)
        top_k_precision, top_k_recall = compute_top_k_precision_recall(cancer_drivers_dict['coding_drivers_mutations'],
                                                                       mRNAs_df, gold_standard_cgc)
        top_k_f1 = (2 * top_k_precision * top_k_recall)/(top_k_precision + top_k_recall + 1e-6)
        plot_precision_recall_curves(top_k_precision, top_k_recall, top_k_f1, OUT_DIR)
        performance_metrics_df = pd.DataFrame(zip(np.arange(top_k[-1]), top_k_precision, top_k_recall, top_k_f1), 
                                              columns=['k', 'precision', 'recall', 'f1'])
        performance_metrics_df.to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'performance_metrics_new.csv'))
        # Write results to csv file
        top_k_coding_genes.to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', f'top_k_{k}_validated_genes.csv'))

        num_top_k_validated_genes.append((top_k_coding_genes['In CGC']=='Yes').sum())
    
    num_top_k_validated_genes = np.array(num_top_k_validated_genes).reshape(1, -1)
    num_top_k_validated_genes = pd.DataFrame(num_top_k_validated_genes, columns=['top-50', 'top-100', 'top-150', 'top-200'])
    num_top_k_validated_genes.to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer',"top_k_validated_genes_weighted.csv"))

    logger.info("Number of top-K coding driver genes found")
    logger.info("{0}".format(num_top_k_validated_genes))
    logger.info("Experiment took {0} seconds".format((datetime.datetime.now() - start_time).seconds))