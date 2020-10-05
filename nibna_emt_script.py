# Standard libraries
import datetime
import logging
import os
from os.path import join
from typing import List

# Libraries for matrix computations
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

# Import functions from nibna module
from nibna.network_node_importance import build_graph_from_edge_list, compute_pearson_correlation, nibna
from nibna.utils import compute_mutation_count, configure_logging, find_coding_noncoding_drivers, DATA_DIR, OUT_DIR, NUM_miR

K = 100
    
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


if __name__ =="__main__":
    configure_logging(logger)
    
    # Read Mes data and network
    logger.info(f"Reading gene expression data for EMT analysis")
    mes_df = pd.read_csv(join(DATA_DIR, 'EMT', 'Mes', 'pVal_cancer_data.csv'))
    
    logger.info(f"Reading network edge list for EMT analysis")
    mes_network = pd.read_csv(join(DATA_DIR, 'EMT', 'Mes', 'pVal_cancer_network.csv'))
    
    logger.info("Computing pearson correlation coefficient")
    mes_network = compute_pearson_correlation(mes_network, mes_df)
    
    logger.info("Building graph from cancer network edge list with pearson correlation as edge weights")
    G = build_graph_from_edge_list(mes_network, is_weighted=True)
    
    np.random.seed(10)
    node_importance_df = nibna(G)
    
    threshold_value = 1/len(G.nodes())
    node_threshold = (node_importance_df.importance > threshold_value).sum() # select nodes with importance score > 1/n

    if not os.path.exists(join(OUT_DIR, 'EMT', 'Mes')):
        os.makedirs(join(OUT_DIR, 'EMT', 'Mes'))

    logger.info(f'Node importance threshold: {threshold_value}, threshold: {node_threshold}')
    critical_nodes = node_importance_df.iloc[:node_threshold].copy()
    critical_nodes["Type"] = "coding"
    critical_nodes.loc[critical_nodes.node.isin(mes_df.columns[:NUM_miR]), "Type"] = "non-coding"
    
    logger.info(f'Saving critical nodes')
    critical_nodes.to_csv(join(OUT_DIR, 'EMT', 'Mes', 'critical_nodes.csv'))

    logger.info('Loading gold standard data for validation')
    coding_gold_standard_df = pd.read_csv(join(DATA_DIR, 'EMT', 'Generic_EMT_signature.csv'))
    noncoding_gold_standard_mes_df = pd.read_csv(join(DATA_DIR, 'EMT', 'EMT_processed_miRNAs.csv'))
    
    coding_gold_standard_mes_df = coding_gold_standard_df.loc[coding_gold_standard_df.Epi=='Mes']
    
    validated_coding_candidate_drivers_mes = set(critical_nodes.node[:K]).intersection(set(coding_gold_standard_mes_df.KRT19))
    logger.info(f'Number of validated coding candidate drivers in top-{K} predicted drivers: {len(validated_coding_candidate_drivers_mes)}')
    
    validated_noncoding_candidate_drivers_mes = set(critical_nodes.node[:K]).intersection(noncoding_gold_standard_mes_df.V1)
    logger.info(f'Number of validated noncoding candidate drivers in top-{K} predicted drivers: {len(validated_noncoding_candidate_drivers_mes)}')
    
    top_20_coding_drivers = critical_nodes.loc[critical_nodes.Type=="coding"].iloc[:20, :]
    top_20_noncoding_drivers = critical_nodes.loc[critical_nodes.Type=="non-coding"].iloc[:20, :]
    
    top_20_coding_drivers.to_csv(join(OUT_DIR, 'EMT', 'Mes', 'top_20_coding_drivers.csv'))
    top_20_noncoding_drivers.to_csv(join(OUT_DIR, 'EMT', 'Mes', 'top_20_non_coding_drivers.csv'))
    
    # Calculate p_value
    # Coding
    # Mes
    n_A = K # top-100 estimated coding drivers
    n_B = len(coding_gold_standard_mes_df) # number of signatures
    n_C = 839 + 5168 # number of genes
    n_A_B = len(validated_coding_candidate_drivers_mes) # number of drivers validated
    p_value = 1 - hypergeom.cdf(n_A_B, n_C, n_B, n_A)
    logger.info("p_value for coding drivers of Mes: {0}".format(p_value))
    
    n_A = (critical_nodes.Type=="non-coding").sum() # number of estimated non-coding drivers
    n_B = len(noncoding_gold_standard_mes_df) # number of signatures
    n_C = NUM_miR # number of miRNA
    n_A_B = len(validated_noncoding_candidate_drivers_mes) # number of drivers validated
    p_value = 1 - hypergeom.cdf(n_A_B, n_C, n_B, n_A)
    logger.info("p_value for non-coding drivers of Mes: {0}".format(p_value))
