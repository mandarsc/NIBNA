# Standard libraries
import datetime
import logging
from os.path import join
import pickle
from typing import List

# Libraries for matrix computations
import numpy as np
import pandas as pd

# Import functions from nibna module
from nibna.network_node_importance import build_graph_from_edge_list, compute_pearson_correlation, nibna
from nibna.utils import compute_mutation_count, configure_logging, find_coding_noncoding_drivers, NUM_miR

DATA_DIR = "/Users/mandar.chaudhary/Research Wednesday/NIBNA/Data/"
OUT_DIR = "/Users/mandar.chaudhary/Research Wednesday/NIBNA/Output/"
SUBTYPE_DATA_DIR = "/Users/mandar.chaudhary/Research Wednesday/CancerDriver/Data/Output/CancerSubtype"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


if __name__ =="__main__":
    configure_logging(logger)
    
    cancer_subtypes = ["Basal", "Her2", "LumA", "LumB", "Normal"]
    cancer_network_filename = 'pVal_cancer_network.csv'
    
    protein_mutations = pd.read_csv(join(DATA_DIR, 'protein_affecting_mutations.csv'))
    mutation_subtypes = pd.read_csv(join(DATA_DIR, 'mutation_subtypes.csv'))

    subtype_coding_drivers_no_mutations_dict = dict()
    subtype_noncoding_drivers_dict = dict()
    
    subtype_coding_drivers_no_mutations = set()
    subtype_noncoding_drivers = set()
    common_coding_drivers_no_mutations = set()
    common_noncoding_drivers = set()

    for subtype in cancer_subtypes:
        logger.info(f'Cancer subtype: {subtype}')
        # Read Mes data and network
        filename = subtype.lower() + '_cancer_data.csv'
        
        logger.info(f"Reading {subtype} cancer network edge list")
        cancer_subtype_df = pd.read_csv(join(SUBTYPE_DATA_DIR, subtype, filename))
        
        logger.info(f"Reading {subtype} cancer gene expression data")
        cancer_subtype_network = pd.read_csv(join(SUBTYPE_DATA_DIR, subtype, cancer_network_filename))
        logger.info('Dataframe size: {0}'.format(cancer_subtype_df.shape))
        
        logger.info("Computing pearson correlation coefficient")
        subtype_network = compute_pearson_correlation(cancer_subtype_network, cancer_subtype_df)

        logger.info("Building graph from cancer network edge list with pearson correlation as edge weights")
        G = build_graph_from_edge_list(cancer_subtype_network, is_weighted=True)

        np.random.seed(10)
        node_importance_df = nibna(G)

        threshold_value = 1/len(G.nodes())
        node_threshold = (node_importance_df.importance > threshold_value).sum() # select nodes with importance score > 1/n
        
        logger.info(f'Node importance threshold: {threshold_value}, threshold: {node_threshold}')
        critical_nodes = node_importance_df.iloc[:node_threshold].copy()

        logger.info(f'Saving critical nodes')
        critical_nodes.to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'critical_nodes.csv'))
        
        logger.info('Identifying coding drivers with mutations, coding drivers without mutations and non-coding drivers')
        cancer_drivers_dict = find_coding_noncoding_drivers(critical_nodes, protein_mutations,
                                                            cancer_subtype_df.columns[:NUM_miR])

        coding_candidate_cancer_drivers_mutations = compute_mutation_count(cancer_drivers_dict['coding_drivers_mutations'],
                                                                          mutation_subtypes, subtype)
        logger.info(f"Coding candidate drivers with mutations: {len(cancer_drivers_dict['coding_drivers_mutations'])}")
        logger.info(f"Coding candidate drivers without mutations: {len(cancer_drivers_dict['coding_drivers_no_mutations'])}")
        logger.info(f"Non-coding candidate drivers without mutations: {len(cancer_drivers_dict['noncoding_drivers'])}")

        coding_candidate_cancer_drivers_mutations.to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'coding_candidate_drivers_mutations.csv'))
        cancer_drivers_dict['coding_drivers_no_mutations'].to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'coding_candidate_drivers_no_mutations.csv'))
        cancer_drivers_dict['noncoding_drivers'].to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'noncoding_candidate_drivers.csv'))
                    
        subtype_coding_drivers_no_mutations_dict[subtype] = cancer_drivers_dict['coding_drivers_no_mutations']
        subtype_noncoding_drivers_dict[subtype] = cancer_drivers_dict['noncoding_drivers']
                    
        # Identify coding drivers without mutations that are common to atleast two cancer subtypes
        common_coding_drivers_no_mutations = common_coding_drivers_no_mutations.union(subtype_coding_drivers_no_mutations.intersection(cancer_drivers_dict['coding_drivers_no_mutations'].node.unique()))      
        subtype_coding_drivers_no_mutations = subtype_coding_drivers_no_mutations.union(cancer_drivers_dict['coding_drivers_no_mutations'].node.unique())
        subtype_coding_drivers_no_mutations = subtype_coding_drivers_no_mutations.difference(common_coding_drivers_no_mutations)

        # Identify noncoding drivers that are common to atleast two cancer subtypes
        common_noncoding_drivers = common_noncoding_drivers.union(subtype_noncoding_drivers.intersection(cancer_drivers_dict['noncoding_drivers'].node.unique()))
        subtype_noncoding_drivers = subtype_noncoding_drivers.union(cancer_drivers_dict['noncoding_drivers'].node.unique())
        subtype_noncoding_drivers = subtype_noncoding_drivers.difference(common_noncoding_drivers)


    # Filter coding drivers without mutations and noncoding drivers that are common to atleast two cancer subtypes
    for subtype in cancer_subtypes:
        coding_drivers_no_mutations = subtype_coding_drivers_no_mutations_dict[subtype]
        noncoding_drivers = subtype_noncoding_drivers_dict[subtype]
        
        coding_drivers_no_mutations = coding_drivers_no_mutations.loc[~coding_drivers_no_mutations.node.isin(common_coding_drivers_no_mutations)]
        noncoding_drivers = noncoding_drivers.loc[~noncoding_drivers.node.isin(common_noncoding_drivers)]
        coding_drivers_no_mutations.to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'coding_candidate_drivers_no_mutations_no_overlap.csv'))
        noncoding_drivers.to_csv(join(OUT_DIR, 'CancerSubtype', subtype, 'noncoding_candidate_drivers_no_overlap.csv'))
