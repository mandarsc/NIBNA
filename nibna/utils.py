import logging
from os.path import join
from typing import Dict, List

# Libraries for matrix computations
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as matlib
import pandas as pd

NUM_miR = 1719
NUM_mR = 20101
NUM_TF = 839


def compute_importance_score_threshold(importance_scores):  
    n_points = len(importance_scores)
    all_coord = np.vstack((range(n_points), importance_scores)).T
    np.array([range(n_points), importance_scores])
    
    first_point = all_coord[0]
    line_vec = all_coord[-1] - all_coord[0]
    line_vec_norm = line_vec / np.sqrt(np.sum(line_vec**2))
    vec_from_first = all_coord - first_point
    scalar_product = np.sum(vec_from_first * matlib.repmat(line_vec_norm, n_points, 1), axis=1)
    vec_from_first_parallel = np.outer(scalar_product, line_vec_norm)
    vec_to_line = vec_from_first - vec_from_first_parallel
    dist_to_line = np.sqrt(np.sum(vec_to_line ** 2, axis=1))
    idx_of_best_point = np.argmax(dist_to_line)    
    return idx_of_best_point


def configure_logging(logger):
    # create console handler and set level to info
    logger_handler = logging.StreamHandler() 
    logger_handler.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    logger_handler.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(logger_handler)
    

def find_coding_noncoding_drivers(critical_nodes: pd.DataFrame, protein_mutations: pd.DataFrame, miR_genes: List[str]) -> Dict[str, pd.DataFrame]:
    """
    """
    cancer_drivers_dict = dict()

    candidate_cancer_drivers = critical_nodes.merge(protein_mutations, how='left', left_on='node', right_on='symbol').copy()
    candidate_cancer_drivers["Type"] = "coding"
    candidate_cancer_drivers.loc[candidate_cancer_drivers.node.isin(miR_genes), "Type"] = "non-coding"


    coding_candidate_cancer_drivers = candidate_cancer_drivers.loc[candidate_cancer_drivers.Type=='coding'].copy()
    coding_candidate_cancer_drivers = coding_candidate_cancer_drivers.sort_values(by='importance', ascending=False)

    coding_candidate_cancer_drivers_mutations = coding_candidate_cancer_drivers.loc[~pd.isna(coding_candidate_cancer_drivers.weight)].copy()
    coding_candidate_cancer_drivers_no_mutations = coding_candidate_cancer_drivers.loc[pd.isna(coding_candidate_cancer_drivers.weight)].copy()

    # Candidate non-coding cancer drivers
    noncoding_candidate_cancer_drivers = candidate_cancer_drivers.loc[candidate_cancer_drivers.Type=='non-coding'].copy()
    
    cancer_drivers_dict['coding_drivers_mutations'] = coding_candidate_cancer_drivers_mutations
    cancer_drivers_dict['coding_drivers_no_mutations'] = coding_candidate_cancer_drivers_no_mutations
    cancer_drivers_dict['noncoding_drivers'] = noncoding_candidate_cancer_drivers
    return cancer_drivers_dict

    
def compute_mutation_count(coding_drivers_mutations: pd.DataFrame, mutation_subtypes:pd.DataFrame, cancer_subtype: str):
    """
    """
    subtype_mutation = []

    for gene in coding_drivers_mutations.node.values:
        temp = mutation_subtypes.loc[mutation_subtypes.symbol == gene,].copy()
        if all(temp.Call.isna()):
            subtype_mutation.append('nan')
            continue
        cancer_subtype_counts = temp.Call.value_counts()
        subtype_mutation.append(cancer_subtype_counts.index[np.argmax(cancer_subtype_counts)])

    coding_drivers_mutations['Subtype'] = subtype_mutation
    coding_drivers_mutations = coding_drivers_mutations.loc[coding_drivers_mutations.Subtype == cancer_subtype,]
    return coding_drivers_mutations
    
def plot_node_importance(node_importance: np.array, threshold: int, network_name: str, out_dir: str):
    """
    """
    plt.scatter(np.arange(len(node_importance)), node_importance, marker='o', linestyle='-', color='blue')
    plt.xlabel('Genes')
    plt.ylabel('Node Importance score')
    plt.axvline(threshold, color='red', linestyle='--', label=f'threshold={threshold}')
    plt.legend(loc='upper right')
    plt.title(f'Node importance in {network_name} network')
    plt.savefig(join(out_dir, f'{network_name}_node_importance.jpg'))
    plt.close()
    
def plot_precision_recall_curves(top_k_precision: np.array, top_k_recall: np.array, top_k_f1: np.array, out_dir: str):
    """
    Plot precision, recall and f1 scores. The plots are stored in the out_dir.
    Args:
        top_k_precision: precision scores for different top-k thresholds
        top_k_recall: recall scores for different top-k thresholds
        top_k_f1: f1 scores for different top-k thresholds
    """
    plt.plot(top_k_precision, color='blue', marker='+', linestyle='-', linewidth=1)
    plt.xlabel('Top N genes')
    plt.ylabel('Precision according to CGC')
    plt.title('Precision Comparison')
    plt.savefig(join(out_dir, 'precision.jpg'))
    plt.close()
    plt.plot(top_k_recall, color='blue', marker='+', linestyle='-', linewidth=1)
    plt.xlabel('Top N genes')
    plt.ylabel('Recall according to CGC')
    plt.title('Recall Comparison')
    plt.savefig(join(out_dir, 'recall.jpg'))
    plt.close()
    plt.plot(top_k_f1, color='blue', marker='+', linestyle='-', linewidth=1)
    plt.xlabel('Top N genes')
    plt.ylabel('F1 Score according to CGC')
    plt.title('F1 Score Comparison')
    plt.savefig(join(out_dir, 'f1_score.jpg'))
    plt.close()
