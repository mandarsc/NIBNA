import logging
from os.path import join
from typing import Dict, List

# Libraries for matrix computations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NUM_miR = 1719
NUM_mR = 20101
NUM_TF = 839


def configure_logging(logger):
    """
    This function intializes a logger to record logging statements.
    """
    # create console handler and set level to info
    logger_handler = logging.StreamHandler() 
    logger_handler.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    logger_handler.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(logger_handler)
    

def process_go_data(go_df: pd.DataFrame) -> pd.DataFrame:
    go_enriched_terms_df = go_df.iloc[:10].copy()
    go_enriched_terms_df['genes_split'] = go_enriched_terms_df.Genes.apply(lambda x: x.split(';'))
    go_enriched_terms_df['enriched_terms'] = go_enriched_terms_df.Term.str.extract(r'\(([^)]+)')
    return go_enriched_terms_df


def get_top_20_drivers(enriched_terms_df: pd.DataFrame) -> List[str]:
    gene_count_dict = defaultdict(int)
    for gene_list in enriched_terms_df.genes_split:
        for gene in gene_list:
            gene_count_dict[gene] += 1

    gene_count_df = pd.DataFrame.from_dict(gene_count_dict, orient='index')
    gene_count_df = gene_count_df.reset_index()
    gene_count_df.columns = ['gene', 'count']
    gene_count_df = gene_count_df.sort_values(by='count', ascending=False)
    top_20_drivers = gene_count_df.gene.head(n=20).tolist()
    return top_20_drivers


def get_unique_genes(enriched_terms_df: pd.DataFrame):
    unique_genes_per_term = enriched_terms_df.genes_split.apply(lambda x: set(x))
    unique_genes = set()
    for genes in unique_genes_per_term:
        unique_genes = unique_genes.union(genes)
    return unique_genes


def add_genes_as_columns(genes, df):
    for gene in genes:
        df[gene] = 0
        gene_in_encriched_term = df.genes_split.apply(lambda x: gene in x)
        df.loc[gene_in_encriched_term, gene] = 1

        
def get_top_drivers_enriched_terms(df, top_20_drivers):
    go_terms = df[['enriched_terms'] + top_20_drivers].T.copy()
    go_terms.columns = go_terms.iloc[0]
    go_terms = go_terms.drop('enriched_terms', axis=0)
    go_terms = go_terms.astype('int')
    return go_terms
    
    
def find_coding_noncoding_drivers(critical_nodes: pd.DataFrame, protein_mutations: pd.DataFrame, miR_genes: List[str]) -> Dict[str, pd.DataFrame]:
    """
    This function merges the critical nodes identified using NIBNA with proteins and divides them into coding and non-coding drivers.
    Args:
        critical_nodes (pd.DataFrame): Pandas dataframe containing critical nodes.
        protein_mutations (pd.DataFrame): Pandas dataframe containing proteints with mutations.
        miR_genes (List[str]): List of miRNA genes.
    Returns:
        Dict[str, pd.DataFrame]: Dictionary containing coding drivers with mutations, coding drivers without mutations and non-coding drivers.
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
    This function computes the frequency of count of the predicted coding drivers with mutations and selects the gene with the highest frequency.
    Args:
        coding_drivers_mutations(pd.DatFrame): Pandas dataframe containing predicted coding drivers with mutations.
        mutation_subtype (pd.DataFrame): Pandas dataframe containing mutation subtypes.
        cancer_subtype (str): String specifying the cancer subtype.
    Returns:
        Pandas dataframe containing coding drivers with mutations having the highest frequency count for the given cancer subtype.
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
    Plot node importance scores of genes and save it in specified directory.
    Args:
        node_importance (np.array): Numpy array containing node importance scores in decrasing order.
        threshold (int): Integer specifying last node with importance score greater than the threshold value.
        network_name (str): String specifying the condition-specific network for which node importance is computed.
        out_dir (str): String specifying output directory to save the plot.
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


def plot_heatmap(go_terms: pd.DataFrame, out_dir: str, go_process: str = "GO Biological Process"):
    rcParams['figure.figsize'] = 10, 7
    sns.heatmap(go_terms, cmap="Blues", linewidths=1, linecolor='black', cbar=False, vmin=0, vmax=1)
    plt.xlabel('Enriched Terms', fontsize=18)
    plt.ylabel('Predicted Cancer Drivers', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xticks(rotation=-60)
    plt.title(go_process, fontsize=20)
    plt.tight_layout()
    plt.savefig(join(out_dir, go_process + ' Heatmap.jpg'), bbox_inches = "tight")
    plt.close()