import logging
from os.path import join

# Libraries for matrix computations
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as matlib


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

