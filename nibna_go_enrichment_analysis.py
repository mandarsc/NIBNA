from collections import defaultdict
import logging
from os.path import join

import pandas as pd

from nibna.utils import add_genes_as_columns, get_top_20_drivers, get_top_drivers_enriched_terms, get_unique_genes, plot_heatmap, process_go_data, DATA_DIR, OUT_DIR


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

OUT_DIR = "/Users/mandar.chaudhary/Research Wednesday/NIBNA/Output/"


if __name__ == "__main__":
    
    logger.info('Reading GO biological processes and GO molecular functions enrichment profiles')
    go_biological_df = pd.read_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'GO_Biological_Process_2018_table.txt'), sep='\t')
    go_molecular_df = pd.read_csv(join(OUT_DIR, 'CancerDriver', 'Cancer', 'GO_Molecular_Function_2018_table.txt'), sep='\t')
    
    biological_10_enriched_terms = process_go_data(go_biological_df)
    molecular_10_enriched_terms = process_go_data(go_molecular_df)

    top_20_drivers_bio = get_top_20_drivers(biological_10_enriched_terms)
    top_20_drivers_mol = get_top_20_drivers(molecular_10_enriched_terms)
    
    unqiue_genes_bio = get_unique_genes(biological_10_enriched_terms)
    unqiue_genes_mol = get_unique_genes(molecular_10_enriched_terms)
    
    add_genes_as_columns(unqiue_genes_bio, biological_10_enriched_terms)
    add_genes_as_columns(unqiue_genes_mol, molecular_10_enriched_terms)
    
    go_bio_terms = get_top_drivers_enriched_terms(biological_10_enriched_terms, top_20_drivers_bio)
    go_mol_terms = get_top_drivers_enriched_terms(molecular_10_enriched_terms, top_20_drivers_mol)
    
    biological_10_enriched_terms[['enriched_terms', 'Term']].to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer',
                                                                    'biological_enriched_terms.csv'))

    molecular_10_enriched_terms[['enriched_terms', 'Term']].to_csv(join(OUT_DIR, 'CancerDriver', 'Cancer',
                                                                    'molecular_enriched_terms.csv'))
    
    plot_heatmap(go_bio_terms, out_dir = join(OUT_DIR, 'CancerDriver', 'Cancer'))
    plot_heatmap(go_mol_terms, out_dir = join(OUT_DIR, 'CancerDriver', 'Cancer'), go_process = "GO Molecular Function")