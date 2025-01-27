import pandas as pd
from scipy.sparse import csr_matrix
from anndata import AnnData
import numpy as np



def preprocess_single_cell_datasets(data_paths, experiment_label="WTA"):
    """
    Preprocess single-cell datasets: loads, transposes, makes into sparse matrix 

    Parameters:
    - data_paths
    - experiment labels for different datasets

    Returns:
    - AnnData: Combined AnnData object with all datasets processed
    """

    adata_list = []
    
    for i, path in enumerate(data_paths):
        # Load dataset and skip unnecessary rows pandas needed for the skipping
        df = pd.read_csv(path, skiprows=7, index_col=0)
        
        # Transpose
        df = df.T
        
        # sparse matrix
        sparse_data = csr_matrix(df.values)
        
        # Create AnnData object
        adata = AnnData(X=sparse_data, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
        
        # Add labels
        adata.obs['sample'] = f"{experiment_label}{i+1}"
        

        adata_list.append(adata)
    
    # Combine all datasets into one AnnData object
    adata_combined = adata_list[0].concatenate(*adata_list[1:], batch_key="sample")
    
    return adata_combined
