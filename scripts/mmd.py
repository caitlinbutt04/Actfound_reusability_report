# calculates the Maximum Mean Discrepancy between the training and finetuning datasets
# code was modified from the implementation provided by Klarner et al. at https://github.com/leojklarner/Q-SAVI
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from scipy.spatial.distance import squareform

# load the data
chembl_data = pd.read_csv('chembl_processed_chembl32.csv')
bdb_data = pd.read_csv('all_bdb_polymer.csv')
finetune_data = pd.read_csv('reproducibility_report_np_dataset_cleaned.csv')

# keep only the SMILES column
chembl_data = chembl_data[['SMILES']]
bdb_data = bdb_data[['SMILES']]
finetune_data = finetune_data[['SMILES', 'Organism']]

# label datasets for identification later
chembl_data['dataset'] = 'ChEMBL'
bdb_data['dataset'] = 'BindingDB'
finetune_data['dataset'] = 'NPs'

# select 10,000 random compounds from the training data
chembl_data = chembl_data.sample(n = 10000, random_state = 1)
bdb_data = bdb_data.sample(n = 10000, random_state = 1)

def generate_fps(df):
    """
    Generates Morgan Fingerprints
    
    Parameters
    ---------
    df: DataFrame
        DataFrame containing SMILES strings in column 'SMILES'
        
    Returns
    -------
    df: DataFrame
        DatFrame containing the fingerprints in the column 'morgan'
    """
    df["mol"] = [Chem.MolFromSmiles(smiles) for smiles in df["SMILES"]]
    num_before = len(df)
    df = df.dropna(subset=["mol"])
    num_after = len(df)
    delta = num_before - num_after
    if delta > 0:
        print(f"Removed {delta} bad records")
    fpgen = rdFingerprintGenerator.GetMorganGenerator()
    df["morgan"] = [fpgen.GetFingerprint(mol) for mol in df["mol"]]
    return df
    
def generate_sim_matrix(fps):
    """
    Generates pair-wise distances and similarities
    
    Parameters
    ---------
    fps: DataFrame
        DatFrame containing the fingerprints in the column 'morgan'
        
    Returns
    -------
    dists: list
        list containing the similarity matrix
    """
    all_fps = fps['morgan'].values

    # generate pair-wise distances and similarities
    dists = []
    for i in range(len(all_fps)):
        dists.extend(
            DataStructs.BulkTanimotoSimilarity(
                all_fps[i],
                all_fps[(i+1):],
                returnDistance=True
            )
        )
    dists = squareform(dists)
    return dists
    
def distributional_similarity(train_size, sim_mat):
    """
    Implements an unbiased estimator of the maximum mean discrepancy statistic 
    for a given training index and kernel similarity matrix.
    
    Parameters
    ---------
    train_size: int
        length of the training dataset
    
    sim_mat: list
        list containing the similarity matrix
        
    Returns
    -------
    np.sqrt(max(mmd_squared, 0))
        the MMD statistic
    """
    n = sim_mat.shape[0]
    train_idx = np.arange(train_size)
    test_idx = np.arange(train_size, n)

    # train-train mean similarity
    train_mean = (sim_mat[np.ix_(train_idx, train_idx)].sum() - train_size) / (train_size * (train_size - 1))
    # test-test mean similarity
    test_mean = (sim_mat[np.ix_(test_idx, test_idx)].sum() - len(test_idx)) / (len(test_idx) * (len(test_idx) - 1))
    # train-test mean similarity
    train_test_mean = sim_mat[np.ix_(train_idx, test_idx)].mean()

    mmd_squared = train_mean + test_mean - 2 * train_test_mean

    return np.sqrt(max(mmd_squared, 0))

# generate fps
chembl_data = generate_fps(chembl_data)
bdb_data = generate_fps(bdb_data)
finetune_data = generate_fps(finetune_data)
    
# compute MMD between ChEMBL data and finetuning data
all_fps = pd.concat([chembl_data, finetune_data])
precomputed_kernels = generate_sim_matrix(all_fps)
train_size = len(chembl_data)
mmd_stat = distributional_similarity(train_size, 1 - precomputed_kernels)
print(f"MMD statistic between ChEMBL and finetuning data: {mmd_stat}")

# compute MMD between BindingDB data and finetuning data
all_fps = pd.concat([bdb_data, finetune_data])
precomputed_kernels = generate_sim_matrix(all_fps)
train_size = len(bdb_data)
mmd_stat = distributional_similarity(train_size, 1 - precomputed_kernels)
print(f"MMD statistic between BindingDB and finetuning data: {mmd_stat}")