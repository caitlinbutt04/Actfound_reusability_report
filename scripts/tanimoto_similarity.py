import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

chembl_data = pd.read_csv('chembl_processed_chembl32.csv')
bdb_data = pd.read_csv('bdb_all_polymers.csv')
finetune_data = pd.read_csv('reproducibility_report_np_dataset_cleaned.csv')

# keep only the SMILES column
chembl_data = chembl_data[['SMILES']]
bdb_data = bdb_data[['SMILES']]
finetune_data = finetune_data[['SMILES', 'Organism']]

def prepare_dataset(df):
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


def find_most_similar(df_test: pd.DataFrame,df_train: pd.DataFrame) -> pd.DataFrame:
    """
    Generates the maximum Tanimoto similarity between each fingerprint in the
    df_test and fingerprints in df_train

    Parameters
    ---------
    df_test: DataFrame
        DataFrame containing the fine-tuning data
    
    df_train: DataFrame
        DataFrame containing the training data
        
    Returns
    -------
    df: DataFrame
        DatFrame containing the tanimoto similaritites in the column 'most_similar'
    """
    train_fps = df_train["morgan"].values
    sim_list = [
        max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) if train_fps.size > 0 else np.nan
        for fp in df_test["morgan"]]
    df_test = df_test.copy()
    df_test["most_similar"] = sim_list
    return df_test
    
chembl_data = prepare_dataset(chembl_data)
bdb_data = prepare_dataset(bdb_data)
finetune_data = prepare_dataset(finetune_data)

# label datasets
chembl_data['dataset'] = 'ChEMBL'
bdb_data['dataset'] = 'BindingDB'
finetune_data['dataset'] = 'NPs'

# evaluate similarity between the fine-tuning and ChEBML training data
finetune_chembl_data = find_most_similar(finetune_data, chembl_data)
finetune_chembl_data.to_csv('max_tan_sim_of_finetune_in_chembl_train.csv', index = False)

# evaluate similarity between the fine-tuning and BindingDB training data
finetune_bdb_data = find_most_similar(finetune_data, bdb_data)
finetune_bdb_data.to_csv('max_tan_sim_of_finetune_in_bdb_train.csv', index = False)