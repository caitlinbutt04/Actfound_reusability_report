# calculates the scaffold novelty of the finetuning data
#code was modified from the implementation provided by Petrov and Bender at https://github.com/PangeAI/SCINS.git 
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

chembl_data = pd.read_csv('chembl_processed_chembl32.csv')
bdb_data = pd.read_csv('bdb_all_polymers.csv')
finetune_data = pd.read_csv('reproducibility_report_np_dataset_cleaned.csv')

# keep only the SMILES column
chembl_data = chembl_data[['SMILES']]
bdb_data = bdb_data[['SMILES']]
finetune_data = finetune_data[['SMILES', 'Organism']]

# remove invalid SMILES from finetune data 
finetune_data["mol"] = [Chem.MolFromSmiles(smiles) for smiles in finetune_data["SMILES"]]
finetune_data = finetune_data.dropna(subset=["mol"])

# label datasets for identification later
chembl_data['dataset'] = 'ChEMBL'
bdb_data['dataset'] = 'BindingDB'
finetune_data['dataset'] = 'NPs'

def get_bemis_murcko_scaffolds(df):
    """
    Generates the Bemis Murcko scaffold
    
    Parameters
    ---------
    df: DataFrame
        DataFrame containing SMILES strings in column 'SMILES'
        
    Returns
    -------
    df: DataFrame
        DatFrame containing the scaffolds in the column 'murcko_scaffolds'
    """
    scaffolds = []
    for smiles in df['SMILES']:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(smiles)
        if len(scaffold) == 0:
            scaffold = smiles
            scaffolds.append(scaffold)
        else:
            scaffolds.append(scaffold)

    df['murcko_scaffold'] = scaffolds
    return df

def scaled_shannon_entropy(counts):
    """
    Computes the scaled Shannon entropy (SSE)
    
    Parameters
    ---------
    counts: DataFrame
        DataFrame of scaffold counts
        
    Returns
    -------
    SSE: int
        scaled Shannon entropy
    """
    # total compounds
    P = counts.sum()
    probs = counts / P
    
    # shannon entropy (base 2)
    SE = -(probs * np.log2(probs)).sum()
    
    # scaled shannon entropy
    SSE = SE / np.log2(len(counts))
    
    return SSE

def get_sse(df, most_pop):
    """
    Computes the scaled Shannon entropy (SSE)
    
    Parameters
    ---------
    df: DataFrame
        DataFrame containing the scaffolds in the column 'murcko_scaffold'
    
    most_pop: list
        list of the amount of most populated scaffolds you would like to 
        compute the SSE on
        
    Returns
    -------
    sse_values: DataFrame
        DataFrame of SSE values for each most populated scaffold
    """
    scaffold_counts = df['murcko_scaffold'].value_counts()
    
    sse_values = []
    for i in most_pop:
        top_scaffolds = scaffold_counts.head(i)
        sse_value = scaled_shannon_entropy(top_scaffolds)
        print(f"SSE of top {i} scaffolds: {sse_value:.4f}")
        sse_values.append({
            'dataset': df['dataset'][0],
            'number of scaffolds': i,
            'SSE': sse_value
        })
    sse_values = pd.DataFrame(sse_values)
    
    return sse_values

chembl_data = get_bemis_murcko_scaffolds(chembl_data)
bdb_data = get_bemis_murcko_scaffolds(bdb_data)
finetune_data = get_bemis_murcko_scaffolds(finetune_data)

# save scaffolds
chembl_data.to_csv('/data/walker_lab/caitlinbutt/actfound_data/chembl/chembl_processed_chembl32_scaffolds.csv', index = False)
bdb_data.to_csv('/data/walker_lab/caitlinbutt/actfound_data/BDB/all_bdb_polymer_scaffolds.csv', index = False)
finetune_data.to_csv('/home/buttc/ActFound/reproducibility_report_np_dataset_cleaned_scaffolds.csv', index = False)

# compute scaffold novelty between finetuning and ChEMBL data
unique_scaffolds = finetune_data.drop_duplicates(subset = ['murcko_scaffold'])
novel_scaffolds = unique_scaffolds[~unique_scaffolds['murcko_scaffold'].isin(chembl_data['murcko_scaffold'])]
perc_novel = len(novel_scaffolds) / len(unique_scaffolds) * 100
print(f"Percentage of novel scaffolds between finetuning and ChEMBL data: {perc_novel}")

# compute scaffold novelty between finetuning and BindingDB data
unique_scaffolds = finetune_data.drop_duplicates(subset = ['murcko_scaffold'])
novel_scaffolds = unique_scaffolds[~unique_scaffolds['murcko_scaffold'].isin(bdb_data['murcko_scaffold'])]
perc_novel = len(novel_scaffolds) / len(unique_scaffolds) * 100
print(f"Percentage of novel scaffolds between finetuning and BindingDB data: {perc_novel}")

# compute the SSE values for each dataset
chembl_sse_values = get_sse(chembl_data, most_pop = [5, 10, 20, 50, 100])
bdb_sse_values = get_sse(bdb_data, most_pop = [5, 10, 20, 50, 100])
finetune_sse_values = get_sse(finetune_data, most_pop = [5, 10, 20, 50, 100])

# concat and save datasets
all_sse_values = pd.concat([chembl_sse_values, bdb_sse_values, finetune_sse_values])
all_sse_values.to_csv('sse_values.csv', index = False)