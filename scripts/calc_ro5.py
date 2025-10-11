# calculates Lipinksi's Rule of 5 properties
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms
import seaborn as sns

chembl_data = pd.read_csv('chembl_processed_chembl32.csv')
bdb_data = pd.read_csv('bdb_all_polymers.csv')
finetune_data = pd.read_csv('reproducibility_report_np_dataset_cleaned.csv')

# keep only the SMILES column
chembl_data = chembl_data[['SMILES']]
bdb_data = bdb_data[['SMILES']]
finetune_data = finetune_data[['SMILES']]

def calc_Ro5(df):
    """
    Calculates the Lipinski's Rule of 5 properties of SMILES string in a pandas DataFrame.

    Parameters
    ---------
        df: DataFrame
            DataFrame containing the SMILES strings.

    Returns
    -------
        ro5_df: DataFrame
            DataFrame with the calculated properties.
    """
    ro5_prop = []
    for smile in df['SMILES']:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            ro5_prop.append({
                'SMILES': smile,
                'MolWt': MolWt(mol),
                'LogP': MolLogP(mol),
                'TPSA': TPSA(mol),
                'HBD': NumHDonors(mol),
                'HBA': NumHAcceptors(mol),
                'num_heavy_atoms': CalcNumHeavyAtoms(mol)
            })
    ro5_df = pd.DataFrame(ro5_prop)

    return ro5_df
    
finetune_ro5_data = calc_Ro5(finetune_data)
chembl_ro5_data = calc_Ro5(chembl_data)
bdb_ro5_data = calc_Ro5(bdb_data)

# save the data
finetune_ro5_data.to_csv('reproducibility_report_np_dataset_cleaned_ro5.csv', index = False)
chembl_ro5_data.to_csv('chembl_processed_chembl32_ro5.csv', index = False)
bdb_ro5_data.to_csv('all_bdb_polymer_ro5.csv', index = False)