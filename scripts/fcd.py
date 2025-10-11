# calculates the Frechnet ChemNet Distance between the training and finetuning datasets
# code was modified from the implementation provided by Preuer et al. at https://github.com/bioinf-jku/FCD 
import pandas as pd
from fcd import get_fcd, canonical_smiles

chembl_data = pd.read_csv('/data/walker_lab/caitlinbutt/actfound_data/chembl/chembl_processed_chembl32_repeat_assays_removed.csv')
bdb_data = pd.read_csv('/data/walker_lab/caitlinbutt/actfound_data/BDB/bdb_all_polymers_assays_removed.csv')
finetune_data = pd.read_csv('/home/buttc/ActFound/reproducibility_report_np_dataset_cleaned.csv')

# keep only the SMILES column
chembl_data = chembl_data[['SMILES']]
bdb_data = bdb_data[['SMILES']]
finetune_data = finetune_data[['SMILES']]

# select 10,000 random compounds from the training data
chembl_data = chembl_data.sample(n = 10000, random_state = 1)
bdb_data = bdb_data.sample(n = 10000, random_state = 1)

# label datasets for identification later
chembl_data['dataset'] = 'ChEMBL'
bdb_data['dataset'] = 'BindingDB'
finetune_data['dataset'] = 'NPs'

# get canonical smiles and filter invalid ones
can_chembl = [w for w in canonical_smiles(chembl_data['SMILES']) if w is not None]
can_bdb = [w for w in canonical_smiles(bdb_data['SMILES']) if w is not None]
can_finetune = [w for w in canonical_smiles(finetune_data['SMILES']) if w is not None]
    
# calculate FCD score
fcd_score = get_fcd(smiles1 = can_chembl, smiles2 = can_finetune)
print("FCD score between finetuning and ChEMBL data: ", fcd_score)

fcd_score = get_fcd(smiles1 = can_bdb, smiles2 = can_finetune)
print("FCD score between finetuning and BindingDB data: ", fcd_score)