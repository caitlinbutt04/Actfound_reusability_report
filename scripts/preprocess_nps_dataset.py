# preprocesses the antibacterial NPs dataset
import pandas as pd

np_df = pd.read_csv('ethnobotany_and_the_role_of_plant_natural_products_in_antibiotic_drug_discovery_compounds.csv')

# remove trailing whitespace in all columns
for col in np_df.columns:
    if np_df[col].dtype == 'object':  # Check if the column contains strings
        np_df[col] = np_df[col].str.rstrip()

# keep only MIC values with ug/mL units
np_df = np_df[np_df['Test type (IC50 or MIC)'] == 'MIC']
np_df = np_df[np_df['unit'] == 'µg/mL']

# average MIC values for duplicate compound/organism pairs
np_df_grouped = (
    np_df.groupby(['SMILES', 'Organism'], as_index=False)["# value"]
      .mean()
)

# remove assays with less than 20 compounds
organism_counts = np_df_grouped['Organism'].value_counts()
valid_organisms = organism_counts[organism_counts >= 20].index
np_df_filtered = np_df_grouped[np_df_grouped['Organism'].isin(valid_organisms)]

# save results
np_df_filtered.to_csv('reproducibility_report_np_dataset_cleaned.csv', index = False)