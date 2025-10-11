# code used to analyze the assays and compounds present in the natural products dataset
# searches ChEMBL for matching assays in two ways:
#   1. searches ChEMBL for the references in the natural products dataset  
#   2. searches ChEMBL for growth inhibitory assays containing the compounds in the natural products dataset. 

import pandas as pd
from chembl_webresource_client.new_client import new_client

# load the dataset containing the antibacterial compound references for the NPs dataset
doi_df = pd.read_csv('ethnobotany_and_the_role_of_plant_natural_products_in_antibiotic_drug_discovery_references.csv')

# load the dataset containing the NPs antibacterial compounds
np_df = pd.read_csv('reproducibility_report_np_dataset_cleaned.csv')

def find_ref(doi_list):
    """
    Finds the references that have been deposited into ChEMBL based on 
    their DOI
    
    Parameters
    ---------
    doi_list: list
        list of DOIs
        
    Returns
    -------
    documents: DataFrame
        The documents found in ChEMBL
        
    Note
    ----
    This function requires chembl_webresource_client to be installed.
    """
    documents_api = new_client.document
    documents = documents_api.get(doi=doi_list)
    documents = pd.DataFrame.from_records(documents)
    
    return documents
    
def remove_matching_ref(ref_df, np_df):
    """
    Removes the compounds in the references that have been deposited into ChEMBL
    from the NPs dataset
    
    Parameters
    ---------
    ref_df: DataFrame
        DataFrame containing the DOIs of the references
    np_df: DataFrame
        DataFrame containing the antibacterial NPs
        
    Returns
    -------
    documents: DataFrame
        DataFrame containing the documents found in ChEMBL
    """
    ref_df.dropna(subset = ['DOI'], inplace=True)
    doi_id = ref_df['DOI'].tolist()
    documents = find_ref(doi_id)
    ref_in_chembl = pd.merge(documents, ref_df, left_on='doi', right_on='DOI')
    np_df_ref_removed = np_df[~np_df['Reference'].isin(ref_in_chembl['Reference Number*'])]
    
    return np_df_ref_removed
    
def find_chembl_id(np_df):
    """
    Finds the molecule chembl ids of the antibacterial NPs
    
    Parameters
    ---------
    np_df: DataFrame
        DataFrame containing the antibacterial NPs
        
    Returns
    -------
    ids_df: DataFrame
        DataFrame containing the molecule chembl ids of the NPs found in ChEMBL
        
    Note
    ----
    This function requires chembl_webresource_client to be installed.
    """
    similarity = new_client.similarity
    
    # remove duplicate SMILES
    smiles_df = np_df[['SMILES']]
    smiles_df.dropna(subset = ['SMILES'], inplace = True)
    smiles_df.drop_duplicates(subset=['SMILES'], inplace = True)
    smiles = smiles_df['SMILES'].tolist()
    
    ids_log = []
    for smile in smiles:
        try:
            res = similarity.filter(smiles = smile, similarity = 100).only(['molecule_chembl_id', 'similarity', 'molecule_structures'])
            for i in res:
                if i is None:
                    continue
                ids_log.append({
                    'query_smiles': smile,
                    'molecule_chembl_id': i['molecule_chembl_id'],
                    'similarity': i['similarity'],
                    'canonical_smiles': i['molecule_structures']['canonical_smiles']
                })
        except Exception as e:
            print(f"Error processing SMILES: {smile} - {e}")

    ids_df = pd.DataFrame(ids_log)
    
    return ids_df
    
def find_assays(ids_df):
    """
    Finds the assays associated with the molecule chembl ids of the antibacterial 
    NPs found in ChEMBL
    
    Parameters
    ---------
    ids_df: DataFrame
        DataFrame containing the molecule chembl ids of the NPs found in ChEMBL
        
    Returns
    -------
    activity_df: DataFrame
        DataFrame containing the activities of the NPs found in ChEMBL
        
    Note
    ----
    This function requires chembl_webresource_client to be installed.
    """
    activity = new_client.activity
    
    ids_list = ids_df['molecule_chembl_id'].tolist()
    activity_log = []
    for ids in ids_list:
        activities = activity.filter(molecule_chembl_id = ids).filter(standard_type = 'MIC', standard_units = 'ug.mL-1', standard_relation = '=')
        for i in activities:
            activity_log.append({'action_type': i['action_type'],
                                 'activity_comment': i['activity_comment'],
                                 'activity_id': i['activity_id'],
                                 'activity_properties': i['activity_properties'],
                                 'assay_chembl_id': i['assay_chembl_id'],
                                 'assay_description': i['assay_description'],
                                 'assay_type': i['assay_type'],
                                 'assay_variant_accession': i['assay_variant_accession'],
                                 'assay_variant_mutation': i['assay_variant_mutation'],
                                 'bao_endpoint': i['bao_endpoint'],
                                 'bao_format': i['bao_format'],
                                 'bao_label': i['bao_label'],
                                 'canonical_smiles': i['canonical_smiles'],
                                 'data_validity_comment': i['data_validity_comment'],
                                 'data_validity_description': i['data_validity_description'],
                                 'document_chembl_id': i['document_chembl_id'],
                                 'document_journal': i['document_journal'],
                                 'document_year': i['document_year'],
                                 'ligand_efficiency': i['ligand_efficiency'],
                                 'molecule_chembl_id': i['molecule_chembl_id'],
                                 'molecule_pref_name': i['molecule_pref_name'],
                                 'parent_molecule_chembl_id': i['parent_molecule_chembl_id'],
                                 'pchembl_value': i['pchembl_value'],
                                 'potential_duplicate': i['potential_duplicate'],
                                 'qudt_units': i['qudt_units'],
                                 'record_id': i['record_id'],
                                 'relation': i['relation'],
                                 'src_id': i['src_id'],
                                 'standard_flag': i['standard_flag'],
                                 'standard_relation': i['standard_relation'],
                                 'standard_text_value': i['standard_text_value'],
                                 'standard_type': i['standard_type'],
                                 'standard_units': i['standard_units'],
                                 'standard_upper_value': i['standard_upper_value'],
                                 'standard_value': i['standard_value'],
                                 'target_chembl_id': i['target_chembl_id'],
                                 'target_organism': i['target_organism'],
                                 'target_pref_name': i['target_pref_name'],
                                 'target_tax_id': i['target_tax_id'],
                                 'text_value': i['text_value'],
                                 'toid': i['toid'],
                                 'type': i['type'],
                                 'units': i['units'],
                                 'uo_units': i['uo_units'],
                                 'upper_value': i['upper_value'],
                                 'value': i['value']})
    
    activity_df = pd.DataFrame(activity_log)
    
    return activity_df
    
def remove_matching_assays(np_df):
    """
    Removes compound/target organism pairs from the NPs dataset that match 
    compound/target organism pairs found in growth inhibitory assays in ChEMBL
    
    Parameters
    ---------
    np_df: DataFrame
        DataFrame containing the antibacterial NPs 
        
    Returns
    -------
    np_df_activities_removed: DataFrame
        DataFrame containing the compound/target organism pairs not found in 
        ChEMBL
    """
    ids_df = find_chembl_id(np_df)
    activity_df = find_assays(ids_df)
    
    # keep activities from target organisms present in the NPs dataset
    activity_df = activity_df[activity_df['bao_label'] == 'organism-based format']
    activity_df = activity_df[activity_df['target_organism'].isin(np_df['Organism'])]
    activity_df = activity_df[['molecule_chembl_id', 'target_organism']]
    activity_df = activity_df.drop_duplicates(subset=['molecule_chembl_id', 'target_organism'])
    
    # add the molecule chembl ids to the NPs dataset
    ids_df = ids_df[['query_smiles', 'molecule_chembl_id']]
    np_df_with_ids = pd.merge(ids_df, np_df, left_on='query_smiles', right_on='SMILES')
    
    # keep the compound/target organism pairs not found in ChEMBL
    np_df_activities_removed = pd.merge(np_df_with_ids, activity_df, left_on=['molecule_chembl_id', 'Organism'], right_on=['molecule_chembl_id', 'target_organism'], how = 'left', indicator = True)
    np_df_activities_removed = np_df_activities_removed[np_df_activities_removed['_merge'] == 'left_only']
    np_df_activities_removed = np_df_activities_removed.drop(columns=['_merge'])
    np_df_activities_removed = np_df_activities_removed.drop_duplicates(subset=['SMILES', 'Organism'])

    return np_df_activities_removed

# remove the compounds in the references that have been deposited into ChEMBL
np_df_ref_removed = remove_matching_ref(doi_df, np_df)
np_df_ref_removed.to_csv('reproducibility_report_np_dataset_with_ref_removed.csv', index = False)
 
# remove compound/organism pairs that match pairs found in growth inhibitory assays in ChEMBL   
np_df_activities_removed = remove_matching_assays(np_df)
np_df_activities_removed.to_csv('reproducibility_report_np_dataset_with_assays_removed.csv', index = False)