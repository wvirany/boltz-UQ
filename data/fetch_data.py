"""Fetch Runs-N-Poses dataset from Polaris and extract sequences and similarity metrics."""

import pandas as pd
import polaris as po

# Load the dataset from the Polaris
dataset = po.load_dataset("plinder-org/runs-n-poses-dataset")

# Pull all metadata rows into a dataframe for easy analysis
records = []
for i in range(len(dataset.rows)):
    if i % 10 == 0: print(i)
    dp = dataset[i]
    meta = dp['plinder_metadata']
    records.append({
        'group_key':          meta['group_key'],
        'system_id':          meta['system_id'],
        'entry_pdb_id':       meta['entry_pdb_id'],
        'ligand_ccd_code':    meta['ligand_ccd_code'],
        'ligand_smiles':      meta['ligand_smiles'],
        'num_protein_chains': meta['num_protein_chains'],
        'num_ligand_chains':  meta['num_ligand_chains'],
        'ligand_is_proper':   meta['ligand_is_proper'],
        'similarity_to_train':      dp['similarity_to_train'],          # sucos_shape_pocket_qcov
        'tanimoto':                 meta['tanimoto'],                   # 2D Morgan vs closest train ligand
        'protein_seqsim_max':       meta['protein_seqsim_max'],         # sequence identity vs closest train protein
    })

df = pd.DataFrame(records)
df.to_csv("metadata.csv", index=False)

print(f"Total rows:             {len(df)}")
print(f"Unique PDB entries:     {df['entry_pdb_id'].nunique()}")
print(f"Unique system_ids:      {df['system_id'].nunique()}")
print(f"Unique ligands (CCD):   {df['ligand_ccd_code'].nunique()}")
print(f"Unique SMILES:          {df['ligand_smiles'].nunique()}")
print()
print("Ligands per system_id:")
print(df.groupby('system_id')['ligand_ccd_code'].count().describe())
print()
print("Similarity distributions:")
print(df[['similarity_to_train', 'tanimoto', 'protein_seqsim_max']].describe())