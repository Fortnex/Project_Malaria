from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd

# Load references
references = []
with open('mol2mol.smi', 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue  # Skip blank lines and comments
        smi = line.split()[0]
        mol = Chem.MolFromSmiles(smi)
        if mol:
            references.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2))

# Load molecules and their properties (use your scored results)
df = pd.read_csv('scoring_results.csv')  # <-- update to your file name

def max_tanimoto(smiles, references):
    mol = Chem.MolFromSmiles(smiles)
    if not mol or not references:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    return max(DataStructs.TanimotoSimilarity(fp, ref_fp) for ref_fp in references)

# If your scored CSV already has a 'Tanimoto' column, you can skip this step
df['Tanimoto'] = df['SMILES'].apply(lambda smi: max_tanimoto(smi, references))
filtered = df[df['Tanimoto'] >= 0.6]

# Drop the 'Input_SMILES' column if it exists
if 'Input_SMILES' in filtered.columns:
    filtered = filtered.drop(columns=['Input_SMILES'])

# ---- Property Filters ----
filtered = filtered[
    (filtered['Molecular weight'] <= 550) &
    (filtered['QED'] >= 0.75) &
    (filtered['SA score'] <= 4) &
    (filtered['SlogP (RDKit)'] >= -1) &
    (filtered['SlogP (RDKit)'] <= 5)
]

print(f"Filtered {len(filtered)} molecules with all property thresholds out of {len(df)} total.")

# Save result
filtered.to_csv('filtered_scoring_results.csv', index=False)
